#include "fdtd.h"

#include <assert.h>
#include <execinfo.h>
#include <math.h>
#include <mpi.h>
#include <signal.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "init.h"
#include "sending.h"
#include "utilities.h"

const char* source_type_keywords[] = {[SINE] = "sine", [AUDIO] = "audio"};

const char* output_type_keywords[] = {[CUTX] = "cut_x",
                                      [CUTY] = "cut_y",
                                      [CUTZ] = "cut_z",
                                      [ALL] = "all",
                                      [POINT] = "point"};

const char* output_source_keywords[] = {[PRESSURE] = "pressure",
                                        [VELOCITYX] = "velocity_x",
                                        [VELOCITYY] = "velocity_y",
                                        [VELOCITYZ] = "velocity_z"};

void print_source(source_t* source) {
    printf(" Source infos:\n\n");

    if (source->type == AUDIO) {
        double duration = (double)source->numsamples / source->sampling;

        printf("          type: audio data file\n");
        printf("      sampling: %d Hz\n", source->sampling);
        printf("      duration: %g\n", duration);

    } else {
        printf("          type: sine wave\n");
        printf("     frequency: %g Hz\n", source->data[0]);
    }

    printf("    position x: %g\n", source->posx);
    printf("    position y: %g\n", source->posy);
    printf("    position z: %g\n\n", source->posy);
}

void print_output(output_t* output) {
    switch (output->source) {
        case PRESSURE:
            printf("      pressure: ");
            break;
        case VELOCITYX:
            printf("    velocity X: ");
            break;
        case VELOCITYY:
            printf("    velocity Y: ");
            break;
        case VELOCITYZ:
            printf("    velocity Z: ");
            break;

        default:
            break;
    }

    switch (output->type) {
        case ALL:
            printf("complete dump");
            break;
        case CUTX:
            printf("cut along the x axis at %g", output->posx);
            break;
        case CUTY:
            printf("cut along the y axis at %g", output->posy);
            break;
        case CUTZ:
            printf("cut along the z axis at %g", output->posz);
            break;
        case POINT:
            printf("single point at %g %g %g", output->posx, output->posy,
                   output->posz);
            break;

        default:
            break;
    }

    printf(" to file %s\n", output->filename);
}

/******************************************************************************
 * Output file functions                                                      *
 ******************************************************************************/

int write_output(output_t* output, data_t* data, int step, double time) {
    if (output == NULL || data == NULL) {
        DEBUG_PRINT("NULL pointer passed as argument");
        return 1;
    }

    output_type_t type = output->type;

    if (type == ALL) {
        return write_data(output->fp, data, step, time);
    }

    int m, n, p;
    closest_index(&data->grid, output->posx, output->posy, output->posz, &m, &n,
                  &p);

    int startm = (type == CUTX || type == POINT) ? m : 0;
    int startn = (type == CUTY || type == POINT) ? n : 0;
    int startp = (type == CUTZ || type == POINT) ? p : 0;

    int endm = (type == CUTX || type == POINT) ? m + 1 : NUMNODESX(data);
    int endn = (type == CUTY || type == POINT) ? n + 1 : NUMNODESY(data);
    int endp = (type == CUTZ || type == POINT) ? p + 1 : NUMNODESZ(data);

    data_t* tmpdata = allocate_data(&output->grid);

    for (m = startm; m < endm; m++) {
        for (n = startn; n < endn; n++) {
            for (p = startp; p < endp; p++) {
                int tmpm = m - startm;
                int tmpn = n - startn;
                int tmpp = p - startp;

                SETVALUE(tmpdata, tmpm, tmpn, tmpp, GETVALUE(data, m, n, p));
            }
        }
    }

    int writeok = (write_data(output->fp, tmpdata, step, time) == 0);

    free(tmpdata->vals);
    free(tmpdata);

    if (writeok == 0) {
        DEBUG_PRINT("Failed to write output data");
        return 1;
    }

    return 0;
}

int open_outputfile(output_t* output, grid_t* simgrid) {
    if (output == NULL || simgrid == NULL) {
        DEBUG_PRINT("Invalid NULL pointer in argment");
        return 1;
    }

    grid_t grid;

    output_type_t type = output->type;

    grid.numnodesx = (type == POINT || type == CUTX) ? 1 : simgrid->numnodesx;
    grid.numnodesy = (type == POINT || type == CUTY) ? 1 : simgrid->numnodesy;
    grid.numnodesz = (type == POINT || type == CUTZ) ? 1 : simgrid->numnodesz;

    grid.xmin = (type == POINT || type == CUTX) ? output->posx : simgrid->xmin;
    grid.xmax = (type == POINT || type == CUTX) ? output->posx : simgrid->xmax;

    grid.ymin = (type == POINT || type == CUTY) ? output->posy : simgrid->ymin;
    grid.ymax = (type == POINT || type == CUTY) ? output->posy : simgrid->ymax;

    grid.zmin = (type == POINT || type == CUTZ) ? output->posz : simgrid->zmin;
    grid.zmax = (type == POINT || type == CUTZ) ? output->posz : simgrid->zmax;

    FILE* fp;
    if ((fp = create_datafile(grid, output->filename)) == NULL) {
        DEBUG_PRINTF("Failed to open output file: '%s'", output->filename);
        return 1;
    }

    output->grid = grid;
    output->fp = fp;

    return 0;
}

/******************************************************************************
 * Parameter file functions                                                   *
 ******************************************************************************/

int read_audiosource(char* filename, source_t* source) {
    FILE* fp;
    if ((fp = fopen(filename, "rb")) == NULL) {
        DEBUG_PRINTF("Could not open source file '%s'", filename);
        return 1;
    }

    fseek(fp, 0, SEEK_END);
    size_t filesize = ftell(fp);
    rewind(fp);

    int numsamples = (filesize - sizeof(int)) / sizeof(double);

    int sampling;
    if (fread(&sampling, sizeof(int), 1, fp) != 1) {
        DEBUG_PRINT("Failed to read source data");
        fclose(fp);
        return 1;
    }

    double* data;
    if ((data = malloc(numsamples * sizeof(double))) == NULL) {
        DEBUG_PRINT("Failed to allocate memory for source data");
        return 1;
    }

    int readok = (fread(data, sizeof(double), numsamples, fp) == (long unsigned int)numsamples);

    fclose(fp);

    if (readok == 0) {
        DEBUG_PRINT("Failed to read source data");
        return 1;
    }

    source->data = data;
    source->numsamples = numsamples;
    source->sampling = sampling;

    return 0;
}

int read_outputparam(FILE* fp, output_t* output, int coords[]) {
    if (fp == NULL || output == NULL) {
        DEBUG_PRINT("NULL passed as argement");
        return 1;
    }

    char typekeyword[BUFSZ_SMALL];
    char sourcekeyword[BUFSZ_SMALL];
    char filename[BUFSZ_LARGE];
    char formatted_filename[BUFSZ_HUMONGOUS];

    double posxyz[3] = {0.0, 0.0, 0.0};

    if (fscanf(fp, BUFFMT_SMALL, typekeyword) != 1 ||
        fscanf(fp, BUFFMT_SMALL, sourcekeyword) != 1 ||
        fscanf(fp, BUFFMT_LARGE, filename) != 1) {
        DEBUG_PRINT("Failed to read an output parameter");
        return 1;
    }

    output_type_t type = CUTX;
    while (type < OUTPUT_TYPE_END &&
           strcmp(output_type_keywords[type], typekeyword) != 0) {
        type++;
    }

    if (type == OUTPUT_TYPE_END) {
        DEBUG_PRINTF("Invalid keyword: '%s'", typekeyword);
        return 1;
    }

    output_source_t source = PRESSURE;
    while (source < OUTPUT_SOURCE_END &&
           strcmp(output_source_keywords[source], sourcekeyword) != 0) {
        source++;
    }

    if (source == OUTPUT_SOURCE_END) {
        DEBUG_PRINTF("Invalid keyword: '%s'", sourcekeyword);
        return 1;
    }

    int readok = 1;
    switch (type) {
        case CUTX:
            readok = (fscanf(fp, "%lf", &posxyz[0]) == 1);
            break;
        case CUTY:
            readok = (fscanf(fp, "%lf", &posxyz[1]) == 1);
            break;
        case CUTZ:
            readok = (fscanf(fp, "%lf", &posxyz[2]) == 1);
            break;
        case ALL:
            break;

        case POINT:
            readok =
                (fscanf(fp, "%lf %lf %lf", &posxyz[0], &posxyz[1], &posxyz[2]) == 3);
            break;

        default:
            break;
    }

    if (readok == 0) {
        DEBUG_PRINT("Failed to read an output parameter");
        return 1;
    }

    sprintf(formatted_filename, "%d_%d_%d_%s", coords[0], coords[1], coords[2], filename);
    output->filename = copy_string(formatted_filename);
    output->type = type;
    output->source = source;
    output->posx = posxyz[0];
    output->posy = posxyz[1];
    output->posz = posxyz[2];

    return 0;
}

int read_sourceparam(FILE* fp, source_t* source) {
    char typekeyword[BUFSZ_SMALL];
    char filename[BUFSZ_LARGE];

    double freq, posx, posy, posz;

    if (fscanf(fp, BUFFMT_SMALL, typekeyword) != 1) {
        DEBUG_PRINT("Failed to read the source parameter");
        return 1;
    }

    source_type_t type = SINE;
    while (type < SOURCE_TYPE_END &&
           strcmp(source_type_keywords[type], typekeyword) != 0) {
        type++;
    }

    if (type == SOURCE_TYPE_END) {
        DEBUG_PRINTF("Invalid keyword: '%s'", typekeyword);
        return 1;
    }

    int readok = 1;
    switch (type) {
        case SINE:
            readok = (fscanf(fp, "%lf", &freq) == 1);
            break;
        case AUDIO:
            readok = (fscanf(fp, BUFFMT_LARGE, filename) == 1);
            break;

        default:
            break;
    }

    if (readok == 0 || fscanf(fp, "%lf %lf %lf", &posx, &posy, &posz) != 3) {
        DEBUG_PRINT("Failed to read the source parameter");
        return 1;
    }

    switch (type) {
        case AUDIO:
            read_audiosource(filename, source);
            break;
        case SINE: {
            if ((source->data = malloc(sizeof(double))) == NULL) {
                DEBUG_PRINT("Failed to allocate memory");
                return 1;
            }

            source->data[0] = freq;
            source->numsamples = 1;

            break;
        }

        default:
            break;
    }

    source->type = type;
    source->posx = posx;
    source->posy = posy;
    source->posz = posz;

    return 0;
}

int read_paramfile(parameters_t* params, const char* filename, int coords[]) {
    if (params == NULL || filename == NULL) {
        DEBUG_PRINT("Invalid print_out params or filename");
        return 1;
    }

    int outrate, numoutputs = 0;

    double dx, dt, maxt;

    char cin_filename[BUFSZ_LARGE];
    char rhoin_filename[BUFSZ_LARGE];

    source_t source;
    output_t* outputs = NULL;

    if ((outputs = malloc(sizeof(output_t) * MAX_OUTPUTS)) == NULL) {
        DEBUG_PRINT("Failed to allocate memory");
        return 1;
    }

    FILE* fp;
    if ((fp = fopen(filename, "r")) == NULL) {
        DEBUG_PRINTF("Could not open parameter file '%s'", filename);
        return 1;
    }

    int readok =
        ((fscanf(fp, "%lf", &dx) == 1) && (fscanf(fp, "%lf", &dt) == 1) &&
         (fscanf(fp, "%lf", &maxt) == 1) && (fscanf(fp, "%d", &outrate) == 1) &&
         (fscanf(fp, BUFFMT_LARGE, cin_filename) == 1) &&
         (fscanf(fp, BUFFMT_LARGE, rhoin_filename) == 1));

    readok = (readok != 0 && read_sourceparam(fp, &source) == 0 &&
              fscanf(fp, " ") == 0);

    while (readok != 0 && numoutputs < MAX_OUTPUTS && feof(fp) == 0) {
        readok = (read_outputparam(fp, &outputs[numoutputs++], coords) == 0 &&
                  fscanf(fp, " ") == 0);
    }

    fclose(fp);

    if (readok == 0) {
        DEBUG_PRINT("Failed to read parameter file");
        free(outputs);
        return 1;
    }

    if (numoutputs == 0) {
        free(outputs);
        outputs = NULL;

    } else if ((outputs = realloc(outputs, sizeof(output_t) * numoutputs)) ==
               NULL) {
        DEBUG_PRINT("Failed to allocate memory");
        return 1;
    }

    params->dx = dx;
    params->dt = dt;
    params->maxt = maxt;
    params->outrate = outrate;
    params->cin_filename = copy_string(cin_filename);
    params->rhoin_filename = copy_string(rhoin_filename);
    params->source = source;
    params->numoutputs = numoutputs;
    params->outputs = outputs;

    return 0;
}

/******************************************************************************
 * Simulation related functions                                               *
 ******************************************************************************/

int interpolate_inputmaps(simulation_data_t* simdata, grid_t* simgrid,
                          data_t* cin, data_t* rhoin) {
    if (simdata == NULL || cin == NULL) {
        DEBUG_PRINT("Invalid NULL simdata or cin");
        return 1;
    }

    if ((simdata->c = allocate_data(simgrid)) == NULL ||
        (simdata->rho = allocate_data(simgrid)) == NULL ||
        (simdata->rhohalf = allocate_data(simgrid)) == NULL) {
        DEBUG_PRINT("Failed to allocate memory");
        return 1;
    }

    double dx = simdata->params.dx;
    double dxd2 = simdata->params.dx / 2;

    for (int m = 0; m < simgrid->numnodesx; m++) {
        for (int n = 0; n < simgrid->numnodesy; n++) {
            for (int p = 0; p < simgrid->numnodesz; p++) {
                double x = simgrid->xmin + m * dx;
                double y = simgrid->ymin + n * dx;
                double z = simgrid->zmin + p * dx;

                double c_approx = trilinear_interpolate(cin, x, y, z);
                double rho_approx = trilinear_interpolate(rhoin, x, y, z);
                SETVALUE(simdata->c, m, n, p, c_approx);
                SETVALUE(simdata->rho, m, n, p, rho_approx);

                x += dxd2;
                y += dxd2;
                z += dxd2;

                double rhohalf_approx = trilinear_interpolate(rhoin, x, y, z);
                SETVALUE(simdata->rhohalf, m, n, p, rhohalf_approx);
            }
        }
    }

    return 0;
}

void apply_source(simulation_data_t* simdata, int step) {
    source_t* source = &simdata->params.source;

    double posx = source->posx;
    double posy = source->posy;
    double posz = source->posz;

    double xmin = simdata->pold->grid.xmin;
    double xmax = simdata->pold->grid.xmax;
    double ymin = simdata->pold->grid.ymin;
    double ymax = simdata->pold->grid.ymax;
    double zmin = simdata->pold->grid.zmin;
    double zmax = simdata->pold->grid.zmax;

    if (posx < xmin || posx > xmax || posy < ymin || posy > ymax || posz < zmin || posz > zmax) {
        // Nothing to do
        return;
    }

    double t = step * simdata->params.dt;

    int m, n, p;
    closest_index(&simdata->pold->grid, posx, posy, posz, &m, &n, &p);

    if (source->type == SINE) {
        double freq = source->data[0];

        SETVALUE(simdata->pold, m, n, p, sin(2 * M_PI * freq * t));

    } else if (source->type == AUDIO) {
        int sample = MIN((int)(t * source->sampling), source->numsamples);

        SETVALUE(simdata->pold, m, n, p, simdata->params.source.data[sample]);
    }
}

void update_pressure(simulation_data_t* simdata) {
    const double dtdx = simdata->params.dt / simdata->params.dx;

    const int numnodesx = NUMNODESX(simdata->pold);
    const int numnodesy = NUMNODESY(simdata->pold);
    const int numnodesz = NUMNODESZ(simdata->pold);

    for (int m = 0; m < numnodesx; m++) {
        for (int n = 0; n < numnodesy; n++) {
            for (int p = 0; p < numnodesz; p++) {
                double c = GETVALUE(simdata->c, m, n, p);
                double rho = GETVALUE(simdata->rho, m, n, p);

                double rhoc2dtdx = rho * c * c * dtdx;

                double dvx = GETVALUE(simdata->vxold, m, n, p);
                double dvy = GETVALUE(simdata->vyold, m, n, p);
                double dvz = GETVALUE(simdata->vzold, m, n, p);

                dvx -= m > 0 ? GETVALUE(simdata->vxold, m - 1, n, p) : *read_from_buffer(simdata->v_buf_old, -1, n, p);
                dvy -= n > 0 ? GETVALUE(simdata->vyold, m, n - 1, p) : *read_from_buffer(simdata->v_buf_old, m, -1, p);
                dvz -= p > 0 ? GETVALUE(simdata->vzold, m, n, p - 1) : *read_from_buffer(simdata->v_buf_old, m, n, -1);

                double prev_p = GETVALUE(simdata->pold, m, n, p);

                SETVALUE(simdata->pnew, m, n, p,
                         prev_p - rhoc2dtdx * (dvx + dvy + dvz));
            }
        }
    }
}

void update_velocities(simulation_data_t* simdata) {
    const double dtdx = simdata->params.dt / simdata->params.dx;

    const int numnodesx = NUMNODESX(simdata->vxold);
    const int numnodesy = NUMNODESY(simdata->vxold);
    const int numnodesz = NUMNODESZ(simdata->vxold);

    for (int m = 0; m < numnodesx; m++) {
        for (int n = 0; n < numnodesy; n++) {
            for (int p = 0; p < numnodesz; p++) {
                double dtdxrho = dtdx / GETVALUE(simdata->rhohalf, m, n, p);

                double p_mnq = GETVALUE(simdata->pold, m, n, p);

                double p_x = m + 1 >= numnodesx ? GETVALUE(simdata->pold, m + 1, n, p) : *read_from_buffer(simdata->p_buf_old, numnodesx, n, p);
                double p_y = n + 1 >= numnodesy ? GETVALUE(simdata->pold, m, n + 1, p) : *read_from_buffer(simdata->p_buf_old, m, numnodesy, p);
                double p_z = p + 1 >= numnodesz ? GETVALUE(simdata->pold, m, n, p + 1) : *read_from_buffer(simdata->p_buf_old, m, n, numnodesz);

                double dpx = p_x - p_mnq;
                double dpy = p_y - p_mnq;
                double dpz = p_z - p_mnq;

                double prev_vx = GETVALUE(simdata->vxold, m, n, p);
                double prev_vy = GETVALUE(simdata->vyold, m, n, p);
                double prev_vz = GETVALUE(simdata->vzold, m, n, p);

                SETVALUE(simdata->vxnew, m, n, p, prev_vx - dtdxrho * dpx);
                SETVALUE(simdata->vynew, m, n, p, prev_vy - dtdxrho * dpy);
                SETVALUE(simdata->vznew, m, n, p, prev_vz - dtdxrho * dpz);
            }
        }
    }
}

void swap_p_timesteps(simulation_data_t* simdata) {
    data_t* tmpp = simdata->pold;
    simdata->pold = simdata->pnew;
    simdata->pnew = tmpp;
}

void swap_v_timesteps(simulation_data_t* simdata) {
    data_t* tmpvx = simdata->vxold;
    data_t* tmpvy = simdata->vyold;
    data_t* tmpvz = simdata->vzold;

    simdata->vxold = simdata->vxnew;
    simdata->vxnew = tmpvx;
    simdata->vyold = simdata->vynew;
    simdata->vynew = tmpvy;
    simdata->vzold = simdata->vznew;
    simdata->vznew = tmpvz;
}

// Segfault handler
void handler(int sig) {
    void* array[10];
    size_t size;

    // get void*'s for all entries on the stack
    size = backtrace(array, 10);

    // print out all the frames to stderr
    printf("Error: signal %d:\n", sig);
    backtrace_symbols_fd(array, size, STDERR_FILENO);
    exit(1);
}

int main(int argc, const char* argv[]) {
    signal(SIGSEGV, handler);  // install our handler
    if (argc < 2) {
        printf("\nUsage: ./fdtd <param_file>\n\n");
        exit(1);
    }

    MPI_Init(NULL, NULL);

    int num_processes;
    MPI_Comm_size(MPI_COMM_WORLD, &num_processes);
    DEBUG_PRINTF("Using %d processes, %d dims", num_processes, NUM_DIMS);

    // dims must be initialized to 0 or their upper bound to not throw
    int dims[NUM_DIMS] = {0, 0, 0};
    MPI_Dims_create(num_processes, NUM_DIMS, dims);

    DEBUG_PRINTF("Created a %dx%dx%d simulation", dims[0], dims[1], dims[2]);

    int periods[NUM_DIMS] = {0, 0, 0};
    MPI_Comm cart_comm;
    MPI_Cart_create(MPI_COMM_WORLD, NUM_DIMS, dims, periods, 0, &cart_comm);

    int cart_rank;
    int coords[3];

    MPI_Comm_rank(cart_comm, &cart_rank);
    MPI_Cart_coords(cart_comm, cart_rank, NUM_DIMS, coords);

    simulation_data_t simdata;
    init_simulation(&simdata, argv[1], coords, dims);

    MPI_Request px_send_req,
        py_send_req,
        pz_send_req,
        vx_send_req,
        vy_send_req,
        vz_send_req,
        px_recv_req,
        py_recv_req,
        pz_recv_req,
        vx_recv_req,
        vy_recv_req,
        vz_recv_req;

    int numtimesteps = floor(simdata.params.maxt / simdata.params.dt);

    double start = GET_TIME();
    for (int tstep = 0; tstep <= numtimesteps; tstep++) {
        apply_source(&simdata, tstep);

        if (simdata.params.outrate > 0 && (tstep % simdata.params.outrate) == 0) {
            for (int i = 0; i < simdata.params.numoutputs; i++) {
                data_t* output_data = NULL;

                switch (simdata.params.outputs[i].source) {
                    case PRESSURE:
                        output_data = simdata.pold;
                        break;
                    case VELOCITYX:
                        output_data = simdata.vxold;
                        break;
                    case VELOCITYY:
                        output_data = simdata.vyold;
                        break;
                    case VELOCITYZ:
                        output_data = simdata.vzold;
                        break;

                    default:
                        break;
                }

                double time = tstep * simdata.params.dt;
                write_output(&simdata.params.outputs[i], output_data, tstep, time);
            }
        }

        if (tstep > 0 && tstep % (numtimesteps / 10) == 0) {
            printf("step %8d/%d", tstep, numtimesteps);

            if (tstep != numtimesteps) {
                double elapsed_sofar = GET_TIME() - start;
                double timeperstep_sofar = elapsed_sofar / tstep;

                double eta = (numtimesteps - tstep) * timeperstep_sofar;

                printf(" (ETA: %8.3lf seconds)", eta);
            }

            printf("\n");
            fflush(stdout);
        }

        rotate_v_recv_request(tstep, &vx_recv_req, &vy_recv_req, &vz_recv_req, simdata.v_buf_old, simdata.v_buf_new, cart_comm);
        update_pressure(&simdata);
        swap_p_timesteps(&simdata);
        copy_send_p_data_to_buffers(&simdata);
        rotate_p_send_request(tstep, &px_send_req, &py_send_req, &pz_send_req, simdata.p_buf_old, simdata.p_buf_new, cart_comm);

        rotate_p_recv_request(tstep, &px_recv_req, &py_recv_req, &pz_recv_req, simdata.p_buf_old, simdata.p_buf_new, cart_comm);
        update_velocities(&simdata);
        swap_v_timesteps(&simdata);
        copy_send_v_data_to_buffers(&simdata);
        rotate_v_send_request(tstep, &vx_send_req, &vy_send_req, &vz_send_req, simdata.v_buf_old, simdata.v_buf_new, cart_comm);
    }

    DEBUG_PRINT("All good, time to clean up");

    double elapsed = GET_TIME() - start;
    double numupdates =
        (double)NUMNODESTOT(simdata.pold->grid) * (numtimesteps + 1);
    double updatespers = numupdates / elapsed / 1e6;

    printf("\nElapsed %.6lf seconds (%.3lf Mupdates/s)\n\n", elapsed,
           updatespers);

    finalize_simulation(&simdata);
    DEBUG_PRINT("Only thing missing is MPI finalize");
    MPI_Finalize();

    return 0;
}
