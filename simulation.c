#include "simulation.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "output.h"
#include "utilities.h"

void init_simulation(simulation_data_t* simdata, const char* params_filename, int coords[], int dims[], int cart_rank) {
    int x_coord = coords[0];
    int y_coord = coords[1];
    int z_coord = coords[2];

    int x_dim = dims[0];
    int y_dim = dims[1];
    int z_dim = dims[2];

    if (read_paramfile(&simdata->params, params_filename, coords) != 0) {
        printf("Failed to read parameters. Aborting...\n\n");
        exit(1);
    }

    grid_t rhoin_grid;
    grid_t cin_grid;
    grid_t sim_grid;

    int rho_numstep;
    int c_numstep;

    FILE* rhofp =
        open_datafile(&rhoin_grid, &rho_numstep, simdata->params.rhoin_filename);
    FILE* cfp =
        open_datafile(&cin_grid, &c_numstep, simdata->params.cin_filename);

    if (rhofp == NULL || rho_numstep <= 0) {
        printf("Failed to open the density map file. Aborting...\n\n");
        exit(1);
    }

    if (cfp == NULL || c_numstep <= 0) {
        printf("Failed to open the speed map file. Aborting...\n\n");
        exit(1);
    }

    if (rhoin_grid.xmin != cin_grid.xmin || rhoin_grid.ymin != cin_grid.ymin ||
        rhoin_grid.zmin != cin_grid.zmin || rhoin_grid.xmax != cin_grid.xmax ||
        rhoin_grid.ymax != cin_grid.ymax || rhoin_grid.zmax != cin_grid.zmax) {
        printf("Grids for the density and speed are not the same. Aborting...\n\n");
        exit(1);
    }

    data_t* rho_map = read_data(rhofp, &rhoin_grid, NULL, NULL);
    data_t* c_map = read_data(cfp, &cin_grid, NULL, NULL);

    if (rho_map == NULL || c_map == NULL) {
        printf("Failed to read data from input maps. Aborting...\n\n");
        exit(1);
    }

    fclose(rhofp);
    fclose(cfp);

    // Note: The way this is done makes rounding sometimes remove a single row from the simulation
    // We deemed this a non-problem, since 1 row out of a 100 is very few.
    int totalnumnodesx = MAX(floor((rhoin_grid.xmax - rhoin_grid.xmin) / (simdata->params.dx)), 1);
    int totalnumnodesy = MAX(floor((rhoin_grid.ymax - rhoin_grid.ymin) / (simdata->params.dx)), 1);
    int totalnumnodesz = MAX(floor((rhoin_grid.zmax - rhoin_grid.zmin) / (simdata->params.dx)), 1);
    DEBUG_PRINTF("Total num nodes on grid is %dx%dx%d", totalnumnodesx, totalnumnodesy, totalnumnodesz);

    sim_grid.numnodesx = totalnumnodesx / x_dim;
    sim_grid.numnodesy = totalnumnodesy / y_dim;
    sim_grid.numnodesz = totalnumnodesz / z_dim;
    // Compensate for errors caused by division
    totalnumnodesx = sim_grid.numnodesx * x_dim;
    totalnumnodesy = sim_grid.numnodesy * y_dim;
    totalnumnodesz = sim_grid.numnodesz * z_dim;

    sim_grid.xmin = rhoin_grid.xmin + (sim_grid.numnodesx) * x_coord * simdata->params.dx;
    sim_grid.xmax = rhoin_grid.xmin + (sim_grid.numnodesx) * (x_coord + 1) * simdata->params.dx;
    sim_grid.ymin = rhoin_grid.ymin + (sim_grid.numnodesy) * y_coord * simdata->params.dx;
    sim_grid.ymax = rhoin_grid.ymin + (sim_grid.numnodesy) * (y_coord + 1) * simdata->params.dx;
    sim_grid.zmin = rhoin_grid.zmin + (sim_grid.numnodesz) * z_coord * simdata->params.dx;
    sim_grid.zmax = rhoin_grid.zmin + (sim_grid.numnodesz) * (z_coord + 1) * simdata->params.dx;

    // Write buffer handling for collection
    grid_t write_buffer_grid;
    if (cart_rank == 0) {
        write_buffer_grid.numnodesx = totalnumnodesx;
        write_buffer_grid.numnodesy = totalnumnodesy;
        write_buffer_grid.numnodesz = totalnumnodesz;

        write_buffer_grid.xmin = rhoin_grid.xmin;
        write_buffer_grid.xmax = rhoin_grid.xmin + totalnumnodesx * simdata->params.dx;
        write_buffer_grid.ymin = rhoin_grid.ymin;
        write_buffer_grid.ymax = rhoin_grid.ymin + totalnumnodesy * simdata->params.dx;
        write_buffer_grid.zmin = rhoin_grid.zmin;
        write_buffer_grid.zmax = rhoin_grid.zmin + totalnumnodesz * simdata->params.dx;

        int numnodestot = NUMNODESTOT(write_buffer_grid);

        if ((simdata->write_data = allocate_data(&write_buffer_grid)) == NULL ||
            (simdata->write_data_buffer = malloc(numnodestot * sizeof(double))) == NULL) {
            printf("Failed to allocate memory for collection. Aborting...\n\n");
            fflush(stdout);
            exit(1);
        }
        fill_data(simdata->write_data, 0.0);
        for (int i = 0; i < numnodestot; i++) {
            simdata->write_data_buffer[i] = 0;
        }
    }

    if (interpolate_inputmaps(simdata, &sim_grid, c_map, rho_map) != 0) {
        printf(
            "Error while converting input map to simulation grid. Aborting...\n\n");
        exit(1);
    }

    if (simdata->params.outrate > 0 && simdata->params.outputs != NULL && cart_rank == 0) {
        for (int i = 0; i < simdata->params.numoutputs; i++) {
            char* outfilei = simdata->params.outputs[i].filename;

            for (int j = 0; j < i; j++) {
                char* outfilej = simdata->params.outputs[j].filename;

                if (strcmp(outfilei, outfilej) == 0) {
                    printf("Duplicate output file: '%s'. Aborting...\n\n", outfilei);
                    exit(1);
                }
            }
        }

        for (int i = 0; i < simdata->params.numoutputs; i++) {
            output_t* output = &simdata->params.outputs[i];

            if (open_outputfile(output, &write_buffer_grid) != 0) {
                printf("Failed to open output file: '%s'. Aborting...\n\n",
                       output->filename);
                exit(1);
            }
        }
    }

    if ((simdata->pold = allocate_data(&sim_grid)) == NULL ||
        (simdata->pnew = allocate_data(&sim_grid)) == NULL ||
        (simdata->vxold = allocate_data(&sim_grid)) == NULL ||
        (simdata->vxnew = allocate_data(&sim_grid)) == NULL ||
        (simdata->vyold = allocate_data(&sim_grid)) == NULL ||
        (simdata->vynew = allocate_data(&sim_grid)) == NULL ||
        (simdata->vzold = allocate_data(&sim_grid)) == NULL ||
        (simdata->vznew = allocate_data(&sim_grid)) == NULL) {
        printf("Failed to allocate memory. Aborting...\n\n");
        exit(1);
    }

    fill_data(simdata->pold, 0.0);
    fill_data(simdata->pnew, 0.0);

    fill_data(simdata->vxold, 0.0);
    fill_data(simdata->vynew, 0.0);
    fill_data(simdata->vyold, 0.0);
    fill_data(simdata->vznew, 0.0);
    fill_data(simdata->vzold, 0.0);

    if ((simdata->p_send_buf = allocate_buffer(&sim_grid)) == NULL ||
        (simdata->p_send_buf_intransmit = allocate_buffer(&sim_grid)) == NULL ||
        (simdata->p_recv_buf = allocate_buffer(&sim_grid)) == NULL ||
        (simdata->p_recv_buf_intransmit = allocate_buffer(&sim_grid)) == NULL ||
        (simdata->v_send_buf = allocate_buffer(&sim_grid)) == NULL ||
        (simdata->v_send_buf_intransmit = allocate_buffer(&sim_grid)) == NULL ||
        (simdata->v_recv_buf = allocate_buffer(&sim_grid)) == NULL ||
        (simdata->v_recv_buf_intransmit = allocate_buffer(&sim_grid)) == NULL) {
        printf("Failed to allocate buffer memory. Aborting...\n\n");
        exit(1);
    }

    fill_buffers(simdata->p_send_buf, 0.0);
    fill_buffers(simdata->p_send_buf_intransmit, 0.0);
    fill_buffers(simdata->p_recv_buf, 0.0);
    fill_buffers(simdata->p_recv_buf_intransmit, 0.0);
    fill_buffers(simdata->v_send_buf, 0.0);
    fill_buffers(simdata->v_send_buf_intransmit, 0.0);
    fill_buffers(simdata->v_recv_buf, 0.0);
    fill_buffers(simdata->v_recv_buf_intransmit, 0.0);

    printf("\n");
    printf(" Grid spacing: %g\n", simdata->params.dx);
    printf("  Grid size X: %d\n", sim_grid.numnodesx);
    printf("  Grid size Y: %d\n", sim_grid.numnodesy);
    printf("  Grid size Z: %d\n", sim_grid.numnodesz);
    printf("    Time step: %g\n", simdata->params.dt);
    printf(" Maximum time: %g\n\n", simdata->params.maxt);

    if (simdata->params.outrate > 0 && simdata->params.outputs) {
        int outsampling =
            (int)(1.0 / (simdata->params.outrate * simdata->params.dt));

        printf("     Output rate: every %d step(s)\n", simdata->params.outrate);
        printf(" Output sampling: %d Hz\n\n", outsampling);
        printf(" Output files:\n\n");

        for (int i = 0; i < simdata->params.numoutputs; i++) {
            print_output(&simdata->params.outputs[i]);
        }

        printf("\n");

    } else if (simdata->params.outrate < 0) {
        printf("  Output is disabled (output rate set to 0)\n\n");

    } else {
        printf("  Output is disabled (not output specified)\n\n");
    }

    print_source(&simdata->params.source);

    fflush(stdout);

    free(rho_map->vals);
    free(rho_map);
    free(c_map->vals);
    free(c_map);
}

void finalize_simulation(simulation_data_t* simdata) {
    DEBUG_PRINT("I'm gonna free");

    if (simdata->params.outputs != NULL) {
        for (int i = 0; i < simdata->params.numoutputs; i++) {
            DEBUG_PRINTF("Freeing output %d, %s", i, simdata->params.outputs[i].filename);
            free(simdata->params.outputs[i].filename);

            int cart_rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &cart_rank);
            if (simdata->params.outrate > 0 && cart_rank == 0) {
                fclose(simdata->params.outputs[i].fp);
            }
        }

        free(simdata->params.outputs);
    }

    DEBUG_PRINT("Freed simdata outputs");

    free(simdata->params.source.data);
    free(simdata->params.cin_filename);
    free(simdata->params.rhoin_filename);

    free(simdata->rho->vals);
    free(simdata->rho);
    free(simdata->rhohalf->vals);
    free(simdata->rhohalf);
    free(simdata->c->vals);
    free(simdata->c);

    // TODO: Free buffers
    // free_buffers(simdata->pold);
    free(simdata->pold->vals);
    free(simdata->pold);
    free(simdata->pnew->vals);
    free(simdata->pnew);

    free(simdata->vxold->vals);
    free(simdata->vxold);
    free(simdata->vxnew->vals);
    free(simdata->vxnew);
    free(simdata->vyold->vals);
    free(simdata->vyold);
    free(simdata->vynew->vals);
    free(simdata->vynew);
    free(simdata->vzold->vals);
    free(simdata->vzold);
    free(simdata->vznew->vals);
    free(simdata->vznew);
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

                dvx -= m > 0 ? GETVALUE(simdata->vxold, m - 1, n, p) : *buffer_index(simdata->v_recv_buf, n, p, DIR_X);
                dvy -= n > 0 ? GETVALUE(simdata->vyold, m, n - 1, p) : *buffer_index(simdata->v_recv_buf, m, p, DIR_Y);
                dvz -= p > 0 ? GETVALUE(simdata->vzold, m, n, p - 1) : *buffer_index(simdata->v_recv_buf, m, n, DIR_Z);

                double prev_p = GETVALUE(simdata->pold, m, n, p);

                SETVALUE(simdata->pnew, m, n, p,
                         prev_p - rhoc2dtdx * (dvx + dvy + dvz));
            }
        }
    }
}

void update_velocities(simulation_data_t* simdata, int coords[], int dims[]) {
    const double dtdx = simdata->params.dt / simdata->params.dx;

    const int numnodesx = NUMNODESX(simdata->vxold);
    const int numnodesy = NUMNODESY(simdata->vxold);
    const int numnodesz = NUMNODESZ(simdata->vxold);

    for (int m = 0; m < numnodesx; m++) {
        for (int n = 0; n < numnodesy; n++) {
            for (int p = 0; p < numnodesz; p++) {
                double dtdxrho = dtdx / GETVALUE(simdata->rhohalf, m, n, p);

                double p_mnq = GETVALUE(simdata->pold, m, n, p);

                double p_x = m + 1 < numnodesx ? GETVALUE(simdata->pold, m + 1, n, p) : *buffer_index(simdata->p_recv_buf, n, p, DIR_X);
                double p_y = n + 1 < numnodesy ? GETVALUE(simdata->pold, m, n + 1, p) : *buffer_index(simdata->p_recv_buf, m, p, DIR_Y);
                double p_z = p + 1 < numnodesz ? GETVALUE(simdata->pold, m, n, p + 1) : *buffer_index(simdata->p_recv_buf, m, n, DIR_Z);

                double dpx = p_x - p_mnq;
                if (coords[0] == dims[0] - 1 && m + 1 >= numnodesx) {
                    dpx = 0.0;
                }
                double dpy = p_y - p_mnq;
                if (coords[1] == dims[1] - 1 && n + 1 >= numnodesy) {
                    dpy = 0.0;
                }
                double dpz = p_z - p_mnq;
                if (coords[2] == dims[2] - 1 && p + 1 >= numnodesz) {
                    dpz = 0.0;
                }

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