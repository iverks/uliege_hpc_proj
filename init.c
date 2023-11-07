#include "init.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "fdtd.h"

/******************************************************************************
 * Data functions                                                             *
 ******************************************************************************/

data_t* allocate_data(grid_t* grid) {
    size_t numnodes = NUMNODESTOT(*grid);
    if (numnodes <= 0) {
        DEBUG_PRINTF("Invalid number of nodes (%lu)", numnodes);
        return NULL;
    }

    data_t* data;
    if ((data = malloc(sizeof(data_t))) == NULL) {
        DEBUG_PRINT("Failed to allocate memory");
        free(data);
        return NULL;
    }

    if ((data->vals = malloc(numnodes * sizeof(double))) == NULL) {
        DEBUG_PRINT("Failed to allocate memory");
        free(data->vals);
        free(data);
        return NULL;
    }

    data->grid = *grid;

    return data;
}

void fill_data(data_t* data, double value) {
    if (data == NULL) {
        DEBUG_PRINT("Invalid NULL data");
        return;
    }

    for (int m = 0; m < NUMNODESX(data); m++) {
        for (int n = 0; n < NUMNODESY(data); n++) {
            for (int p = 0; p < NUMNODESZ(data); p++) {
                SETVALUE(data, m, n, p, value);
            }
        }
    }
}

/******************************************************************************
 * Data file functions                                                        *
 ******************************************************************************/

FILE* create_datafile(grid_t grid, char* filename) {
    if (filename == NULL) {
        DEBUG_PRINT("Invalid NULL filename");
        return NULL;
    }

    FILE* fp;
    if ((fp = fopen(filename, "wb")) == NULL) {
        DEBUG_PRINTF("Failed to open file '%s'", filename);
        return NULL;
    }

    if (fwrite(&grid.numnodesx, sizeof(int), 1, fp) != 1 ||
        fwrite(&grid.numnodesy, sizeof(int), 1, fp) != 1 ||
        fwrite(&grid.numnodesz, sizeof(int), 1, fp) != 1 ||
        fwrite(&grid.xmin, sizeof(double), 1, fp) != 1 ||
        fwrite(&grid.xmax, sizeof(double), 1, fp) != 1 ||
        fwrite(&grid.ymin, sizeof(double), 1, fp) != 1 ||
        fwrite(&grid.ymax, sizeof(double), 1, fp) != 1 ||
        fwrite(&grid.zmin, sizeof(double), 1, fp) != 1 ||
        fwrite(&grid.zmax, sizeof(double), 1, fp) != 1) {
        DEBUG_PRINTF("Failed to write header of file '%s'", filename);
        fclose(fp);
        return NULL;
    }

    return fp;
}

FILE* open_datafile(grid_t* grid, int* numsteps, char* filename) {
    if (grid == NULL || filename == NULL) {
        DEBUG_PRINT("Invalid NULL grid or filename");
        return NULL;
    }

    FILE* fp;
    if ((fp = fopen(filename, "rb")) == NULL) {
        DEBUG_PRINTF("Failed to open file '%s'", filename);
        return NULL;
    }

    fseek(fp, 0, SEEK_END);
    size_t file_size = ftell(fp);
    rewind(fp);

    if (fread(&grid->numnodesx, sizeof(int), 1, fp) != 1 ||
        fread(&grid->numnodesy, sizeof(int), 1, fp) != 1 ||
        fread(&grid->numnodesz, sizeof(int), 1, fp) != 1 ||
        fread(&grid->xmin, sizeof(double), 1, fp) != 1 ||
        fread(&grid->xmax, sizeof(double), 1, fp) != 1 ||
        fread(&grid->ymin, sizeof(double), 1, fp) != 1 ||
        fread(&grid->ymax, sizeof(double), 1, fp) != 1 ||
        fread(&grid->zmin, sizeof(double), 1, fp) != 1 ||
        fread(&grid->zmax, sizeof(double), 1, fp) != 1) {
        DEBUG_PRINTF("Failed to read header of file '%s'", filename);
        fclose(fp);
        return NULL;
    }

    size_t numnodestot =
        (size_t)grid->numnodesx * grid->numnodesy * grid->numnodesz;

    size_t values_size = numnodestot * sizeof(double);
    size_t stepindex_size = sizeof(int);
    size_t timestamp_size = sizeof(double);
    size_t header_size = 6 * sizeof(double) + 3 * sizeof(int);

    size_t onetimestep_size = values_size + stepindex_size + timestamp_size;
    size_t alltimestep_size = file_size - header_size;

    if (alltimestep_size % onetimestep_size != 0) {
        DEBUG_PRINTF("Data size is inconsistent with number of nodes (%lu, %lu)",
                     alltimestep_size, onetimestep_size);

        fclose(fp);
        return NULL;
    }

    if (numsteps != NULL) {
        *numsteps = (alltimestep_size / onetimestep_size);
    }

    return fp;
}

data_t* read_data(FILE* fp, grid_t* grid, int* step, double* time) {
    if (fp == NULL) {
        DEBUG_PRINT("Invalid NULL file pointer");
        return NULL;
    }

    double ltime;
    int lstep;

    size_t numnodes = NUMNODESTOT(*grid);

    data_t* data;
    if ((data = allocate_data(grid)) == NULL) {
        DEBUG_PRINT("Failed to allocate data");
        return NULL;
    }

    if (fread(&lstep, sizeof(int), 1, fp) != 1 ||
        fread(&ltime, sizeof(double), 1, fp) != 1 ||
        fread(data->vals, sizeof(double), numnodes, fp) != numnodes) {
        DEBUG_PRINT("Failed to read data");
        free(data);
        return NULL;
    }

    if (step != NULL)
        *step = lstep;
    if (time != NULL)
        *time = ltime;

    return data;
}

int write_data(FILE* fp, data_t* data, int step, double time) {
    if (fp == NULL || data == NULL || data->vals == NULL) {
        DEBUG_PRINT("Invalid NULL data or file pointer");
        return 1;
    }

    size_t numnodes = NUMNODESTOT(data->grid);
    if (numnodes <= 0) {
        DEBUG_PRINTF("Invalid number of nodes (%lu)", numnodes);
        return 1;
    }

    if (fwrite(&step, sizeof(int), 1, fp) != 1 ||
        fwrite(&time, sizeof(double), 1, fp) != 1 ||
        fwrite(data->vals, sizeof(double), numnodes, fp) != numnodes) {
        DEBUG_PRINT("Failed to write data");
        return 1;
    }

    return 0;
}

void init_simulation(simulation_data_t* simdata, const char* params_filename, int coords[]) {
    int x_cart = coords[0];
    int y_cart = coords[1];
    int z_cart = coords[2];
    if (read_paramfile(&simdata->params, params_filename) != 0) {
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

    sim_grid.numnodesx = totalnumnodesx / NUM_NODES_PER_DIM + 2;
    sim_grid.numnodesy = totalnumnodesy / NUM_NODES_PER_DIM + 2;
    sim_grid.numnodesz = totalnumnodesz / NUM_NODES_PER_DIM + 2;

    sim_grid.xmin = rhoin_grid.xmin + (sim_grid.numnodesx - 2) * x_cart * simdata->params.dx;
    sim_grid.xmax = rhoin_grid.xmin + (sim_grid.numnodesx - 2) * (x_cart + 1) * simdata->params.dx;
    sim_grid.ymin = rhoin_grid.ymin + (sim_grid.numnodesy - 2) * y_cart * simdata->params.dx;
    sim_grid.ymax = rhoin_grid.ymin + (sim_grid.numnodesy - 2) * (y_cart + 1) * simdata->params.dx;
    sim_grid.zmin = rhoin_grid.zmin + (sim_grid.numnodesz - 2) * z_cart * simdata->params.dx;
    sim_grid.zmax = rhoin_grid.zmin + (sim_grid.numnodesz - 2) * (z_cart + 1) * simdata->params.dx;

    printf("WALLAH, I have %d * %d * %d nodes\n", sim_grid.numnodesx, sim_grid.numnodesy, sim_grid.numnodesz);
    printf("WALLAH, I go from %f to %f and from %f to %f\n", sim_grid.xmin, sim_grid.xmax, sim_grid.ymin, sim_grid.ymax);

    if (interpolate_inputmaps(simdata, &sim_grid, c_map, rho_map) != 0) {
        printf(
            "Error while converting input map to simulation grid. Aborting...\n\n");
        exit(1);
    }

    if (simdata->params.outrate > 0 && simdata->params.outputs != NULL) {
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

            if (open_outputfile(output, &sim_grid) != 0) {
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

    fill_data(simdata->vynew, 0.0);
    fill_data(simdata->vxold, 0.0);
    fill_data(simdata->vynew, 0.0);
    fill_data(simdata->vyold, 0.0);
    fill_data(simdata->vznew, 0.0);
    fill_data(simdata->vzold, 0.0);

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
    if (simdata->params.outputs != NULL) {
        for (int i = 0; i < simdata->params.numoutputs; i++) {
            free(simdata->params.outputs[i].filename);

            if (simdata->params.outrate > 0) {
                fclose(simdata->params.outputs[i].fp);
            }
        }

        free(simdata->params.outputs);
    }

    free(simdata->params.source.data);
    free(simdata->params.cin_filename);
    free(simdata->params.rhoin_filename);

    free(simdata->rho->vals);
    free(simdata->rho);
    free(simdata->rhohalf->vals);
    free(simdata->rhohalf);
    free(simdata->c->vals);
    free(simdata->c);

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