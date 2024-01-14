#include "data.h"

#include <math.h>
#include <stdlib.h>

#include "utilities.h"

data_t* allocate_data(grid_t* grid) {
    size_t numnodes = NUMNODESTOT(*grid);
    if (numnodes <= 0) {
        DEBUG_PRINTF("Invalid number of nodes (%lu)", numnodes);
        return NULL;
    }

    data_t* data;
    if ((data = malloc(sizeof(data_t))) == NULL) {
        DEBUG_PRINT("Failed to allocate memory for data struct");
        free(data);
        return NULL;
    }

    if ((data->vals = malloc(numnodes * sizeof(double))) == NULL) {
        DEBUG_PRINT("Failed to allocate memory for data vals");
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

// int data_from_inputmap(data_t* new_data, grid_t* simgrid, data_t* old_data, double delta_x, double offset) {
//     if (new_data = allocate_data(simgrid) == NULL)
// }

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

void closest_index(grid_t* grid, double x, double y, double z, int* cx, int* cy,
                   int* cz) {
    int m = (int)((x - grid->xmin) / (grid->xmax - grid->xmin) * grid->numnodesx);
    int n = (int)((y - grid->ymin) / (grid->ymax - grid->ymin) * grid->numnodesy);
    int p = (int)((z - grid->zmin) / (grid->zmax - grid->zmin) * grid->numnodesz);

    *cx = (m < 0) ? 0 : (m > grid->numnodesx - 1) ? grid->numnodesx - 1
                                                  : m;
    *cy = (n < 0) ? 0 : (n > grid->numnodesy - 1) ? grid->numnodesy - 1
                                                  : n;
    *cz = (p < 0) ? 0 : (p > grid->numnodesz - 1) ? grid->numnodesz - 1
                                                  : p;
}

double trilinear_interpolate(data_t* dat, double x, double y, double z) {
    grid_t* grid = &dat->grid;
    double x_on_grid = (x - grid->xmin) / (grid->xmax - grid->xmin) * grid->numnodesx;
    double y_on_grid = (y - grid->ymin) / (grid->ymax - grid->ymin) * grid->numnodesy;
    double z_on_grid = (z - grid->zmin) / (grid->zmax - grid->zmin) * grid->numnodesz;

    int lx = (int)floor(x_on_grid);
    int ly = (int)floor(y_on_grid);
    int lz = (int)floor(z_on_grid);

    // The "bottom" corner must be a full point away from the end
    lx = (lx < 0) ? 0 : lx;
    lx = (lx > grid->numnodesx - 2) ? grid->numnodesx - 2 : lx;
    ly = (ly < 0) ? 0 : ly;
    ly = (ly > grid->numnodesy - 2) ? grid->numnodesy - 2 : ly;
    lz = (lz < 0) ? 0 : lz;
    lz = (lz > grid->numnodesz - 2) ? grid->numnodesz - 2 : lz;

    // Relative distances from "origin" of the cube
    double x_d = (x_on_grid - (double)lx);
    double y_d = (y_on_grid - (double)ly);
    double z_d = (z_on_grid - (double)lz);

    double c_000 = GETVALUE(dat, lx, ly, lz);
    double c_100 = GETVALUE(dat, lx + 1, ly, lz);
    double c_001 = GETVALUE(dat, lx, ly, lz + 1);
    double c_101 = GETVALUE(dat, lx + 1, ly, lz + 1);
    double c_010 = GETVALUE(dat, lx, ly + 1, lz);
    double c_110 = GETVALUE(dat, lx + 1, ly + 1, lz);
    double c_011 = GETVALUE(dat, lx, ly + 1, lz + 1);
    double c_111 = GETVALUE(dat, lx + 1, ly + 1, lz + 1);

    double c_00 = c_000 * (1 - x_d) + c_100 * x_d;
    double c_01 = c_001 * (1 - x_d) + c_101 * x_d;
    double c_10 = c_010 * (1 - x_d) + c_110 * x_d;
    double c_11 = c_011 * (1 - x_d) + c_111 * x_d;

    double c_0 = c_00 * (1 - y_d) + c_10 * y_d;
    double c_1 = c_01 * (1 - y_d) + c_11 * y_d;

    double c = c_0 * (1 - z_d) + c_1 * z_d;

    return c;
}
