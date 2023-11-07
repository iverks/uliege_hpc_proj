#include "utilities.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "fdtd.h"

/******************************************************************************
 * Utilities functions                                                        *
 ******************************************************************************/

char* copy_string(char* str) {
    size_t len;
    if (str == NULL || (len = strlen(str)) == 0) {
        DEBUG_PRINT("NULL of zero length string passed as argument");
        return NULL;
    }

    char* cpy;
    if ((cpy = malloc((len + 1) * sizeof(char))) == NULL) {
        DEBUG_PRINT("Failed to allocate memory");
        return NULL;
    }

    return strcpy(cpy, str);
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
    int lx = (int)floor((x - grid->xmin) / (grid->xmax - grid->xmin) * grid->numnodesx);
    int ly = (int)floor((y - grid->ymin) / (grid->ymax - grid->ymin) * grid->numnodesy);
    int lz = (int)floor((z - grid->zmin) / (grid->zmax - grid->zmin) * grid->numnodesz);

    // The "bottom" corner must be a full point away from the end
    lx = (lx < 0) ? 0 : lx;
    lx = (lx > grid->numnodesx - 2) ? grid->numnodesx - 2 : lx;
    ly = (ly < 0) ? 0 : ly;
    ly = (ly > grid->numnodesy - 2) ? grid->numnodesy - 2 : ly;
    lz = (lz < 0) ? 0 : lz;
    lz = (lz > grid->numnodesz - 2) ? grid->numnodesz - 2 : lz;

    // Relative distances from "origin" of the cube
    double x_d = (x - (double)lx) / (grid->xmax - grid->xmin) * grid->numnodesx;
    double y_d = (y - (double)ly) / (grid->ymax - grid->ymin) * grid->numnodesy;
    double z_d = (z - (double)lz) / (grid->zmax - grid->zmin) * grid->numnodesz;

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

/******************************************************************************
 * Buffer utilities                                                           *
 ******************************************************************************/

double* buffer_index(data_t* dat, int x, int y, buffer_direction_t b_dir, buffer_type_t b_type) {
    int indexof_buffer = b_dir + BUFFER_DIR_TYPE_END * b_type;
    double* buffer = dat->buffers[indexof_buffer];
    int x_dim_size = dat->grid.numnodesx;
    if (b_dir == X_MAX || b_dir == X_MIN) {
        x_dim_size = dat->grid.numnodesy;
    }
    return &buffer[x + y * x_dim_size];
}

void fill_buffers(data_t* dat, double val) {
    for (buffer_direction_t b_dir = 0; b_dir < BUFFER_DIR_TYPE_END; b_dir++) {
        for (buffer_type_t b_type = 0; b_type < BUFFER_TYPE_TYPE_END; b_type++) {
            int numnodesx = dat->grid.numnodesx;
            if (b_dir == X_MAX || b_dir == X_MIN) {
                numnodesx = dat->grid.numnodesy;
            }
            int numnodesy = dat->grid.numnodesz;
            if (b_dir == Z_MAX || b_dir == Z_MIN) {
                numnodesy = dat->grid.numnodesy;
            }
            for (int y = 0; y < numnodesy; y++) {
                for (int x = 0; x < numnodesx; x++) {
                    *buffer_index(dat, x, y, b_dir, b_type) = val;
                }
            }
        }
    }
}