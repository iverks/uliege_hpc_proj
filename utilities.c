#include "utilities.h"

#include <assert.h>
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

double* buffer_index(buffer_t* dat, int x, int y, buffer_direction_t b_dir) {
    int indexof_buffer = b_dir;
    double* buffer = dat->buffers[indexof_buffer];
    int x_dim_size = NUMNODESX(dat);
    if (b_dir == X_MAX || b_dir == X_MIN) {
        x_dim_size = NUMNODESY(dat);
    }
    return &buffer[x + y * x_dim_size];
}

void fill_buffers(buffer_t* dat, double val) {
    for (buffer_direction_t b_dir = 0; b_dir < BUFFER_DIR_TYPE_END; b_dir++) {
        int numnodesx = NUMNODESX(dat);
        if (b_dir == X_MAX || b_dir == X_MIN) {
            numnodesx = NUMNODESY(dat);
        }
        int numnodesy = NUMNODESZ(dat);
        if (b_dir == Z_MAX || b_dir == Z_MIN) {
            numnodesy = NUMNODESY(dat);
        }
        for (int y = 0; y < numnodesy; y++) {
            for (int x = 0; x < numnodesx; x++) {
                *buffer_index(dat, x, y, b_dir) = val;
            }
        }
    }
}

double* read_from_buffer(buffer_t* data, int m, int n, int p) {
    int numnodesx = NUMNODESX(data);
    int numnodesy = NUMNODESY(data);
    int numnodesz = NUMNODESZ(data);
    if (m == -1 || m == numnodesx) {
        buffer_direction_t b_dir = m == -1 ? X_MIN : X_MAX;
        assert(-1 < n && n < numnodesy);
        assert(-1 < p && p < numnodesz);
        return buffer_index(data, n, p, b_dir);
    } else if (n == -1 || n == numnodesy) {
        buffer_direction_t b_dir = n == -1 ? Y_MIN : Y_MAX;
        assert(-1 < n && n < numnodesy);
        assert(-1 < p && p < numnodesz);
        return buffer_index(data, m, p, b_dir);
    } else if (p == -1 || p == numnodesz) {
        buffer_direction_t b_dir = p == -1 ? Z_MIN : Z_MAX;
        assert(-1 < n && n < numnodesy);
        assert(-1 < p && p < numnodesz);
        return buffer_index(data, m, n, b_dir);
    } else {
        DEBUG_PRINT("CANNOT READ FROM BUFFER IN GRID")
        exit(1);
    }
}

void free_buffers(buffer_t* dat) {
    int numbuffers = BUFFER_DIR_TYPE_END;
    for (int buffer_idx = 0; buffer_idx < numbuffers; buffer_idx++) {
        free(dat->buffers[buffer_idx]);
    }
    free(dat->buffers);
}

int get_buffer_size(buffer_t* data, buffer_direction_t b_dir) {
    int numnodesx = NUMNODESX(data);
    if (b_dir == X_MAX || b_dir == X_MIN) {
        numnodesx = NUMNODESY(data);
    }
    int numnodesy = NUMNODESZ(data);
    if (b_dir == Z_MAX || b_dir == Z_MIN) {
        numnodesy = NUMNODESY(data);
    }
    int numnodestot = numnodesx * numnodesy;
    return numnodestot;
}