#pragma once
#include <stdio.h>

#include "simulation.h"
#include "types.h"

#define NUMNODESX(dat) ((dat)->grid.numnodesx)
#define NUMNODESY(dat) ((dat)->grid.numnodesy)
#define NUMNODESZ(dat) ((dat)->grid.numnodesz)

#define XMIN(dat) ((dat)->grid.xmin)
#define YMIN(dat) ((dat)->grid.ymin)
#define ZMIN(dat) ((dat)->grid.zmin)

#define XMAX(dat) ((dat)->grid.xmax)
#define YMAX(dat) ((dat)->grid.ymax)
#define ZMAX(dat) ((dat)->grid.zmax)

#define NUMNODESTOT(grid) \
    ((size_t)(grid).numnodesx * (grid).numnodesy * (grid).numnodesz)

#define INDEX3D(grid, m, n, p) \
    ((size_t)(p) + grid.numnodesz * (n) + grid.numnodesz * grid.numnodesy * (m))

#define GETVALUE(dat, m, n, p) ((dat)->vals[INDEX3D((dat)->grid, m, n, p)])
#define SETVALUE(dat, m, n, p, val) \
    ((dat)->vals[INDEX3D((dat)->grid, m, n, p)] = (val))

/**
 * @brief Allocate a new data object and a 3D array to store values according to
 *        the grid provided in argument
 *
 * @param grid [IN] the grid describing the 3D array to allocate
 * @return data_t* return a new data object or NULL if allocation failed
 */
data_t* allocate_data(grid_t* grid);

/**
 * @brief Fills the 3D array of a data object with the value passed in argument
 *
 * @param data  [IN] the data object to fill
 * @param value [IN] the value used to fill the array
 */
void fill_data(data_t* data, double value);

/**
 * @brief Interpolate the input density and spped maps so that they corresponds
 * to the simulation grid
 *
 * @param simdata [OUT] a simulation data object used to store the interpolated
 * values that will be used during the simulation
 * @param simgrid [IN] the simulation grid, i.e., the grid to which the input
 * data will be converted to
 * @param cin [IN] the input speed map
 * @param rhoin [IN] the input density map
 * @return int 0 if read was a success, returns 1 otherwise
 */
int interpolate_inputmaps(simulation_data_t* simdata, grid_t* simgrid,
                          data_t* cin, data_t* rhoin);

/**
 * @brief create a new data file at the path provided in argument. A data file
 *        consists of a header describing the grid/domain followed by and one
 *        or multiple time step data. This function creates the file and writes
 *        the header.
 *
 * @param grid [IN] the grid for the data that will be written to the file
 * @param filename [IN] the path to the file to create
 * @return FILE* file handle to the newly created file or NULL if creation
 * failed
 */
FILE* create_datafile(grid_t grid, char* filename);

/**
 * @brief Write data of a time step to a data file previously opened with the
 *        create_datafile function.
 *
 * @param fp [IN] file handle of the file opened with the create_datafile
 * function
 * @param data [IN] the data to write
 * @param step [IN] the time step index of the step
 * @param time [IN] the time in the simulation of the step
 * @return int 0 if read was a success, returns 1 otherwise
 */
int write_data(FILE* fp, data_t* data, int step, double time);

/**
 * @brief Return the closest indexes in the grid for a given set of x, y and z
 *        coordinates.
 *
 * @param grid [IN] the computational grid
 * @param x [IN] the x coordinate
 * @param y [IN] the y coordinate
 * @param z [IN] the z coordinate
 * @param cx [OUT] the closest x index
 * @param cy [OUT] the closest y index
 * @param cz [OUT] the closest z index
 */
void closest_index(grid_t* grid, double x, double y, double z, int* cx, int* cy,
                   int* cz);

/**
 * @brief Return an approximation of dat at position x, y, z using trilinear interpolation.
 *
 * @param dat [IN] the grid and the data
 * @param x [IN] the x coordinate
 * @param y [IN] the y coordinate
 * @param z [IN] the z coordinate
 * @return The approximate value
 */
double trilinear_interpolate(data_t* dat, double x, double y, double z);

/**
 * @brief Return an approximation of dat at position x, y, z using trilinear interpolation.
 *
 * @param dat [IN] the grid and the data
 * @param x [IN] the x coordinate
 * @param y [IN] the y coordinate
 * @param z [IN] the z coordinate
 * @return The approximate value
 */
double trilinear_interpolate(data_t* dat, double x, double y, double z);