#pragma once
#include "fdtd.h"

/*********************************************************************************
 * Utilities functions                                                        *
 ******************************************************************************/

/**
 * @brief Copy the string passed in argument
 *
 * @param str [IN] the string to copy
 * @return char* a copy of the string passed in argument or NULL if memory
 * allocation for the copy failed
 */
char* copy_string(char* str);

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
 * @brief Index into a buffer at position x, y. Convention is index in alphabetical order, e.g
 * for the y-direction plane we index at [x, z]
 * @param
 * @return a pointer to the value at the given index
 */
double* buffer_index(buffer_t* dat, int x, int y, buffer_direction_t b_dir);

/**
 * @brief Fill all buffers of data object with value
 */
void fill_buffers(buffer_t* dat, double val);

/**
 * @brief Correctly frees all buffers of data object
 */
void free_buffers(buffer_t* dat);

double* read_from_buffer(buffer_t* data, int m, int n, int p);

int get_buffer_size(buffer_t* data, buffer_direction_t b_dir);