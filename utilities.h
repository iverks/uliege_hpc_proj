#pragma once
#include "fdtd.h"

/******************************************************************************
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
 */
double trilinear_interpolate(data_t* dat, double x, double y, double z);