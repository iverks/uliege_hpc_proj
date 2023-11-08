#pragma once
#include "fdtd.h"

/******************************************************************************
 * Data functions                                                             *
 ******************************************************************************/

/**
 * @brief Allocate a new data object and a 3D array to store values according to
 *        the grid provided in argument
 *
 * @param grid [IN] the grid describing the 3D array to allocate
 * @return data_t* return a new data object or NULL if allocation failed
 */
data_t* allocate_data(grid_t* grid);

/**
 * @brief Allocate a new buffer object and 6 2D arrays to store values according to
 *        the grid provided in argument
 *
 * @param grid [IN] the grid describing the 3D array to allocate buffers for
 * @return buffer_t* return a new data object or NULL if allocation failed
 */
buffer_t* allocate_buffer(grid_t* grid);

/**
 * @brief Fills the 3D array of a data object with the value passed in argument
 *
 * @param data  [IN] the data object to fill
 * @param value [IN] the value used to fill the array
 */
void fill_data(data_t* data, double value);

/******************************************************************************
 * Data file functions                                                        *
 ******************************************************************************/

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
 * @brief Open an existing data file for reading. This function opens the file
 *        and reads the header.
 *
 * @param grid     [OUT] upon return store the values read from the file header
 * @param numsteps [OUT] upon return, the number of time steps stored in the
 * file
 * @param filename  [IN] path of the file to open
 * @return FILE* file handle to the opened file or NULL if opening the file
 * failed
 */
FILE* open_datafile(grid_t* grid, int* numsteps, char* filename);

/**
 * @brief Read the data for one time step from a data file previously opened
 * with the open_datafile function and return the values for the time step in a
 * data object.
 *
 * @param fp [IN] file handle of the file opened with the open_datafile
 * function
 * @param grid [IN] grid of the data to read from the file. Should be identical
 * to the one returned by the open_datafile function
 * @param step [OUT] upon return, the time step index of the step
 * @param time [OUT] upon return, the time in the simulation of the step
 * @return data_t* a data object with the values of the time step read from the
 *                 file or NULL if read failed
 */
data_t* read_data(FILE* fp, grid_t* grid, int* step, double* time);

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
 * @brief Initialize the simulation
 *
 * @param simdata [OUT] a simulation data that will be used to store the data
 * used during the simulation
 * @param params_filename [IN] a path to a parameter file to read
 * @param coords [IN] cartesian coordinates representing the position
 * in the grid this process corresponds to
 */
void init_simulation(simulation_data_t* simdata, const char* params_filename, int coords[]);

/**
 * @brief Finalize the simulation by deallocating the data used for the
 * simulation
 *
 * @param simdata [INOUT] a simulation data object describing the simulation to
 * finalize
 */
void finalize_simulation(simulation_data_t* simdata);
