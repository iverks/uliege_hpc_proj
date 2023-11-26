#pragma once
#include <stdio.h>

#include "output.h"
#include "types.h"

#define MAX_OUTPUTS 32

#define BUFSZ_LARGE 256
#define BUFSZ_HUMONGOUS 512
#define BUFSZ_SMALL 16

#define BUFFMT_LARGE "%255s"
#define BUFFMT_SMALL "%15s"

/**
 * @brief Read an audio source from file
 *
 * @param filename [IN] the file to read the audio from
 * @param source [INOUT] the source object in which to store the data read
 * from the audio file
 * @return int 0 if read was a success, returns 1 otherwise
 */
int read_audiosource(char* filename, source_t* source);

/**
 * @brief Read an output specification parameter from a parameter file
 *
 * @param fp [IN] a file handle to an opened parameter file
 * @param output [OUT] an output oject that will store the output specification
 * read from the parameter file
 * @return int 0 if read was a success, returns 1 otherwise
 */
int read_outputparam(FILE* fp, output_t* output, int coords[]);

/**
 * @brief Read a source specification parameter from a parameter file
 *
 * @param fp [IN] a file handle to an opened parameter file
 * @param source [OUT] a source oject that will store the source specification
 * read from the parameter file
 * @return int 0 if read was a success, returns 1 otherwise
 */
int read_sourceparam(FILE* fp, source_t* source);

/**
 * @brief Read the parameter file specified in argument
 *
 * @param params [OUT] a parameter object into which the parameter read from
 * the file will be stored
 * @param filename [IN] path to the parameter file to read
 * @return int 0 if read was a success, returns 1 otherwise
 */
int read_paramfile(parameters_t* params, const char* filename, int coords[]);

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
