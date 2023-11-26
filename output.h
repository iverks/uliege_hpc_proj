#pragma once

#include "data.h"
#include "parameters.h"
#include "types.h"

/**
 * @brief Print information about the source object passed in argument
 *
 * @param source [IN] the source
 */
void print_source(source_t* source);

/**
 * @brief Print information about the output object passed in argument
 *
 * @param output [IN] the output
 */
void print_output(output_t* output);

/**
 * @brief Write output data to a file
 *
 * @param output [IN] an output object describing which part of the data source
 * needs to be written to the output file
 * @param data  [IN] the data source for the output
 * @param step [IN] the time step index of the step
 * @param time [IN] the time in the simulation of the step
 * @return int 0 if read was a success, returns 1 otherwise
 */
int write_output(output_t* output, data_t* data, int step, double time);

/**
 * @brief Open an output file for writing output data
 * @param output [INOUT] an output object describing which part of the data
 * source needs to be written to the output file. Upon successful return the
 * output file handle member of the output object (fp) will be set to the opened
 * file
 * @param simgrid [IN] the grid used for the simulation
 * @return int 0 if read was a success, returns 1 otherwise
 */
int open_outputfile(output_t* output, grid_t* simgrid);