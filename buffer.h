#pragma once
#include <mpi.h>

#include "data.h"
#include "simulation.h"
#include "types.h"

/**
 * @brief Allocate a new buffer object and 6 2D arrays to store values according to
 *        the grid provided in argument
 *
 * @param grid [IN] the grid describing the 3D array to allocate buffers for
 * @return buffer_t* return a new data object or NULL if allocation failed
 */
buffer_t* allocate_buffer(grid_t* grid);

void do_mpi_request(int timestep, MPI_Request reqs[], buffer_t* buf, buffer_t* buf_intransmit, send_direction_t send_dir, send_or_rcv_t send_or_rcv, MPI_Comm cart_comm);

void copy_send_p_data_to_buffers(simulation_data_t* simdata);

void copy_send_v_data_to_buffers(simulation_data_t* simdata);

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

int get_buffer_size(buffer_t* data, buffer_direction_t b_dir);