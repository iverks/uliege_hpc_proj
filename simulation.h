#pragma once
#include "data.h"
#include "parameters.h"
#include "types.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/**
 * @brief Initialize the simulation
 *
 * @param simdata [OUT] a simulation data that will be used to store the data
 * used during the simulation
 * @param params_filename [IN] a path to a parameter file to read
 * @param dims [IN] number of slices in each direction
 * in the grid this process corresponds to
 */
void init_simulation(simulation_data_t* simdata, const char* params_filename);

/**
 * @brief Finalize the simulation by deallocating the data used for the
 * simulation
 *
 * @param simdata [INOUT] a simulation data object describing the simulation to
 * finalize
 */
void finalize_simulation(simulation_data_t* simdata);

/**
 * @brief Apply the source to the simulation. The source is applied to the
 * pressure data of the previous step (pold).
 *
 * @param simdata [INOUT] the simulation data. Upon returns the source is
 * applied to the pold member this object
 * @param step [IN] the simulation time step index
 */
void apply_source(simulation_data_t* simdata, int step);

/**
 * @brief Perform the pressure update step
 *
 * @param simdata [INOUT] a simulation data object used to get the input and
 * store result of the update step
 */
void update_pressure(simulation_data_t* simdata);

/**
 * @brief Perform the velocities update step
 *
 * @param simdata [INOUT] a simulation data object used to get the input and
 * store result of the update step
 */
void update_velocities(simulation_data_t* simdata);

/**
 * @brief Swap the time steps data, i.e., make the new time step the old one
 *
 * @param simdata [INOUT] a simulation data object describing the simulation
 */
void swap_p_timesteps(simulation_data_t* simdata);

/**
 * @brief Swap the time steps data, i.e., make the new time step the old one
 *
 * @param simdata [INOUT] a simulation data object describing the simulation
 */
void swap_v_timesteps(simulation_data_t* simdata);