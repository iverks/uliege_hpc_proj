#pragma once

#include "fdtd.h"

void allocate_request_handlers(simulation_data_t data, MPI_Request* n_reqs[]);

void send_data(data_t dat, MPI_Comm cart_comm);
