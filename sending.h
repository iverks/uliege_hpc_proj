#pragma once

#include <mpi.h>

#include "fdtd.h"

void create_v_recv_request(int timestep,
                           MPI_Request* vx_recv_req,
                           MPI_Request* vy_recv_req,
                           MPI_Request* vz_recv_req, buffer_t* v_buf_new, MPI_Comm cart_comm);

void rotate_v_recv_request(int timestep,
                           MPI_Request* vx_recv_req,
                           MPI_Request* vy_recv_req,
                           MPI_Request* vz_recv_req, buffer_t* v_buf_old, buffer_t* v_buf_new, MPI_Comm cart_comm);

void rotate_v_send_request(int timestep,
                           MPI_Request* vx_send_req,
                           MPI_Request* vy_send_req,
                           MPI_Request* vz_send_req, buffer_t* v_buf_old, buffer_t* v_buf_new, MPI_Comm cart_comm);

void create_p_recv_request(int timestep,
                           MPI_Request* px_recv_req,
                           MPI_Request* py_recv_req,
                           MPI_Request* pz_recv_req, buffer_t* p_buf_new, MPI_Comm cart_comm);

void rotate_p_recv_request(int timestep,
                           MPI_Request* px_recv_req,
                           MPI_Request* py_recv_req,
                           MPI_Request* pz_recv_req, buffer_t* p_buf_old, buffer_t* p_buf_new, MPI_Comm cart_comm);

void rotate_p_send_request(int timestep,
                           MPI_Request* px_send_req,
                           MPI_Request* py_send_req,
                           MPI_Request* pz_send_req, buffer_t* p_buf_old, buffer_t* p_buf_new, MPI_Comm cart_comm);

void copy_send_p_data_to_buffers(simulation_data_t* simdata);

void copy_send_v_data_to_buffers(simulation_data_t* simdata);
