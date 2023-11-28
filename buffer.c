#include "buffer.h"

#include <assert.h>
#include <stdlib.h>

#include "utilities.h"

buffer_t* allocate_buffer(grid_t* grid) {
    buffer_t* data;

    if ((data = malloc(sizeof(buffer_t))) == NULL) {
        DEBUG_PRINT("Failed to allocate memory for data struct");
        free(data);
        return NULL;
    }

    data->grid = *grid;
    // Buffer code
    int numbuffers = BUFFER_DIR_TYPE_END;
    if ((data->buffers = malloc(numbuffers * sizeof(double*))) == NULL) {
        DEBUG_PRINT("Failed to allocate memory for buffer pointers");
        free(data->buffers);
        free(data);
        return NULL;
    }

    int success = 1;
    for (buffer_direction_t b_dir = 0; b_dir < BUFFER_DIR_TYPE_END; b_dir++) {
        int numnodesx = NUMNODESX(data);
        if (b_dir == DIR_X) {
            numnodesx = NUMNODESY(data);
        }
        int numnodesy = NUMNODESZ(data);
        if (b_dir == DIR_Z) {
            numnodesy = NUMNODESY(data);
        }
        int buffer_idx = b_dir;
        int numnodestot = numnodesx * numnodesy;
        if ((data->buffers[buffer_idx] = malloc(numnodestot * sizeof(double))) == NULL) {
            success = 0;
        }
    }

    if (!success) {
        DEBUG_PRINT("Failed to allocate memory for buffers");
        for (int buffer_idx = 0; buffer_idx < numbuffers; buffer_idx++) {
            free(data->buffers[buffer_idx]);
        }
        free(data->buffers);
        free(data);
        return NULL;
    }

    return data;
}

#define P_TRANSFER_TAG 0
#define V_TRANSFER_TAG 1

const buffer_direction_t b_dirs[3] = {DIR_X, DIR_Y, DIR_Z};

void do_mpi_request(int timestep, MPI_Request reqs[], buffer_t* buf, buffer_t* buf_intransmit, send_direction_t send_dir, send_or_rcv_t send_or_rcv, MPI_Comm cart_comm) {
    if (timestep != 0) {
        MPI_Waitall(3, reqs, MPI_STATUSES_IGNORE);
    }

    // Buffer swap
    {
        buffer_t* tmp = buf;
        buf = buf_intransmit;
        buf_intransmit = tmp;
    }

    for (int i = 0; i < 3; i++) {
        int data_size = get_buffer_size(buf_intransmit, b_dirs[i]);
        int neighbour, dont_care;
        if (send_dir == SEND_POSITIVE) {
            MPI_Cart_shift(cart_comm, i, 1, &dont_care, &neighbour);
        } else {
            MPI_Cart_shift(cart_comm, i, 1, &neighbour, &dont_care);
        }
        if (send_or_rcv == SEND) {
            MPI_Isend(buf_intransmit->buffers[b_dirs[i]], data_size, MPI_DOUBLE, neighbour, timestep, MPI_COMM_WORLD, &reqs[i]);
        } else {
            MPI_Irecv(buf_intransmit->buffers[b_dirs[i]], data_size, MPI_DOUBLE, neighbour, timestep, MPI_COMM_WORLD, &reqs[i]);
        }
    }
}

// Remember we only send v towards +
void copy_send_v_data_to_buffers(simulation_data_t* simdata) {
    const int numnodesx = NUMNODESX(simdata->vxold);
    const int numnodesy = NUMNODESY(simdata->vxold);
    const int numnodesz = NUMNODESZ(simdata->vxold);
    // We send vx towards x, vy towards y, vz towards z
    for (int n = 0; n < numnodesy; n++) {
        for (int p = 0; p < numnodesz; p++) {
            double prev_vx = GETVALUE(simdata->vxold, numnodesx - 1, n, p);
            *buffer_index(simdata->v_send_buf, n, p, DIR_X) = prev_vx;
        }
    }

    for (int m = 0; m < numnodesx; m++) {
        for (int p = 0; p < numnodesz; p++) {
            double prev_vy = GETVALUE(simdata->vyold, m, numnodesy - 1, p);
            *buffer_index(simdata->v_send_buf, m, p, DIR_Y) = prev_vy;
        }
    }

    for (int m = 0; m < numnodesx; m++) {
        for (int n = 0; n < numnodesy; n++) {
            double prev_vz = GETVALUE(simdata->vzold, m, n, numnodesz - 1);
            *buffer_index(simdata->v_send_buf, m, n, DIR_Z) = prev_vz;
        }
    }
}

// Remember we send p towards -
void copy_send_p_data_to_buffers(simulation_data_t* simdata) {
    const int numnodesx = NUMNODESX(simdata->vxold);
    const int numnodesy = NUMNODESY(simdata->vxold);
    const int numnodesz = NUMNODESZ(simdata->vxold);
    for (int n = 0; n < numnodesy; n++) {
        for (int p = 0; p < numnodesz; p++) {
            double prev_p = GETVALUE(simdata->pold, 0, n, p);
            *buffer_index(simdata->p_send_buf, n, p, DIR_X) = prev_p;
        }
    }

    int n = -1;
    for (int m = 0; m < numnodesx; m++) {
        for (int p = 0; p < numnodesz; p++) {
            double prev_p = GETVALUE(simdata->pold, m, n + 1, p);
            *buffer_index(simdata->p_send_buf, m, p, DIR_Y) = prev_p;
        }
    }

    int p = -1;
    for (int m = 0; m < numnodesx; m++) {
        for (int n = 0; n < numnodesy; n++) {
            double prev_p = GETVALUE(simdata->pold, m, n, p + 1);
            *buffer_index(simdata->p_send_buf, m, n, DIR_Z) = prev_p;
        }
    }
}

double* buffer_index(buffer_t* dat, int x, int y, buffer_direction_t b_dir) {
    int indexof_buffer = b_dir;
    double* buffer = dat->buffers[indexof_buffer];
    int x_dim_size = NUMNODESX(dat);
    if (b_dir == DIR_X) {
        x_dim_size = NUMNODESY(dat);
    }
    return &buffer[x + y * x_dim_size];
}

void fill_buffers(buffer_t* dat, double val) {
    for (buffer_direction_t b_dir = 0; b_dir < BUFFER_DIR_TYPE_END; b_dir++) {
        int numnodesx = NUMNODESX(dat);
        if (b_dir == DIR_X) {
            numnodesx = NUMNODESY(dat);
        }
        int numnodesy = NUMNODESZ(dat);
        if (b_dir == DIR_Z) {
            numnodesy = NUMNODESY(dat);
        }
        for (int y = 0; y < numnodesy; y++) {
            for (int x = 0; x < numnodesx; x++) {
                *buffer_index(dat, x, y, b_dir) = val;
            }
        }
    }
}

void free_buffers(buffer_t* dat) {
    int numbuffers = BUFFER_DIR_TYPE_END;
    for (int buffer_idx = 0; buffer_idx < numbuffers; buffer_idx++) {
        free(dat->buffers[buffer_idx]);
    }
    free(dat->buffers);
}

int get_buffer_size(buffer_t* data, buffer_direction_t b_dir) {
    int numnodesx = NUMNODESX(data);
    if (b_dir == DIR_X) {
        numnodesx = NUMNODESY(data);
    }
    int numnodesy = NUMNODESZ(data);
    if (b_dir == DIR_Z) {
        numnodesy = NUMNODESY(data);
    }
    int numnodestot = numnodesx * numnodesy;
    return numnodestot;
}