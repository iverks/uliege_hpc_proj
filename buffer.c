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
        if (b_dir == X_MAX || b_dir == X_MIN) {
            numnodesx = NUMNODESY(data);
        }
        int numnodesy = NUMNODESZ(data);
        if (b_dir == Z_MAX || b_dir == Z_MIN) {
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

const buffer_direction_t b_dirs_min[3] = {X_MIN, Y_MIN, Z_MIN};
const buffer_direction_t b_dirs_max[3] = {X_MAX, Y_MAX, Z_MAX};

void create_v_recv_request(int timestep,
                           MPI_Request* vx_recv_req,
                           MPI_Request* vy_recv_req,
                           MPI_Request* vz_recv_req, buffer_t* v_buf_new, MPI_Comm cart_comm) {
    MPI_Request* reqs[3] = {vx_recv_req, vy_recv_req, vz_recv_req};
    for (int i = 0; i < 3; i++) {
        int data_size = get_buffer_size(v_buf_new, b_dirs_min[i]);
        int neighbour, dont_care;
        MPI_Cart_shift(cart_comm, i, 1, &neighbour, &dont_care);
        MPI_Irecv(v_buf_new->buffers[b_dirs_min[i]], data_size, MPI_DOUBLE, neighbour, timestep, MPI_COMM_WORLD, reqs[i]);
    }
}

// Remember we recv v from -
void rotate_v_recv_request(int timestep,
                           MPI_Request* vx_recv_req,
                           MPI_Request* vy_recv_req,
                           MPI_Request* vz_recv_req, buffer_t* v_buf_old, buffer_t* v_buf_new, MPI_Comm cart_comm) {
    // Await old reqs
    MPI_Request wait_reqs[3] = {*vx_recv_req, *vy_recv_req, *vz_recv_req};
    if (timestep != 0) {
        MPI_Waitall(3, wait_reqs, MPI_STATUSES_IGNORE);
    }

    // Shuffle buffers, be careful to not include send buffers.
    for (int i = 0; i < 3; i++) {
        double* tmp = v_buf_old->buffers[b_dirs_min[i]];
        v_buf_old->buffers[b_dirs_min[i]] = v_buf_new->buffers[b_dirs_min[i]];
        v_buf_new->buffers[b_dirs_min[i]] = tmp;
    }
    // Create new reqs, here timestep is the timestep from which we send the data.
    create_v_recv_request(timestep, vx_recv_req, vy_recv_req, vz_recv_req, v_buf_new, cart_comm);
}

// Remember we only send v towards +
void rotate_v_send_request(int timestep,
                           MPI_Request* vx_send_req,
                           MPI_Request* vy_send_req,
                           MPI_Request* vz_send_req, buffer_t* v_buf_old, buffer_t* v_buf_new, MPI_Comm cart_comm) {
    // Await old reqs
    MPI_Request* reqs[3] = {vx_send_req, vy_send_req, vz_send_req};
    MPI_Request wait_reqs[3] = {*reqs[0], *reqs[1], *reqs[2]};
    if (timestep != 0) {
        MPI_Waitall(3, wait_reqs, MPI_STATUSES_IGNORE);
    }

    // Shuffle buffers, be careful to not include send buffers.
    for (int i = 0; i < 3; i++) {
        double* tmp = v_buf_old->buffers[b_dirs_max[i]];
        v_buf_old->buffers[b_dirs_max[i]] = v_buf_new->buffers[b_dirs_max[i]];
        v_buf_new->buffers[b_dirs_max[i]] = tmp;
    }
    // Create new reqs
    for (int i = 0; i < 3; i++) {
        int data_size = get_buffer_size(v_buf_new, b_dirs_max[i]);
        int neighbour, dont_care;
        MPI_Cart_shift(cart_comm, i, 1, &dont_care, &neighbour);
        MPI_Isend(v_buf_new->buffers[b_dirs_max[i]], data_size, MPI_DOUBLE, neighbour, timestep, MPI_COMM_WORLD, reqs[i]);
    }
}

void create_p_recv_request(int timestep,
                           MPI_Request* px_recv_req,
                           MPI_Request* py_recv_req,
                           MPI_Request* pz_recv_req, buffer_t* p_buf_new, MPI_Comm cart_comm) {
    MPI_Request* reqs[3] = {px_recv_req, py_recv_req, pz_recv_req};

    for (int i = 0; i < 3; i++) {
        int data_size = get_buffer_size(p_buf_new, b_dirs_max[i]);
        int neighbour, dont_care;
        MPI_Cart_shift(cart_comm, i, 1, &dont_care, &neighbour);
        MPI_Irecv(p_buf_new->buffers[b_dirs_max[i]], data_size, MPI_DOUBLE, neighbour, timestep, MPI_COMM_WORLD, reqs[i]);
    }
}

// Remember we recv p from +
void rotate_p_recv_request(int timestep,
                           MPI_Request* px_recv_req,
                           MPI_Request* py_recv_req,
                           MPI_Request* pz_recv_req, buffer_t* p_buf_old, buffer_t* p_buf_new, MPI_Comm cart_comm) {
    // Await old reqs
    MPI_Request* reqs[3] = {px_recv_req, py_recv_req, pz_recv_req};
    MPI_Request wait_reqs[3] = {*reqs[0], *reqs[1], *reqs[2]};
    // int cart_rank;
    // MPI_Comm_rank(cart_comm, &cart_rank);
    // DEBUG_PRINTF("Waiting for p at ts %d, %d", timestep, cart_rank);
    MPI_Waitall(3, wait_reqs, MPI_STATUSES_IGNORE);
    // DEBUG_PRINTF("DONE    for p at ts %d, %d", timestep, cart_rank);

    // Shuffle buffers, be careful to not include send buffers.
    for (int i = 0; i < 3; i++) {
        double* tmp = p_buf_old->buffers[b_dirs_max[i]];
        p_buf_old->buffers[b_dirs_max[i]] = p_buf_new->buffers[b_dirs_max[i]];
        p_buf_new->buffers[b_dirs_max[i]] = tmp;
    }
    // Create new reqs
    create_p_recv_request(timestep + 1, px_recv_req, py_recv_req, pz_recv_req, p_buf_new, cart_comm);
}

// Remember we send p towards -
void rotate_p_send_request(int timestep,
                           MPI_Request* px_send_req,
                           MPI_Request* py_send_req,
                           MPI_Request* pz_send_req, buffer_t* p_buf_old, buffer_t* p_buf_new, MPI_Comm cart_comm) {
    // Await old reqs
    MPI_Request* reqs[3] = {px_send_req, py_send_req, pz_send_req};
    MPI_Request wait_reqs[3] = {*reqs[0], *reqs[1], *reqs[2]};
    if (timestep != 0) {
        MPI_Waitall(3, wait_reqs, MPI_STATUSES_IGNORE);
    }

    // Shuffle buffers, be careful to not include send buffers.
    for (int i = 0; i < 3; i++) {
        double* tmp = p_buf_old->buffers[b_dirs_min[i]];
        p_buf_old->buffers[b_dirs_min[i]] = p_buf_new->buffers[b_dirs_min[i]];
        p_buf_new->buffers[b_dirs_min[i]] = tmp;
    }
    // Create new reqs
    for (int i = 0; i < 3; i++) {
        int data_size = get_buffer_size(p_buf_new, b_dirs_min[i]);
        int neighbour, dont_care;
        MPI_Cart_shift(cart_comm, i, 1, &neighbour, &dont_care);
        MPI_Isend(p_buf_new->buffers[b_dirs_min[i]], data_size, MPI_DOUBLE, neighbour, timestep, MPI_COMM_WORLD, reqs[i]);
    }
}

// Remember we only send v towards +
void copy_send_v_data_to_buffers(simulation_data_t* simdata) {
    const int numnodesx = NUMNODESX(simdata->vxold);
    const int numnodesy = NUMNODESY(simdata->vxold);
    const int numnodesz = NUMNODESZ(simdata->vxold);
    // We send vx towards x, vy towards y, vz towards z
    int m = numnodesx;
    for (int n = 0; n < numnodesy; n++) {
        for (int p = 0; p < numnodesz; p++) {
            double prev_vx = GETVALUE(simdata->vxold, m - 1, n, p);
            *read_from_buffer(simdata->v_buf_old, m, n, p) = prev_vx;
        }
    }

    int n = numnodesy;
    for (int m = 0; m < numnodesx; m++) {
        for (int p = 0; p < numnodesz; p++) {
            double prev_vy = GETVALUE(simdata->vyold, m, n - 1, p);
            *read_from_buffer(simdata->v_buf_old, m, n, p) = prev_vy;
        }
    }

    int p = numnodesz;
    for (int m = 0; m < numnodesx; m++) {
        for (int n = 0; n < numnodesy; n++) {
            double prev_vz = GETVALUE(simdata->vzold, m, n, p - 1);
            *read_from_buffer(simdata->v_buf_old, m, n, p) = prev_vz;
        }
    }
}

// Remember we send p towards -
void copy_send_p_data_to_buffers(simulation_data_t* simdata) {
    const int numnodesx = NUMNODESX(simdata->vxold);
    const int numnodesy = NUMNODESY(simdata->vxold);
    const int numnodesz = NUMNODESZ(simdata->vxold);
    int m = -1;
    for (int n = 0; n < numnodesy; n++) {
        for (int p = 0; p < numnodesz; p++) {
            double prev_p = GETVALUE(simdata->pold, m + 1, n, p);
            *read_from_buffer(simdata->p_buf_old, m, n, p) = prev_p;
        }
    }

    int n = -1;
    for (int m = 0; m < numnodesx; m++) {
        for (int p = 0; p < numnodesz; p++) {
            double prev_p = GETVALUE(simdata->pold, m, n + 1, p);
            *read_from_buffer(simdata->p_buf_old, m, n, p) = prev_p;
        }
    }

    int p = -1;
    for (int m = 0; m < numnodesx; m++) {
        for (int n = 0; n < numnodesy; n++) {
            double prev_p = GETVALUE(simdata->pold, m, n, p + 1);
            *read_from_buffer(simdata->p_buf_old, m, n, p) = prev_p;
        }
    }
}

double* buffer_index(buffer_t* dat, int x, int y, buffer_direction_t b_dir) {
    int indexof_buffer = b_dir;
    double* buffer = dat->buffers[indexof_buffer];
    int x_dim_size = NUMNODESX(dat);
    if (b_dir == X_MAX || b_dir == X_MIN) {
        x_dim_size = NUMNODESY(dat);
    }
    return &buffer[x + y * x_dim_size];
}

void fill_buffers(buffer_t* dat, double val) {
    for (buffer_direction_t b_dir = 0; b_dir < BUFFER_DIR_TYPE_END; b_dir++) {
        int numnodesx = NUMNODESX(dat);
        if (b_dir == X_MAX || b_dir == X_MIN) {
            numnodesx = NUMNODESY(dat);
        }
        int numnodesy = NUMNODESZ(dat);
        if (b_dir == Z_MAX || b_dir == Z_MIN) {
            numnodesy = NUMNODESY(dat);
        }
        for (int y = 0; y < numnodesy; y++) {
            for (int x = 0; x < numnodesx; x++) {
                *buffer_index(dat, x, y, b_dir) = val;
            }
        }
    }
}

double* read_from_buffer(buffer_t* data, int m, int n, int p) {
    int numnodesx = NUMNODESX(data);
    int numnodesy = NUMNODESY(data);
    int numnodesz = NUMNODESZ(data);
    if (m == -1 || m == numnodesx) {
        buffer_direction_t b_dir = m == -1 ? X_MIN : X_MAX;
        assert(-1 < n && n < numnodesy);
        assert(-1 < p && p < numnodesz);
        return buffer_index(data, n, p, b_dir);
    } else if (n == -1 || n == numnodesy) {
        buffer_direction_t b_dir = n == -1 ? Y_MIN : Y_MAX;
        assert(-1 < m && m < numnodesx);
        assert(-1 < p && p < numnodesz);
        return buffer_index(data, m, p, b_dir);
    } else if (p == -1 || p == numnodesz) {
        buffer_direction_t b_dir = p == -1 ? Z_MIN : Z_MAX;
        assert(-1 < m && m < numnodesx);
        assert(-1 < n && n < numnodesy);
        return buffer_index(data, m, n, b_dir);
    } else {
        DEBUG_PRINT("CANNOT READ FROM BUFFER IN GRID");
        exit(1);
    }
    // unreachable
    return NULL;
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
    if (b_dir == X_MAX || b_dir == X_MIN) {
        numnodesx = NUMNODESY(data);
    }
    int numnodesy = NUMNODESZ(data);
    if (b_dir == Z_MAX || b_dir == Z_MIN) {
        numnodesy = NUMNODESY(data);
    }
    int numnodestot = numnodesx * numnodesy;
    return numnodestot;
}