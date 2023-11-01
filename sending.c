#include "fdtd.h"

void allocate_request_handlers(simulation_data_t data, MPI_Request* n_reqs[]) {
    // There is no way we are supposed to send
}

void send_data(data_t dat, MPI_Comm cart_comm) {
    int n1, n2, n3, n4, n5, n6;
    MPI_Cart_shift(cart_comm, 0, 1, &n1, &n2);
    MPI_Cart_shift(cart_comm, 1, 1, &n3, &n4);
    MPI_Cart_shift(cart_comm, 2, 1, &n5, &n6);

    for (int z : [ 1, dat.grid.numnodesz - 2 ]) {
        for (int x = 1; x < dat.grid.numnodesx - 1; x++) {
            for (int y = 1; y < dat.grid.numnodesy - 1; y++) {
                int pos = x + y * dat.grid.numnodesx;
                // Sending NORTH should have the same number as recieving from SOUTH
                MPI_Isend(GETVALUE(dat, x, y, z), 1, MPI_DOUBLE, n1, pos, MPI_COMM_WORLD, MPI_)
            }
        }
    }
}
