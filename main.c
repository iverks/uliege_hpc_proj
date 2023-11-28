#include <assert.h>
#include <execinfo.h>
#include <math.h>
#include <signal.h>
#include <unistd.h>

#include "buffer.h"
#include "data.h"
#include "output.h"
#include "parameters.h"
#include "simulation.h"
#include "types.h"
#include "utilities.h"

#define NUM_DIMS 3

// Segfault handler
void handler(int sig) {
    void* array[10];
    size_t size;

    // get void*'s for all entries on the stack
    size = backtrace(array, 10);

    // print out all the frames to stderr
    printf("Error: signal %d:\n", sig);
    backtrace_symbols_fd(array, size, STDERR_FILENO);
    exit(1);
}

int main(int argc, const char* argv[]) {
    signal(SIGSEGV, handler);  // install our handler
    if (argc < 2) {
        printf("\nUsage: ./fdtd <param_file>\n\n");
        exit(1);
    }

    MPI_Init(NULL, NULL);

    int num_processes;
    MPI_Comm_size(MPI_COMM_WORLD, &num_processes);
    DEBUG_PRINTF("Using %d processes, %d dims", num_processes, NUM_DIMS);

    // dims must be initialized to 0 or their upper bound to not throw
    int dims[NUM_DIMS] = {0, 0, 0};
    MPI_Dims_create(num_processes, NUM_DIMS, dims);

    DEBUG_PRINTF("Created a %dx%dx%d simulation", dims[0], dims[1], dims[2]);

    int periods[NUM_DIMS] = {0, 0, 0};
    MPI_Comm cart_comm;
    MPI_Cart_create(MPI_COMM_WORLD, NUM_DIMS, dims, periods, 0, &cart_comm);

    int cart_rank;
    int coords[3];

    MPI_Comm_rank(cart_comm, &cart_rank);
    MPI_Cart_coords(cart_comm, cart_rank, NUM_DIMS, coords);

    simulation_data_t simdata;
    init_simulation(&simdata, argv[1], coords, dims, cart_rank);

    MPI_Request px_send_req,
        py_send_req,
        pz_send_req,
        vx_send_req,
        vy_send_req,
        vz_send_req,
        px_recv_req,
        py_recv_req,
        pz_recv_req,
        vx_recv_req,
        vy_recv_req,
        vz_recv_req;

    int numtimesteps = floor(simdata.params.maxt / simdata.params.dt);

    double start = GET_TIME();
    // Prepare to recieve already in loop idx 0
    create_p_recv_request(0,
                          &px_recv_req,
                          &py_recv_req,
                          &pz_recv_req, simdata.p_buf_new, cart_comm);

    for (int tstep = 0; tstep <= numtimesteps; tstep++) {
        apply_source(&simdata, tstep);

        if (simdata.params.outrate > 0 && (tstep % simdata.params.outrate) == 0) {
            for (int i = 0; i < simdata.params.numoutputs; i++) {
                data_t* output_data = NULL;

                switch (simdata.params.outputs[i].source) {
                    case PRESSURE:
                        output_data = simdata.pold;
                        break;
                    case VELOCITYX:
                        output_data = simdata.vxold;
                        break;
                    case VELOCITYY:
                        output_data = simdata.vyold;
                        break;
                    case VELOCITYZ:
                        output_data = simdata.vzold;
                        break;

                    default:
                        break;
                }

                double time = tstep * simdata.params.dt;
                MPI_Gather(output_data->vals, NUMNODESTOT(output_data->grid), MPI_DOUBLE, simdata.write_data_buffer, NUMNODESTOT(output_data->grid), MPI_DOUBLE, 0, MPI_COMM_WORLD);

                if (cart_rank == 0) {
                    for (size_t j = 0; j < NUMNODESTOT(simdata.write_data->grid); j++) {
                        // Which rank am i recieveing from
                        int src_idx = j / NUMNODESTOT(output_data->grid);
                        int src_coords[3];
                        MPI_Cart_coords(cart_comm, src_idx, NUM_DIMS, src_coords);
                        int src_x = src_coords[0];
                        int src_y = src_coords[1];
                        int src_z = src_coords[2];
                        // Which index in that given rank
                        int local_idx = j % NUMNODESTOT(output_data->grid);
                        // ((size_t)(z) + grid.numnodesz * (y) + grid.numnodesz * grid.numnodesy * (x))
                        int local_x = local_idx / (NUMNODESY(output_data) * NUMNODESZ(output_data));
                        int rest = local_idx % (NUMNODESY(output_data) * NUMNODESZ(output_data));
                        int local_y = rest / NUMNODESZ(output_data);
                        int local_z = local_idx % NUMNODESZ(output_data);

                        int big_x = local_x + src_x * NUMNODESX(output_data);
                        int big_y = local_y + src_y * NUMNODESY(output_data);
                        int big_z = local_z + src_z * NUMNODESZ(output_data);

                        assert(big_x < NUMNODESX(simdata.write_data));
                        assert(big_y < NUMNODESY(simdata.write_data));
                        assert(big_z < NUMNODESZ(simdata.write_data));

                        assert(j < NUMNODESTOT(simdata.write_data->grid));

                        double data_to_write = simdata.write_data_buffer[j];
                        SETVALUE(simdata.write_data, big_x, big_y, big_z, data_to_write);
                    }
                    write_output(&simdata.params.outputs[i], simdata.write_data, tstep, time);
                    fill_data(simdata.write_data, 0.0);
                };
            }
        }

        if (tstep > 0 && tstep % (numtimesteps / 10) == 0) {
            printf("job %d: step %8d/%d", cart_rank, tstep, numtimesteps);

            if (tstep != numtimesteps) {
                double elapsed_sofar = GET_TIME() - start;
                double timeperstep_sofar = elapsed_sofar / tstep;

                double eta = (numtimesteps - tstep) * timeperstep_sofar;

                printf(" (ETA: %8.3lf seconds)", eta);
            }

            printf("\n");
            fflush(stdout);
        }

        rotate_v_recv_request(tstep, &vx_recv_req, &vy_recv_req, &vz_recv_req, simdata.v_buf_old, simdata.v_buf_new, cart_comm);
        update_pressure(&simdata);
        swap_p_timesteps(&simdata);
        copy_send_p_data_to_buffers(&simdata);
        rotate_p_send_request(tstep, &px_send_req, &py_send_req, &pz_send_req, simdata.p_buf_old, simdata.p_buf_new, cart_comm);

        rotate_p_recv_request(tstep, &px_recv_req, &py_recv_req, &pz_recv_req, simdata.p_buf_old, simdata.p_buf_new, cart_comm);
        update_velocities(&simdata, coords, dims);
        swap_v_timesteps(&simdata);
        copy_send_v_data_to_buffers(&simdata);
        rotate_v_send_request(tstep, &vx_send_req, &vy_send_req, &vz_send_req, simdata.v_buf_old, simdata.v_buf_new, cart_comm);
    }

    DEBUG_PRINT("All good, time to clean up");

    double elapsed = GET_TIME() - start;
    double numupdates =
        (double)NUMNODESTOT(simdata.pold->grid) * (numtimesteps + 1);
    double updatespers = numupdates / elapsed / 1e6;

    printf("\nElapsed %.6lf seconds (%.3lf Mupdates/s)\n\n", elapsed,
           updatespers);

    finalize_simulation(&simdata);
    DEBUG_PRINT("Only thing missing is MPI finalize");
    MPI_Finalize();

    return 0;
}
