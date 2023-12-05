#include <assert.h>
#include <execinfo.h>
#include <math.h>
#include <signal.h>
#include <unistd.h>

#include "data.h"
#include "output.h"
#include "parameters.h"
#include "simulation.h"
#include "types.h"
#include "utilities.h"

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

    simulation_data_t simdata;
    init_simulation(&simdata, argv[1]);

    int numtimesteps = floor(simdata.params.maxt / simdata.params.dt);

    double start = GET_TIME();
    // Prepare to recieve already in loop idx 0

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

                write_output(&simdata.params.outputs[i], output_data, tstep, time);
            };
        }

        if (tstep > 0 && tstep % (numtimesteps / 10) == 0) {
            printf("  step %8d/%d", tstep, numtimesteps);

            if (tstep != numtimesteps) {
                double elapsed_sofar = GET_TIME() - start;
                double timeperstep_sofar = elapsed_sofar / tstep;

                double eta = (numtimesteps - tstep) * timeperstep_sofar;

                printf(" (ETA: %8.3lf seconds)", eta);
            }

            printf("\n");
            fflush(stdout);
        }

        update_pressure(&simdata);
        swap_p_timesteps(&simdata);

        update_velocities(&simdata);
        swap_v_timesteps(&simdata);
    }

    DEBUG_PRINT("All good, time to clean up");

    double elapsed = GET_TIME() - start;
    double numupdates =
        (double)NUMNODESTOT(simdata.pold->grid) * (numtimesteps + 1);
    double updatespers = numupdates / elapsed / 1e6;

    printf("\nElapsed %.6lf seconds (%.3lf Mupdates/s)\n\n", elapsed,
           updatespers);

    finalize_simulation(&simdata);

    return 0;
}
