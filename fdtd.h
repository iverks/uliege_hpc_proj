#pragma once

#include <mpi.h>
#include <stdio.h>

#define NUM_DIMS 3
// #define NUM_NODES 27
// #define NUM_NODES_PER_DIM 3

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#if NDEBUG

#define DEBUG_PRINTF(fmt, ...)
#define DEBUG_PRINT(msg)

#else

#define DEBUG_PRINTF(fmt, ...)                                           \
    printf("[DEBUG][%s:%d] " fmt "\n", __FILE__, __LINE__, __VA_ARGS__); \
    fflush(stdout);
#define DEBUG_PRINT(msg)                                    \
    printf("[DEBUG][%s:%d] %s\n", __FILE__, __LINE__, msg); \
    fflush(stdout);

#endif

#define MAX_OUTPUTS 32
#define BUFSZ_LARGE 256
#define BUFSZ_HUMONGOUS 512
#define BUFSZ_SMALL 16

#define BUFFMT_LARGE "%255s"
#define BUFFMT_SMALL "%15s"

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

#define NUMNODESX(dat) ((dat)->grid.numnodesx)
#define NUMNODESY(dat) ((dat)->grid.numnodesy)
#define NUMNODESZ(dat) ((dat)->grid.numnodesz)

#define XMIN(dat) ((dat)->grid.xmin)
#define YMIN(dat) ((dat)->grid.ymin)
#define ZMIN(dat) ((dat)->grid.zmin)

#define XMAX(dat) ((dat)->grid.xmax)
#define YMAX(dat) ((dat)->grid.ymax)
#define ZMAX(dat) ((dat)->grid.zmax)

#define NUMNODESTOT(grid) \
    ((size_t)(grid).numnodesx * (grid).numnodesy * (grid).numnodesz)

#define INDEX3D(grid, m, n, p) \
    ((size_t)(p) + grid.numnodesz * (n) + grid.numnodesz * grid.numnodesy * (m))

#define GETVALUE(dat, m, n, p) ((dat)->vals[INDEX3D((dat)->grid, m, n, p)])
#define SETVALUE(dat, m, n, p, val) \
    ((dat)->vals[INDEX3D((dat)->grid, m, n, p)] = (val))

#ifdef _OPENMP

#include <omp.h>

#define GET_TIME() (omp_get_wtime())

#else

#define GET_TIME() ((double)clock() / CLOCKS_PER_SEC)

#endif

typedef enum source_type {
    SINE = 0,
    AUDIO,
    SOURCE_TYPE_END

} source_type_t;

typedef enum output_source {
    PRESSURE = 0,
    VELOCITYX,
    VELOCITYY,
    VELOCITYZ,
    OUTPUT_SOURCE_END

} output_source_t;

typedef enum output_type {
    CUTX = 0,
    CUTY,
    CUTZ,
    ALL,
    POINT,
    OUTPUT_TYPE_END

} output_type_t;

typedef enum buffer_direction {
    X_MIN = 0,
    Y_MIN,
    Z_MIN,
    X_MAX,
    Y_MAX,
    Z_MAX,
    BUFFER_DIR_TYPE_END
} buffer_direction_t;

typedef struct grid {
    int numnodesx;
    int numnodesy;
    int numnodesz;

    double xmin;
    double ymin;
    double zmin;

    double xmax;
    double ymax;
    double zmax;

} grid_t;

typedef struct source {
    source_type_t type;

    double posx;
    double posy;
    double posz;

    int sampling;
    int numsamples;

    double* data;

} source_t;

typedef struct output {
    output_type_t type;
    output_source_t source;

    char* filename;

    double posx;
    double posy;
    double posz;

    grid_t grid;

    FILE* fp;

} output_t;

typedef struct parameters {
    double dx;
    double dt;
    double maxt;

    int outrate;
    int numoutputs;

    char* cin_filename;
    char* rhoin_filename;

    source_t source;
    output_t* outputs;

} parameters_t;

typedef struct data {
    grid_t grid;

    double* vals;
} data_t;

typedef struct buffer {
    grid_t grid;

    double** buffers;
} buffer_t;

typedef struct simulation_data {
    parameters_t params;

    data_t *c, *rho, *rhohalf;

    data_t *pold, *pnew;
    data_t *vxold, *vxnew;
    data_t *vyold, *vynew;
    data_t *vzold, *vznew;

    buffer_t *p_buf_old, *p_buf_new;
    buffer_t *v_buf_old, *v_buf_new;
} simulation_data_t;

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

/******************************************************************************
 * Output file functions                                                      *
 ******************************************************************************/

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

/******************************************************************************
 * Parameters file functions                                                  *
 ******************************************************************************/

/**
 * @brief Read an audio source from file
 *
 * @param filename [IN] the file to read the audio from
 * @param source [INOUT] the source object in which to store the data read
 * from the audio file
 * @return int 0 if read was a success, returns 1 otherwise
 */
int read_audiosource(char* filename, source_t* source);

/**
 * @brief Read an output specification parameter from a parameter file
 *
 * @param fp [IN] a file handle to an opened parameter file
 * @param output [OUT] an output oject that will store the output specification
 * read from the parameter file
 * @return int 0 if read was a success, returns 1 otherwise
 */
int read_outputparam(FILE* fp, output_t* output, int coords[]);

/**
 * @brief Read a source specification parameter from a parameter file
 *
 * @param fp [IN] a file handle to an opened parameter file
 * @param source [OUT] a source oject that will store the source specification
 * read from the parameter file
 * @return int 0 if read was a success, returns 1 otherwise
 */
int read_sourceparam(FILE* fp, source_t* source);

/**
 * @brief Read the parameter file specified in argument
 *
 * @param params [OUT] a parameter object into which the parameter read from
 * the file will be stored
 * @param filename [IN] path to the parameter file to read
 * @return int 0 if read was a success, returns 1 otherwise
 */
int read_paramfile(parameters_t* params, const char* filename, int coords[]);

/******************************************************************************
 * Simulation-related functions                                               *
 ******************************************************************************/

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
 * @brief Interpolate the input density and spped maps so that they corresponds
 * to the simulation grid
 *
 * @param simdata [OUT] a simulation data object used to store the interpolated
 * values that will be used during the simulation
 * @param simgrid [IN] the simulation grid, i.e., the grid to which the input
 * data will be converted to
 * @param cin [IN] the input speed map
 * @param rhoin [IN] the input density map
 * @return int 0 if read was a success, returns 1 otherwise
 */
int interpolate_inputmaps(simulation_data_t* simdata, grid_t* simgrid,
                          data_t* cin, data_t* rhoin);

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
