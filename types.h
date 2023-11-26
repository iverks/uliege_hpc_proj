#pragma once
#include <stdlib.h>

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

typedef enum buffer_direction {
    X_MIN = 0,
    Y_MIN,
    Z_MIN,
    X_MAX,
    Y_MAX,
    Z_MAX,
    BUFFER_DIR_TYPE_END
} buffer_direction_t;

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

    data_t* write_data;         // Structured
    double* write_data_buffer;  // Unstructured

    buffer_t *p_buf_old, *p_buf_new;
    buffer_t *v_buf_old, *v_buf_new;
} simulation_data_t;