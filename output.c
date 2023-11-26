#include "output.h"

#include <stdlib.h>

void print_source(source_t* source) {
    printf(" Source infos:\n\n");

    if (source->type == AUDIO) {
        double duration = (double)source->numsamples / source->sampling;

        printf("          type: audio data file\n");
        printf("      sampling: %d Hz\n", source->sampling);
        printf("      duration: %g\n", duration);

    } else {
        printf("          type: sine wave\n");
        printf("     frequency: %g Hz\n", source->data[0]);
    }

    printf("    position x: %g\n", source->posx);
    printf("    position y: %g\n", source->posy);
    printf("    position z: %g\n\n", source->posy);
}

void print_output(output_t* output) {
    switch (output->source) {
        case PRESSURE:
            printf("      pressure: ");
            break;
        case VELOCITYX:
            printf("    velocity X: ");
            break;
        case VELOCITYY:
            printf("    velocity Y: ");
            break;
        case VELOCITYZ:
            printf("    velocity Z: ");
            break;

        default:
            break;
    }

    switch (output->type) {
        case ALL:
            printf("complete dump");
            break;
        case CUTX:
            printf("cut along the x axis at %g", output->posx);
            break;
        case CUTY:
            printf("cut along the y axis at %g", output->posy);
            break;
        case CUTZ:
            printf("cut along the z axis at %g", output->posz);
            break;
        case POINT:
            printf("single point at %g %g %g", output->posx, output->posy,
                   output->posz);
            break;

        default:
            break;
    }

    printf(" to file %s\n", output->filename);
}

int write_output(output_t* output, data_t* data, int step, double time) {
    if (output == NULL || data == NULL) {
        DEBUG_PRINT("NULL pointer passed as argument");
        return 1;
    }

    output_type_t type = output->type;

    if (type == ALL) {
        return write_data(output->fp, data, step, time);
    }

    int m, n, p;
    closest_index(&data->grid, output->posx, output->posy, output->posz, &m, &n,
                  &p);

    int startm = (type == CUTX || type == POINT) ? m : 0;
    int startn = (type == CUTY || type == POINT) ? n : 0;
    int startp = (type == CUTZ || type == POINT) ? p : 0;

    int endm = (type == CUTX || type == POINT) ? m + 1 : NUMNODESX(data);
    int endn = (type == CUTY || type == POINT) ? n + 1 : NUMNODESY(data);
    int endp = (type == CUTZ || type == POINT) ? p + 1 : NUMNODESZ(data);

    data_t* tmpdata = allocate_data(&output->grid);

    for (m = startm; m < endm; m++) {
        for (n = startn; n < endn; n++) {
            for (p = startp; p < endp; p++) {
                int tmpm = m - startm;
                int tmpn = n - startn;
                int tmpp = p - startp;

                SETVALUE(tmpdata, tmpm, tmpn, tmpp, GETVALUE(data, m, n, p));
            }
        }
    }

    int writeok = (write_data(output->fp, tmpdata, step, time) == 0);

    free(tmpdata->vals);
    free(tmpdata);

    if (writeok == 0) {
        DEBUG_PRINT("Failed to write output data");
        return 1;
    }

    return 0;
}

int open_outputfile(output_t* output, grid_t* simgrid) {
    if (output == NULL || simgrid == NULL) {
        DEBUG_PRINT("Invalid NULL pointer in argment");
        return 1;
    }

    grid_t grid;

    output_type_t type = output->type;

    grid.numnodesx = (type == POINT || type == CUTX) ? 1 : simgrid->numnodesx;
    grid.numnodesy = (type == POINT || type == CUTY) ? 1 : simgrid->numnodesy;
    grid.numnodesz = (type == POINT || type == CUTZ) ? 1 : simgrid->numnodesz;

    grid.xmin = (type == POINT || type == CUTX) ? output->posx : simgrid->xmin;
    grid.xmax = (type == POINT || type == CUTX) ? output->posx : simgrid->xmax;

    grid.ymin = (type == POINT || type == CUTY) ? output->posy : simgrid->ymin;
    grid.ymax = (type == POINT || type == CUTY) ? output->posy : simgrid->ymax;

    grid.zmin = (type == POINT || type == CUTZ) ? output->posz : simgrid->zmin;
    grid.zmax = (type == POINT || type == CUTZ) ? output->posz : simgrid->zmax;

    FILE* fp;
    if ((fp = create_datafile(grid, output->filename)) == NULL) {
        DEBUG_PRINTF("Failed to open output file: '%s'", output->filename);
        return 1;
    }

    output->grid = grid;
    output->fp = fp;

    return 0;
}
