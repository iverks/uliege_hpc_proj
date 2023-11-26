#include "parameters.h"

#include <stdlib.h>
#include <string.h>

#include "utilities.h"

const char* source_type_keywords[] = {[SINE] = "sine", [AUDIO] = "audio"};

const char* output_type_keywords[] = {[CUTX] = "cut_x",
                                      [CUTY] = "cut_y",
                                      [CUTZ] = "cut_z",
                                      [ALL] = "all",
                                      [POINT] = "point"};

const char* output_source_keywords[] = {[PRESSURE] = "pressure",
                                        [VELOCITYX] = "velocity_x",
                                        [VELOCITYY] = "velocity_y",
                                        [VELOCITYZ] = "velocity_z"};

int read_audiosource(char* filename, source_t* source) {
    FILE* fp;
    if ((fp = fopen(filename, "rb")) == NULL) {
        DEBUG_PRINTF("Could not open source file '%s'", filename);
        return 1;
    }

    fseek(fp, 0, SEEK_END);
    size_t filesize = ftell(fp);
    rewind(fp);

    int numsamples = (filesize - sizeof(int)) / sizeof(double);

    int sampling;
    if (fread(&sampling, sizeof(int), 1, fp) != 1) {
        DEBUG_PRINT("Failed to read source data");
        fclose(fp);
        return 1;
    }

    double* data;
    if ((data = malloc(numsamples * sizeof(double))) == NULL) {
        DEBUG_PRINT("Failed to allocate memory for source data");
        return 1;
    }

    int readok = (fread(data, sizeof(double), numsamples, fp) == (long unsigned int)numsamples);

    fclose(fp);

    if (readok == 0) {
        DEBUG_PRINT("Failed to read source data");
        return 1;
    }

    source->data = data;
    source->numsamples = numsamples;
    source->sampling = sampling;

    return 0;
}

int read_outputparam(FILE* fp, output_t* output, int coords[]) {
    if (fp == NULL || output == NULL) {
        DEBUG_PRINT("NULL passed as argement");
        return 1;
    }

    char typekeyword[BUFSZ_SMALL];
    char sourcekeyword[BUFSZ_SMALL];
    char filename[BUFSZ_LARGE];
    char formatted_filename[BUFSZ_HUMONGOUS];

    double posxyz[3] = {0.0, 0.0, 0.0};

    if (fscanf(fp, BUFFMT_SMALL, typekeyword) != 1 ||
        fscanf(fp, BUFFMT_SMALL, sourcekeyword) != 1 ||
        fscanf(fp, BUFFMT_LARGE, filename) != 1) {
        DEBUG_PRINT("Failed to read an output parameter");
        return 1;
    }

    output_type_t type = CUTX;
    while (type < OUTPUT_TYPE_END &&
           strcmp(output_type_keywords[type], typekeyword) != 0) {
        type++;
    }

    if (type == OUTPUT_TYPE_END) {
        DEBUG_PRINTF("Invalid keyword: '%s'", typekeyword);
        return 1;
    }

    output_source_t source = PRESSURE;
    while (source < OUTPUT_SOURCE_END &&
           strcmp(output_source_keywords[source], sourcekeyword) != 0) {
        source++;
    }

    if (source == OUTPUT_SOURCE_END) {
        DEBUG_PRINTF("Invalid keyword: '%s'", sourcekeyword);
        return 1;
    }

    int readok = 1;
    switch (type) {
        case CUTX:
            readok = (fscanf(fp, "%lf", &posxyz[0]) == 1);
            break;
        case CUTY:
            readok = (fscanf(fp, "%lf", &posxyz[1]) == 1);
            break;
        case CUTZ:
            readok = (fscanf(fp, "%lf", &posxyz[2]) == 1);
            break;
        case ALL:
            break;

        case POINT:
            readok =
                (fscanf(fp, "%lf %lf %lf", &posxyz[0], &posxyz[1], &posxyz[2]) == 3);
            break;

        default:
            break;
    }

    if (readok == 0) {
        DEBUG_PRINT("Failed to read an output parameter");
        return 1;
    }

    sprintf(formatted_filename, "%d_%d_%d_%s", coords[0], coords[1], coords[2], filename);
    output->filename = copy_string(formatted_filename);
    output->type = type;
    output->source = source;
    output->posx = posxyz[0];
    output->posy = posxyz[1];
    output->posz = posxyz[2];

    return 0;
}

int read_sourceparam(FILE* fp, source_t* source) {
    char typekeyword[BUFSZ_SMALL];
    char filename[BUFSZ_LARGE];

    double freq, posx, posy, posz;

    if (fscanf(fp, BUFFMT_SMALL, typekeyword) != 1) {
        DEBUG_PRINT("Failed to read the source parameter");
        return 1;
    }

    source_type_t type = SINE;
    while (type < SOURCE_TYPE_END &&
           strcmp(source_type_keywords[type], typekeyword) != 0) {
        type++;
    }

    if (type == SOURCE_TYPE_END) {
        DEBUG_PRINTF("Invalid keyword: '%s'", typekeyword);
        return 1;
    }

    int readok = 1;
    switch (type) {
        case SINE:
            readok = (fscanf(fp, "%lf", &freq) == 1);
            break;
        case AUDIO:
            readok = (fscanf(fp, BUFFMT_LARGE, filename) == 1);
            break;

        default:
            break;
    }

    if (readok == 0 || fscanf(fp, "%lf %lf %lf", &posx, &posy, &posz) != 3) {
        DEBUG_PRINT("Failed to read the source parameter");
        return 1;
    }

    switch (type) {
        case AUDIO:
            read_audiosource(filename, source);
            break;
        case SINE: {
            if ((source->data = malloc(sizeof(double))) == NULL) {
                DEBUG_PRINT("Failed to allocate memory");
                return 1;
            }

            source->data[0] = freq;
            source->numsamples = 1;

            break;
        }

        default:
            break;
    }

    source->type = type;
    source->posx = posx;
    source->posy = posy;
    source->posz = posz;

    return 0;
}

int read_paramfile(parameters_t* params, const char* filename, int coords[]) {
    if (params == NULL || filename == NULL) {
        DEBUG_PRINT("Invalid print_out params or filename");
        return 1;
    }

    int outrate, numoutputs = 0;

    double dx, dt, maxt;

    char cin_filename[BUFSZ_LARGE];
    char rhoin_filename[BUFSZ_LARGE];

    source_t source;
    output_t* outputs = NULL;

    if ((outputs = malloc(sizeof(output_t) * MAX_OUTPUTS)) == NULL) {
        DEBUG_PRINT("Failed to allocate memory");
        return 1;
    }

    FILE* fp;
    if ((fp = fopen(filename, "r")) == NULL) {
        DEBUG_PRINTF("Could not open parameter file '%s'", filename);
        return 1;
    }

    int readok =
        ((fscanf(fp, "%lf", &dx) == 1) && (fscanf(fp, "%lf", &dt) == 1) &&
         (fscanf(fp, "%lf", &maxt) == 1) && (fscanf(fp, "%d", &outrate) == 1) &&
         (fscanf(fp, BUFFMT_LARGE, cin_filename) == 1) &&
         (fscanf(fp, BUFFMT_LARGE, rhoin_filename) == 1));

    readok = (readok != 0 && read_sourceparam(fp, &source) == 0 &&
              fscanf(fp, " ") == 0);

    while (readok != 0 && numoutputs < MAX_OUTPUTS && feof(fp) == 0) {
        readok = (read_outputparam(fp, &outputs[numoutputs++], coords) == 0 &&
                  fscanf(fp, " ") == 0);
    }

    fclose(fp);

    if (readok == 0) {
        DEBUG_PRINT("Failed to read parameter file");
        free(outputs);
        return 1;
    }

    if (numoutputs == 0) {
        free(outputs);
        outputs = NULL;

    } else if ((outputs = realloc(outputs, sizeof(output_t) * numoutputs)) ==
               NULL) {
        DEBUG_PRINT("Failed to allocate memory");
        return 1;
    }

    params->dx = dx;
    params->dt = dt;
    params->maxt = maxt;
    params->outrate = outrate;
    params->cin_filename = copy_string(cin_filename);
    params->rhoin_filename = copy_string(rhoin_filename);
    params->source = source;
    params->numoutputs = numoutputs;
    params->outputs = outputs;

    return 0;
}

FILE* open_datafile(grid_t* grid, int* numsteps, char* filename) {
    if (grid == NULL || filename == NULL) {
        DEBUG_PRINT("Invalid NULL grid or filename");
        return NULL;
    }

    FILE* fp;
    if ((fp = fopen(filename, "rb")) == NULL) {
        DEBUG_PRINTF("Failed to open file '%s'", filename);
        return NULL;
    }

    fseek(fp, 0, SEEK_END);
    size_t file_size = ftell(fp);
    rewind(fp);

    if (fread(&grid->numnodesx, sizeof(int), 1, fp) != 1 ||
        fread(&grid->numnodesy, sizeof(int), 1, fp) != 1 ||
        fread(&grid->numnodesz, sizeof(int), 1, fp) != 1 ||
        fread(&grid->xmin, sizeof(double), 1, fp) != 1 ||
        fread(&grid->xmax, sizeof(double), 1, fp) != 1 ||
        fread(&grid->ymin, sizeof(double), 1, fp) != 1 ||
        fread(&grid->ymax, sizeof(double), 1, fp) != 1 ||
        fread(&grid->zmin, sizeof(double), 1, fp) != 1 ||
        fread(&grid->zmax, sizeof(double), 1, fp) != 1) {
        DEBUG_PRINTF("Failed to read header of file '%s'", filename);
        fclose(fp);
        return NULL;
    }

    size_t numnodestot =
        (size_t)grid->numnodesx * grid->numnodesy * grid->numnodesz;

    size_t values_size = numnodestot * sizeof(double);
    size_t stepindex_size = sizeof(int);
    size_t timestamp_size = sizeof(double);
    size_t header_size = 6 * sizeof(double) + 3 * sizeof(int);

    size_t onetimestep_size = values_size + stepindex_size + timestamp_size;
    size_t alltimestep_size = file_size - header_size;

    if (alltimestep_size % onetimestep_size != 0) {
        DEBUG_PRINTF("Data size is inconsistent with number of nodes (%lu, %lu)",
                     alltimestep_size, onetimestep_size);

        fclose(fp);
        return NULL;
    }

    if (numsteps != NULL) {
        *numsteps = (alltimestep_size / onetimestep_size);
    }

    return fp;
}

data_t* read_data(FILE* fp, grid_t* grid, int* step, double* time) {
    if (fp == NULL) {
        DEBUG_PRINT("Invalid NULL file pointer");
        return NULL;
    }

    double ltime;
    int lstep;

    size_t numnodes = NUMNODESTOT(*grid);

    data_t* data;
    if ((data = allocate_data(grid)) == NULL) {
        DEBUG_PRINT("Failed to allocate data");
        return NULL;
    }

    if (fread(&lstep, sizeof(int), 1, fp) != 1 ||
        fread(&ltime, sizeof(double), 1, fp) != 1 ||
        fread(data->vals, sizeof(double), numnodes, fp) != numnodes) {
        DEBUG_PRINT("Failed to read data");
        free(data);
        return NULL;
    }

    if (step != NULL)
        *step = lstep;
    if (time != NULL)
        *time = ltime;

    return data;
}