#include "utilities.h"

#include <stdlib.h>
#include <string.h>

char* copy_string(char* str) {
    size_t len;
    if (str == NULL || (len = strlen(str)) == 0) {
        DEBUG_PRINT("NULL of zero length string passed as argument");
        return NULL;
    }

    char* cpy;
    if ((cpy = malloc((len + 1) * sizeof(char))) == NULL) {
        DEBUG_PRINT("Failed to allocate memory");
        return NULL;
    }

    return strcpy(cpy, str);
}