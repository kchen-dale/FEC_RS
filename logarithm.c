#include <stdio.h>
#include <stdlib.h>
#include "logarithm.h"

short LOG_TABLE[FIELD_SIZE];
unsigned char EXP_TABLE[FIELD_SIZE * 2 - 2];

//cyclic property insures a^0 = a^255 = 1
short* generateLogTable(int polynomial) {
    int generator = 1; // Starting point for the generator
    polynomial = 0x100 | polynomial; //Add the implicit x^8
    short* result = (short*)malloc(FIELD_SIZE * sizeof(short));
    if (result == NULL) {
        perror("Memory allocation failed");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < FIELD_SIZE; i++) {
        result[i] = -1; // -1 means "not set"
    }

    for (int i = 0; i < FIELD_SIZE-1; i++) {
        if (result[generator] != -1) {
            fprintf(stderr, "BUG: duplicate logarithm (bad polynomial?). Polynomial must be primitive.\n");
            free(result); // Clean up before exiting.
            exit(EXIT_FAILURE);
        }
        result[generator] = i;
        generator <<= 1; // Multiply by x in GF. i.e. move to next x^n in the sequence (1, x, x^2, ..., x^255)

        if (generator & FIELD_SIZE) { // Check for overflow
            generator ^= polynomial; // Reduce modulo the primitive polynomial
        }
        

    }
    return result;
}

unsigned char* generateExpTable(short* logTable) {
    // Allocate memory for the result array
    unsigned char* result = (unsigned char*)malloc((FIELD_SIZE * 2 - 2) * sizeof(unsigned char));
    if (result == NULL) {
        perror("Memory allocation failed");
        exit(EXIT_FAILURE);
    }

    // Populate the result array
    for (int i = 1; i < FIELD_SIZE; i++) {
        int log = logTable[i];
        result[log] = (unsigned char)i;
        result[log + FIELD_SIZE - 1] = (unsigned char)i;
    }

    return result;
}

void initializeTables() {
    short* logTable = generateLogTable(POLYNOMIAL);
    unsigned char* expTable = generateExpTable(logTable); // Assume generateExpTable is defined

    for (int i = 0; i < FIELD_SIZE; i++) {
        LOG_TABLE[i] = logTable[i];
    }

    for (int i = 0; i < FIELD_SIZE * 2 - 2; i++) {
        EXP_TABLE[i] = expTable[i];
    }

    // Free temporary memory if generateLogTable and generateExpTable allocate dynamically
    free(logTable);
    free(expTable);
}
