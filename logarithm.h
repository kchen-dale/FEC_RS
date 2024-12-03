#ifndef LOGARITHM_H
#define LOGARITHM_H

#include <stdint.h>

#define FIELD_SIZE 256
// Possible polynomails for GF(2^8): 29, 43, 45, 77, 95, 99, 101, 105,
// 113, 135, 141, 169, 195, 207, 231, and 245.
// Use different polynomials will result in differnet multiplication results.
#define POLYNOMIAL 0x1D

// Declare the global tables as external
extern short LOG_TABLE[FIELD_SIZE];
extern unsigned char EXP_TABLE[FIELD_SIZE * 2 - 2];

// Declare the initialization function
void initializeTables();

#endif // LOGARITHM_H
