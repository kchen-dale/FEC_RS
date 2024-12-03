#include <stdio.h>
#include <stdlib.h>
#include "logarithm.h"
#include "util.h"
#include "rs.h"


//Adds two elements of the field (XOR operation).
unsigned char add(unsigned char a, unsigned char b) {
    return a ^ b;
}

 //Subtracts two elements of the field (XOR operation).
unsigned char subtract(unsigned char a, unsigned char b) {
    return a ^ b;
}

// Multiplies two elements of the field.
unsigned char multiply(unsigned char a, unsigned char b) {
    if (a == 0 || b == 0) {
        return 0;
    } else {
        int logA = LOG_TABLE[a];
        int logB = LOG_TABLE[b];
        int logResult = logA + logB;
        return EXP_TABLE[logResult];
    }
}

//Divides two elements of the field.
unsigned char divide(unsigned char a, unsigned char b) {
    if (a == 0) {
        return 0;
    }
    if (b == 0) {
        fprintf(stderr, "Error: Argument 'divisor' is 0\n");
        exit(EXIT_FAILURE);
    }
    int logA = LOG_TABLE[a];
    int logB = LOG_TABLE[b];
    int logResult = logA - logB;
    if (logResult < 0) {
        logResult += FIELD_SIZE - 1; // Adjust for field wrap-around
    }
    return EXP_TABLE[logResult];
}

// Computes a**n.
unsigned char galois_exp(unsigned char a, int n) {
    if (n == 0) {
        return 1;
    } else if (a == 0) {
        return 0;
    } else {
        int logA = LOG_TABLE[a];
        int logResult = logA * n;
        while (logResult >= FIELD_SIZE - 1) {
            logResult -= FIELD_SIZE - 1;
        }
        return EXP_TABLE[logResult];
    }
}

// Function to create a new matrix
Matrix *create_matrix(int rows, int cols) {
    // Allocate memory for the matrix structure
    Matrix *matrix = (Matrix *)malloc(sizeof(Matrix));
    if (!matrix) {
        fprintf(stderr, "Memory allocation failed for Matrix structure.\n");
        exit(EXIT_FAILURE);
    }

    matrix->rows = rows;
    matrix->cols = cols;

    // Allocate memory for the row pointers
    matrix->data = (unsigned char **)malloc(rows * sizeof(unsigned char *));
    if (!matrix->data) {
        fprintf(stderr, "Memory allocation failed for Matrix rows.\n");
        free(matrix);
        exit(EXIT_FAILURE);
    }

    // Allocate memory for each row and initialize to 0
    for (int i = 0; i < rows; i++) {
        matrix->data[i] = (unsigned char *)calloc(cols, sizeof(unsigned char)); // Use calloc to initialize to 0
        if (!matrix->data[i]) {
            fprintf(stderr, "Memory allocation failed for Matrix row %d.\n", i);
            // Free previously allocated rows
            for (int j = 0; j < i; j++) {
                free(matrix->data[j]);
            }
            free(matrix->data);
            free(matrix);
            exit(EXIT_FAILURE);
        }
    }

    return matrix;
}
// Function to set a value in the matrix
void set_matrix_value(Matrix *matrix, int row, int col, unsigned char value) {
    matrix->data[row][col] = value;
}

// Function to free the matrix
void free_matrix(Matrix *matrix) {
    for (int i = 0; i < matrix->rows; i++) {
        free(matrix->data[i]);
    }
    free(matrix->data);
    free(matrix);
}

// Vandermonde matrix function
Matrix *vandermonde(int rows, int cols) {
    Matrix *result = create_matrix(rows, cols);
    for (int r = 0; r < rows; r++) {
        for (int c = 0; c < cols; c++) {
            set_matrix_value(result, r, c, galois_exp((unsigned char)r, c));
        }
    }
    return result;
}

//Create Vandermonde matrix based on config
Matrix *get_vandermonde(void) {
    return vandermonde(N, K);
}

//Create an indentity matrix of input size
Matrix *identity(int size) {
    Matrix *result = create_matrix(size, size);
    for (int i = 0; i < size; i++) {
        set_matrix_value(result, i, i, 1); // Set diagonal elements to 1
    }
    return result;
}

//check if two matrix are equal
int matrix_equals(Matrix *a, Matrix *b) {
    if (a->rows != b->rows || a->cols != b->cols) {
        return 0; // Matrices are not of the same dimensions
    }
    for (int r = 0; r < a->rows; r++) {
        for (int c = 0; c < a->cols; c++) {
            if (a->data[r][c] != b->data[r][c]) {
                return 0; // Matrices are not equal
            }
        }
    }
    return 1; // Matrices are equal
}

//Matrix multiplication in GF
Matrix *matrix_multiply(Matrix *left, Matrix *right) {
    if (left->cols != right->rows) {
        fprintf(stderr, "Matrix multiplication error: incompatible dimensions.\n");
        exit(EXIT_FAILURE); // Exit on error
    }

    Matrix *result = create_matrix(left->rows, right->cols);

    for (int r = 0; r < left->rows; r++) {
        for (int c = 0; c < right->cols; c++) {
            unsigned char value = 0;
            for (int i = 0; i < left->cols; i++) {
                value ^= multiply(left->data[r][i], right->data[i][c]);
            }
            set_matrix_value(result, r, c, value);
        }
    }

    return result;
}

//Get the sub matrix of a given matrix. Result is matrix[rmin:rmax][cmin:cmax]
Matrix *submatrix(Matrix *matrix, int rmin, int rmax, int cmin, int cmax) {
    int rows = rmax - rmin + 1;
    int cols = cmax - cmin + 1;

    Matrix *result = create_matrix(rows, cols);
    for (int r = rmin; r <= rmax; r++) {
        for (int c = cmin; c <= cmax; c++) {
            set_matrix_value(result, r - rmin, c - cmin, matrix->data[r][c]);
        }
    }

    return result;
}

//swap two rows of a given matrix
void swap_rows(Matrix *matrix, int r1, int r2) {
    if (r1 < 0 || r1 >= matrix->rows || r2 < 0 || r2 >= matrix->rows) {
        fprintf(stderr, "Row index out of range\n");
        exit(EXIT_FAILURE);
    }
    unsigned char *tmp = matrix->data[r1];
    matrix->data[r1] = matrix->data[r2];
    matrix->data[r2] = tmp;
}

//Concatenate the right matrix to the right of given matrix
Matrix *augment(Matrix *left, Matrix *right) {
    // Check if the matrices have the same number of rows
    if (left->rows != right->rows) {
        fprintf(stderr, "Matrices must have the same number of rows to augment\n");
        exit(EXIT_FAILURE);
    }

    // Create a new matrix with the combined width of the two matrices
    int rows = left->rows;
    int cols = left->cols + right->cols;
    Matrix *result = create_matrix(rows, cols);

    // Copy data from the left matrix into the result
    for (int r = 0; r < rows; r++) {
        for (int c = 0; c < left->cols; c++) {
            result->data[r][c] = left->data[r][c];
        }
    }

    // Copy data from the right matrix into the result
    for (int r = 0; r < rows; r++) {
        for (int c = 0; c < right->cols; c++) {
            result->data[r][left->cols + c] = right->data[r][c];
        }
    }

    return result;
}

void gaussian_elimination(Matrix *matrix) {
    int rows = matrix->rows;
    int cols = matrix->cols;

    // Clear below the main diagonal and scale the diagonal to 1
    for (int r = 0; r < rows; r++) {
        // Ensure diagonal is non-zero by swapping rows if necessary
        if (matrix->data[r][r] == 0) {
            int swapped = 0;
            for (int rowBelow = r + 1; rowBelow < rows; rowBelow++) {
                if (matrix->data[rowBelow][r] != 0) {
                    swap_rows(matrix, r, rowBelow);
                    swapped = 1;
                    break;
                }
            }
            if (!swapped) {
                fprintf(stderr, "Matrix is singular\n");
                exit(EXIT_FAILURE);
            }
        }

        // Scale the diagonal to 1
        if (matrix->data[r][r] != 1) {
            unsigned char scale = divide(1, matrix->data[r][r]);
            for (int c = 0; c < cols; c++) {
                matrix->data[r][c] = multiply(matrix->data[r][c], scale);
            }
        }

        // Clear everything below the diagonal
        for (int rowBelow = r + 1; rowBelow < rows; rowBelow++) {
            if (matrix->data[rowBelow][r] != 0) {
                unsigned char scale = matrix->data[rowBelow][r];
                for (int c = 0; c < cols; c++) {
                    matrix->data[rowBelow][c] ^= multiply(scale, matrix->data[r][c]);
                }
            }
        }
    }

    // Clear above the main diagonal
    for (int d = 0; d < rows; d++) {
        for (int rowAbove = 0; rowAbove < d; rowAbove++) {
            if (matrix->data[rowAbove][d] != 0) {
                unsigned char scale = matrix->data[rowAbove][d];
                for (int c = 0; c < cols; c++) {
                    matrix->data[rowAbove][c] ^= multiply(scale, matrix->data[d][c]);
                }
            }
        }
    }
}

Matrix *invert(Matrix *matrix) {
    if (matrix->rows != matrix->cols) {
        fprintf(stderr, "Only square matrices can be inverted\n");
        exit(EXIT_FAILURE);
    }

    // Create a working matrix by augmenting with an identity matrix
    Matrix *identity_matrix = identity(matrix->rows);
    Matrix *work = augment(matrix, identity_matrix);
    free_matrix(identity_matrix);

    // Perform Gaussian elimination
    gaussian_elimination(work);

    // Extract the right half as the inverse
    Matrix *inverse = submatrix(work, 0, matrix->rows - 1, matrix->cols, matrix->cols * 2 - 1);
    free_matrix(work);

    return inverse;
}

Matrix *read_data_from_file(const char *filename, int rows, int cols) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Error opening file");
        return NULL;
    }

    // Create the matrix
    Matrix *matrix = create_matrix(rows, cols);

    // Read the values into the matrix
    for (int r = 0; r < rows; r++) {
        for (int c = 0; c < cols; c++) {
            int value;
            if (fscanf(file, "%x", &value) != 1 || value < 0 || value > 255) {
                fprintf(stderr, "Error: Invalid value in file at row %d, col %d\n", r, c);
                free_matrix(matrix);
                fclose(file);
                return NULL;
            }
            matrix->data[r][c] = (unsigned char)value;
        }
    }

    fclose(file);
    return matrix;
}