#ifndef UTIL_H
#define UTIL_H

// Define a matrix structure
typedef struct {
    int rows;
    int cols;
    unsigned char **data;
} Matrix;

unsigned char multiply(unsigned char a, unsigned char b);
unsigned char galois_exp(unsigned char a, int n);
int matrix_equals(Matrix *a, Matrix *b);

Matrix *create_matrix(int rows, int cols);
Matrix *get_vandermonde(void);

Matrix *matrix_multiply(Matrix *left, Matrix *right);
Matrix *submatrix(Matrix *matrix, int rmin, int rmax, int cmin, int cmax);
Matrix *invert(Matrix *matrix);
void free_matrix(Matrix *matrix);

Matrix *read_data_from_file(const char *filename, int rows, int cols);

#endif //UTIL_H