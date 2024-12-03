#include <stdio.h>
#include "logarithm.h"
#include "util.h"
#include "rs.h"

// Drop specific rows according to loss_index. Always return a matrix with K rows.
Matrix *get_matrix_with_drop(Matrix *encoding, int *loss_index, int loss_count) {
    int new_row = 0;
    // Check if loss_count exceeds N - K
    if (loss_count > N - K) {
        printf("Error: Cannot recover more than %d loss packets.\n", N - K);
        return NULL;
    }
    Matrix *loss_matrix = create_matrix(K, encoding->cols);
    for (int r = 0; r < encoding->rows; r++) {
        // Check if the current row index is in the remove_indices array
        int should_remove = 0;
        for (int i = 0; i < loss_count; i++) {
            if (r == loss_index[i]) {
                should_remove = 1;
                break;
            }
        }
        // If not in remove_indices, copy the row to the new matrix
        if (!should_remove) {
            for (int c = 0; c < encoding->cols; c++) {
                loss_matrix->data[new_row][c] = encoding->data[r][c];
            }
            new_row++;
        }
        if(new_row == K) break;
    }
    return loss_matrix;
}

//Create the encoding matrix given a Vandermonde matrix
Matrix *get_encoding_matrix(Matrix *V) {
    Matrix *result;
    Matrix *upper = submatrix(V, 0, K-1, 0, K-1);
    Matrix *upper_inverse = invert(upper);
    free_matrix(upper);
    result = matrix_multiply(V, upper_inverse);
    free_matrix(upper_inverse);

    return result;
}

//Create the decoding matrix based on the loss position
Matrix *get_decoding_matrix(Matrix *encoding, int *loss_index, int loss_count) {
    
    Matrix *result;
    Matrix *loss_matrix = get_matrix_with_drop(encoding, loss_index, loss_count);
    result = invert(loss_matrix);
    free_matrix(loss_matrix);

    return result;
}

Matrix *rs_encoding(Matrix *encoding, Matrix *data) {
    return matrix_multiply(encoding, data);
}

Matrix *rs_decoding(Matrix *decoding, Matrix *data) {
    return matrix_multiply(decoding, data);
}

Matrix *get_encoded_matrix_after_loss(Matrix *encoded_data, int *loss_index, int loss_count) {
    return get_matrix_with_drop(encoded_data, loss_index, loss_count);
}
