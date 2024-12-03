#include <stdio.h>
#include "logarithm.h"
#include "util.h"
#include "rs.h"

int main() {
    //Build tables
    initializeTables();
    // Print the log table for verification
    // for (int i = 0; i < FIELD_SIZE; i++) {
    //     printf("logTable[%d] = %d\n", i, LOG_TABLE[i]);
    // }

    // // Print the exp table for verification
    // for (int i = 0; i < FIELD_SIZE - 1; i++) {
    //     printf("expTable[%d] = %d\n", i, EXP_TABLE[i]);
    // }

    // Example: Multiply 0x57 and 0x83 in GF(2^8)
    //Test multiplication using 0x11D polynomial
    // uint8_t a = 0x1; // Example input
    // uint8_t b = 0x83; // Example input
    // uint8_t result = multiply(a, b);

    // printf("Result of 0x%x * 0x%x in GF(2^8) is: 0x%x\n", a, b, result);
    Matrix *matrix = get_vandermonde();
    Matrix *encoding = get_encoding_matrix(matrix);
    int loss_index[] = {0, 3};
    int loss_count = 2;
    Matrix *decoding = get_decoding_matrix(encoding, loss_index, loss_count);
    Matrix *data = read_data_from_file("data.txt", K, PKT_LEN);
    Matrix *after_encoding;
    Matrix *after_decoding;
    Matrix *encoded_with_drop;
    int is_equal = 0;
    // Print the matrix
    printf("=====Vandermonde Matrix=====\n");
    for (int r = 0; r < matrix->rows; r++) {
        for (int c = 0; c < matrix->cols; c++) {
            printf("%02x ", matrix->data[r][c]);
        }
        printf("\n");
    }
    printf("=====Encoding Matrix=====\n");
    for (int r = 0; r < encoding->rows; r++) {
        for (int c = 0; c < encoding->cols; c++) {
            printf("%02x ", encoding->data[r][c]);
        }
        printf("\n");
    }
    printf("=====Data Matrix=====\n");
    for (int r = 0; r < data->rows; r++) {
        for (int c = 0; c < data->cols; c++) {
            printf("%02x ", data->data[r][c]);
        }
        printf("\n");
    }
    printf("=====Encoded Data Matrix=====\n");
    after_encoding = rs_encoding(encoding, data);
    for (int r = 0; r < after_encoding->rows; r++) {
        for (int c = 0; c < after_encoding->cols; c++) {
            printf("%02x ", after_encoding->data[r][c]);
        }
        printf("\n");
    }
    printf("=====Encoded Data Matrix with Drop=====\n");
    encoded_with_drop = get_encoded_matrix_after_loss(after_encoding, loss_index, loss_count);
    for (int r = 0; r < encoded_with_drop->rows; r++) {
        for (int c = 0; c < encoded_with_drop->cols; c++) {
            printf("%02x ", encoded_with_drop->data[r][c]);
        }
        printf("\n");
    }
    printf("=====Decoding Matrix=====\n");
    for (int r = 0; r < decoding->rows; r++) {
        for (int c = 0; c < decoding->cols; c++) {
            printf("%02x ", decoding->data[r][c]);
        }
        printf("\n");
    }
    printf("=====Data Matrix after Decoding=====\n");
    after_decoding = rs_decoding(decoding, encoded_with_drop);
    for (int r = 0; r < after_decoding->rows; r++) {
        for (int c = 0; c < after_decoding->cols; c++) {
            printf("%02x ", after_decoding->data[r][c]);
        }
        printf("\n");
    }

    //Check if data matrix and after_decoding matrix are equal
    is_equal = matrix_equals(data, after_decoding);
    if(is_equal)
        printf("RS recovered the lost packets successfully!");

    // Free the matrix
    free_matrix(matrix);
    free_matrix(encoding);
    free_matrix(decoding);
    free_matrix(data);
    free_matrix(after_encoding);
    free_matrix(encoded_with_drop);
    free_matrix(after_decoding);

    return 0;
}