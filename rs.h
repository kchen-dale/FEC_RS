#ifndef RS_H
#define RS_H

#include "util.h"

#define K 4 //Number of original packets
#define N 6 //Number of packets after FEC encoding. N - K parity packets are added.
#define PKT_LEN 8 //Number of bytes in each packets

// Vandermonde rows are evaluated at 0..N-1, which must be distinct elements of GF(2^8)
#if K < 1 || K > N || N > 256
#error "Require 1 <= K <= N <= 256"
#endif

Matrix *get_encoded_matrix_after_loss(Matrix *encoded_data, int *loss_index, int loss_count);
Matrix *get_encoding_matrix(Matrix *V);
Matrix *get_decoding_matrix(Matrix *encoding, int *loss_index, int loss_count);
Matrix *rs_encoding(Matrix *encoding, Matrix *data);
Matrix *rs_decoding(Matrix *decoding, Matrix *data);
#endif