#ifndef RS_H
#define RS_H

#include "util.h"

#define K 4 //Number of original packets
#define N 6 //Number of packets after FEC encoding. K - N parity packets are added.
#define PKT_LEN 8 //Number of bytes in each packets

Matrix *get_encoded_matrix_after_loss(Matrix *encoded_data, int *loss_index, int loss_count);
Matrix *get_encoding_matrix(Matrix *V);
Matrix *get_decoding_matrix(Matrix *encoding, int *loss_index, int loss_count);
Matrix *rs_encoding(Matrix *encoding, Matrix *data);
Matrix *rs_decoding(Matrix *decoding, Matrix *data);
#endif