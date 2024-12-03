#include <stdio.h>
#include <stdint.h>

// We need primitive polynomial.
// All primitive polynomials are irreducible, because we simply defined them that way: to be primitive, 
// a polynomial must first be irreducible.
// However, not all irreducible polynomials are primitive.
// Function to multiply two elements in GF(2^8)
uint8_t gf_multiply(uint8_t a, uint8_t b) {
    uint8_t result = 0;   // Initialize the result
    uint8_t mod_poly = 0x1D; // x^8 + x^4 + x^3 + x^2 + 1 (0x11D without the implicit x^8)
    // uint8_t mod_poly = 0b1011;
    while (b) {
        // If the lowest bit of b is set, add (XOR) the current value of a to the result
        if (b & 1) {
            result ^= a;
        }
        
        // Prepare for the next bit of b
        // After each left-shift, check if a exceeds 8 bits (degree ≥8): If so, XOR a with
        // P(x) to reduce
        if (a & 0x80) { // Check if x^7 term in a (highest bit) is set
        // if (a & 0b100) { 
            a = (a << 1) ^ mod_poly;
        } else {
            a <<= 1;
        }
        
        // Move to the next bit of b
        b >>= 1;
    }
    
    return result;
}

int main() {
    // Example: Multiply 0x57 and 0x83 in GF(2^8)
    uint8_t a = 0x57; // Example input
    uint8_t b = 0x83; // Example input
    // uint8_t a = 0b110; // Example input
    // uint8_t b = 0b011; // Example input
    uint8_t result = gf_multiply(a, b);

    printf("Result of 0x%x * 0x%x in GF(2^8) is: 0x%x\n", a, b, result);
    return 0;
}
