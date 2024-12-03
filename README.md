# FEC_RS
This repository provides an implementation of Forward Error Correction (FEC) using Reed-Solomon (RS) codes.

## Features
- **Encoding and Decoding**: Core functions to encode data into RS codewords and decode received codewords, correcting errors and erasures.
- **Galois Field Arithmetic**: Implementation of arithmetic operations over Galois Fields, essential for RS coding.
- **Utility Functions**: Additional helper functions to facilitate encoding and decoding processes.

## File Structure

The repository contains the following key files:

- **`logarithm.c`** and **`logarithm.h`**: Logarithm and anti-logarithm tables for efficient computation in Galois Fields.
- **`rs.c`** and **`rs.h`**: Implementation of RS encoding and decoding algorithms.
- **`util.c`** and **`util.h`**: Utility functions supporting the main RS operations.
- **`test.c`**: Test cases to validate the functionality of the RS code implementation.
- **`data.txt`**: Sample data used in test cases.
- **`Makefile`**: Build instructions for compiling the project.
- **`galois.c`**: Another implementation of GF multiplication using polynomail multiplication. Only for comparison.

## Getting Started

To build and test the project, follow these steps:

   ```bash
   git clone https://github.com/kchen-dale/FEC_RS.git
   cd FEC_RS
   make
   ./program
