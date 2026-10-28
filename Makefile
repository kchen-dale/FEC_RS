# Compiler and flags
CC = gcc
CFLAGS = -Wall -Wextra -g
DEPFLAGS = -MMD -MP

# Target executable
TARGET = program

# Source files
SRCS = test.c logarithm.c util.c rs.c

# Object files
OBJS = $(SRCS:.c=.o)

# Dependency files
DEPS = $(OBJS:.o=.d)

# Default target
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS)

# Include generated dependency files
-include $(DEPS)

# Compile source files into object files with dependency generation
%.o: %.c
	$(CC) $(CFLAGS) $(DEPFLAGS) -c $< -o $@

# Build and run the demo
run: $(TARGET)
	./$(TARGET)

# Standalone bitwise GF(2^8) multiplication demo (for comparison only)
galois: galois.c
	$(CC) $(CFLAGS) -o $@ $<

# Clean target
clean:
	rm -rf $(OBJS) $(DEPS) $(TARGET) galois *.dSYM

# PHONY targets
.PHONY: run clean
