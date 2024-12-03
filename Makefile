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

# Clean target
clean:
	rm -f $(OBJS) $(DEPS) $(TARGET)

# PHONY targets
.PHONY: clean
