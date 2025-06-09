# Makefile for the K&R Astrology Program

# Compiler and flags
CC = cc
CFLAGS = -O      # Basic optimization, K&R C often doesn't have many CFLAGS
                 # For debugging, you might use: CFLAGS = -g
LDFLAGS =        # Linker flags
LIBS = -lm       # Math library for ephemeris.c

# Source files
SRCS = main_app.c aphorism_utils.c ephemeris.c

# Object files (derived from SRCS)
OBJS = $(SRCS:.c=.o)

# Executable name
TARGET = my_astrology_app

# Default target: build the executable
all: $(TARGET)

# Rule to link the executable
# Original rule using intermediate object files (generally more efficient for incremental builds):
# $(TARGET): $(OBJS)
#	$(CC) $(LDFLAGS) -o $(TARGET) $(OBJS) $(LIBS)

# Updated rule to compile and link source files directly for the target,
# similar to the one-step command line compilation you used.
$(TARGET): $(SRCS) aphorism_utils.h ephemeris.h # Dependencies: source files and their main headers
	$(CC) $(CFLAGS) $(LDFLAGS) -o $(TARGET) $(SRCS) $(LIBS)

# Generic rule to compile .c files into .o files (still useful for individual .o files or other targets)
# This uses automatic variables:
#   $< is the first prerequisite (the .c file)
#   $@ is the target name (the .o file)
.c.o:
	$(CC) $(CFLAGS) -c $< -o $@

# Specific dependencies for object files if they include other headers
# This ensures that if a header file changes, the dependent .c files are recompiled.
main_app.o: main_app.c aphorism_utils.h ephemeris.h
aphorism_utils.o: aphorism_utils.c aphorism_utils.h
ephemeris.o: ephemeris.c ephemeris.h

# Clean up build files
clean:
	rm -f $(TARGET) $(OBJS)

# Phony targets are targets that are not actual files
.PHONY: all clean
