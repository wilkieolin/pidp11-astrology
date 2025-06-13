# Makefile for the Astrology Program

# Compiler and flags
CC = cc
CFLAGS = -O      # Basic optimization, K&R C often doesn't have many CFLAGS
                 # For debugging, you might use: CFLAGS = -g
LDFLAGS =        # Linker flags
LIBS = -lm       # Math library for ephemeris.c

# Source files
MAIN_SRCS = main_app.c aphorism_utils.c ephemeris.c

# Source files for the bare application
BARE_SRCS = bare_app.c aphorism_utils.c ephemeris.c

# Object files (derived from SRCS)
# OBJS = $(MAIN_SRCS:.c=.o) # Kept for reference if you revert to .o linking
OBJS = main_app.o aphorism_utils.o ephemeris.o bare_app.o # List all possible .o files

# Executable name
MAIN_TARGET = my_astrology_app
BARE_TARGET = bare_astrology_app

# Default target: build the executable
all: $(MAIN_TARGET)

# Target for the bare application
bare_app: $(BARE_TARGET)

# Rule to link the executable
# Original rule using intermediate object files (generally more efficient for incremental builds):
# $(MAIN_TARGET): $(OBJS)
#	$(CC) $(LDFLAGS) -o $(TARGET) $(OBJS) $(LIBS)

# Rule to compile and link main_app.c and its dependencies
$(MAIN_TARGET): $(MAIN_SRCS) aphorism_utils.h ephemeris.h
	$(CC) $(CFLAGS) $(LDFLAGS) -o $(MAIN_TARGET) $(MAIN_SRCS) $(LIBS)

# Rule to compile and link bare_app.c and its dependencies
$(BARE_TARGET): $(BARE_SRCS) aphorism_utils.h ephemeris.h
	$(CC) $(CFLAGS) $(LDFLAGS) -o $(BARE_TARGET) $(BARE_SRCS) $(LIBS)

# Generic rule to compile .c files into .o files (still useful for individual .o files or other targets)
# This uses automatic variables:
#   $< is the first prerequisite (the .c file)
#   $@ is the target name (the .o file)
.c.o:
	$(CC) $(CFLAGS) -c $< -o $@

# Specific dependencies for object files if they include other headers
# This ensures that if a header file changes, the dependent .c files are recompiled.
main_app.o: main_app.c aphorism_utils.h ephemeris.h
bare_app.o: bare_app.c aphorism_utils.h ephemeris.h
aphorism_utils.o: aphorism_utils.c aphorism_utils.h
ephemeris.o: ephemeris.c ephemeris.h

# Clean up build files
clean:
	rm -f $(MAIN_TARGET) $(BARE_TARGET) $(OBJS)

# Phony targets are targets that are not actual files
.PHONY: all clean bare_app