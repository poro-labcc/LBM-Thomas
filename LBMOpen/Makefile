# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++17 -Wall -O3 -fopenmp -I.

# Find all source files recursively
SRCS := $(shell find . -name '*.cpp')
OBJS := $(SRCS:.cpp=.o)

# Executable name
TARGET = LBMOpen

# Default target
all: $(TARGET)

# Linking the final executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(TARGET)

# Compiling source files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean up
clean:
	rm -f $(OBJS) $(TARGET) velocity/*.vtk force.txt velocity_profile*.txt density_profile*.txt

.PHONY: all clean

