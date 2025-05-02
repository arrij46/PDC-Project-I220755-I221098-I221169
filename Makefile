# To run the MetisSequential program

# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++11 -Wall -O2

# Target executable
TARGET = sequential

# Source files
SRCS = src/MetisSequential.cpp

# Object files
OBJS = $(SRCS:.cpp=.o)

# Default target
all: $(TARGET)

# Build target
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

# Run the program
run: $(TARGET)
	./$(TARGET) MetisPartition/Data.graph 3 11372

# Clean up build files
clean:
	rm -f $(TARGET) $(OBJS)
