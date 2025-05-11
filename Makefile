# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++11 -Wall -O2
OMPFLAGS = -fopenmp

# Source files
SRC_DIR = src
DATA_DIR = Data
METIS_DIR = MetisPartition

SERIAL_SRC = $(SRC_DIR)/Serial.cpp
OMP_SRC = $(SRC_DIR)/serialWithOpenMP.cpp
METIS_SRC = $(SRC_DIR)/MetisSequential.cpp

# Executables
SERIAL_EXE = serial
OMP_EXE = serial_omp
METIS_EXE = sequential

# Default target
all: prepare_data $(SERIAL_EXE) $(OMP_EXE) $(METIS_EXE)

# Compile targets gprof
$(SERIAL_EXE): $(SERIAL_SRC)
	$(CXX) $(CXXFLAGS) -pg -o $@ $<

$(OMP_EXE): $(OMP_SRC)
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) -pg -o $@ $<

$(METIS_EXE): $(METIS_SRC)
	$(CXX) $(CXXFLAGS) -pg -o $@ $<

# Run all programs
run: all
	./$(SERIAL_EXE) $(DATA_DIR)/data.txt
		gprof $(SERIAL_EXE) gmon.out > gprof-serial.txt

	./$(OMP_EXE) $(DATA_DIR)/data.txt
		gprof $(OMP_EXE) gmon.out > gprof-openmp.txt

	./$(METIS_EXE) $(METIS_DIR)/data.graph.part.3 --adjust-indices
		gprof $(METIS_EXE) gmon.out > gprof-sequentialMetis.txt


# Prepare data only if needed
prepare_data:
	@if [ ! -f $(DATA_DIR)/data.txt ]; then \
		echo "Generating data.txt from mtx..."; \
		python3 $(DATA_DIR)/mtx-to-txt.py; \
	fi
	@if [ ! -f $(METIS_DIR)/data.graph ]; then \
		echo "Generating data.graph from mtx..."; \
		python3 $(METIS_DIR)/mtxConvertScript.py; \
	fi
	@if [ ! -f $(METIS_DIR)/data.graph.part.3 ]; then \
		echo "Running gpmetis..."; \
		cd $(METIS_DIR) && gpmetis data.graph 3; \
	fi

# Clean build artifacts
clean:
	rm -f $(SERIAL_EXE) $(OMP_EXE) $(METIS_EXE)
