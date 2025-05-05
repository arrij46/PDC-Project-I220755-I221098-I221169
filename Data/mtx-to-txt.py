#!/usr/bin/env python3
# For converting MTX files to a simple graph format for Serial code
input_file = "Data/soc-twitter-follows.mtx"  
output_file = "Data/data.txt"

import sys
import re
import os.path

def convert_mtx_to_graph(input_file, output_file):
    """
    Convert an MTX file to a simple graph format.
    The MTX file is expected to be in coordinate format.
    """
    print(f"Converting {input_file} to {output_file}...")
    
    # Read the MTX file
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    # Parse header information
    header_found = False
    dimensions_found = False
    num_vertices = 0
    num_edges = 0
    edges = []
    
    for line in lines:
        # Skip comments but capture header
        if line.startswith('%'):
            if "%%MatrixMarket" in line:
                header_found = True
                parts = line.strip().split()
                if len(parts) >= 5:
                    format_type = parts[1]
                    representation = parts[2]
                    data_type = parts[3]
                    pattern = parts[4]
                    
                    if format_type != "matrix" or representation != "coordinate":
                        print(f"Warning: Expected 'matrix coordinate' format, got '{format_type} {representation}'")
            continue
        
        # First non-comment line contains dimensions
        if not dimensions_found:
            dimensions_found = True
            parts = line.strip().split()
            if len(parts) >= 3:
                num_rows = int(parts[0])
                num_cols = int(parts[1])
                num_nonzeros = int(parts[2])
                
                # For a graph, we expect rows = cols
                if num_rows != num_cols:
                    print(f"Warning: Non-square matrix ({num_rows}x{num_cols}). Using max dimension as vertex count.")
                
                num_vertices = max(num_rows, num_cols)
                num_edges = num_nonzeros
                
                print(f"Found matrix with {num_vertices} vertices and {num_edges} edges.")
            else:
                print("Error: Invalid dimension line in MTX file.")
                return False
            continue
        
        # Process edge entries
        parts = line.strip().split()
        if len(parts) >= 2:
            row = int(parts[0])
            col = int(parts[1])
            
            # Extract weight if present, default to 1
            weight = 1
            if len(parts) >= 3:
                try:
                    weight = int(parts[2])
                except ValueError:
                    try:
                        weight = float(parts[2])
                    except ValueError:
                        print(f"Warning: Invalid weight '{parts[2]}', using 1 instead.")
            
            # Adjust indices if they're 1-based (common in MTX)
            # Our output format will use 0-based indices for the C++ code
            row -= 1
            col -= 1
            
            edges.append((row, col, weight))
    
    # Write the simple graph format
    with open(output_file, 'w') as f:
        f.write(f"{num_vertices} {num_edges}\n")
        for edge in edges:
            f.write(f"{edge[0]} {edge[1]} {edge[2]}\n")
    
    print(f"Conversion complete: {num_vertices} vertices, {num_edges} edges.")
    return True

def main():
    
    
    
    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' does not exist.")
        return 1
    
    if convert_mtx_to_graph(input_file, output_file):
        print(f"Successfully converted '{input_file}' to '{output_file}'.")
        return 0
    else:
        print(f"Failed to convert '{input_file}'.")
        return 1

if __name__ == "__main__":
    sys.exit(main())