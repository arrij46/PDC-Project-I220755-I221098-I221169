# To run this script
# python3 ConvertScript.py
# gpmetis test.graph


from collections import defaultdict

# Input and output file paths
input_file = "../Data/wiki-RfA_Neutral.mtx"
output_file = "Data.graph"

# Adjacency list to store the graph
adj_list = defaultdict(set)
n = 0  # Number of nodes

try:
    # Read the input file
    with open(input_file, 'r') as f:
        for line in f:
            if line.startswith('%'):  # Skip comment lines
                continue
            parts = line.strip().split()
            if len(parts) == 3:  # Ensure the line has three parts
                u, v, _ = map(int, parts)
                if u == v:
                    continue  # Skip self-loops
                adj_list[u].add(v)
                adj_list[v].add(u)  # Undirected graph
                n = max(n, u, v)  # Update the maximum node number

    # Count edges (each undirected edge counted once)
    edges_count = sum(len(neighbors) for neighbors in adj_list.values()) // 2

    # Write the output file
    with open(output_file, 'w') as f:
        f.write(f"{n} {edges_count}\n")  # Write number of nodes and edges
        for i in range(1, n + 1):
            neighbors = sorted(adj_list[i])  # Ensure neighbors are sorted
            f.write(" ".join(map(str, neighbors)) + "\n")  # Write adjacency list

    print(f"Graph successfully converted and saved to {output_file}")

except FileNotFoundError:
    print(f"Error: The file {input_file} was not found.")
except Exception as e:
    print(f"An error occurred: {e}")