from collections import defaultdict

def convert_to_metis(input_file, output_file):
    with open(input_file, 'r') as f:
        first_line = f.readline()
        num_vertices, num_edges = map(int, first_line.strip().split())

        adjacency = defaultdict(list)

        for line in f:
            if line.strip():
                u, v, w = map(int, line.strip().split())
                # Convert to 1-indexed
                u += 1
                v += 1
                adjacency[u].append((v, w))
                adjacency[v].append((u, w))  # since undirected

    with open(output_file, 'w') as f:
        f.write(f"{num_vertices} {num_edges} 1\n")  # '1' for edge weights

        for i in range(1, num_vertices + 1):
            neighbors = adjacency.get(i, [])
            line = ' '.join(f"{v} {w}" for v, w in sorted(neighbors))
            f.write(line + "\n")

    print(f"Graph written in METIS format to {output_file}")

if __name__ == "__main__":
    convert_to_metis("Data/graph.txt", "graph1.graph")
    