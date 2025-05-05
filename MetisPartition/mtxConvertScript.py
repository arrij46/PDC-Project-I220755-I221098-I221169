# MTX to .graph 
# python3 ConvertScript.py
# gpmetis data.graph

from collections import defaultdict

input_file = "Data/soc-twitter-follows.mtx"
output_file = "MetisPartition/data.graph"

adj_list = defaultdict(set)
n = 0  # Number of nodes
header_skipped = False

try:
    with open(input_file, 'r') as f:
        for line in f:
            if line.startswith('%'):
                continue
            parts = line.strip().split()

            # First non-comment line is the header
            if not header_skipped:
                nrows, ncols, nedges = map(int, parts)
                n = max(nrows, ncols)
                header_skipped = True
                continue

            if len(parts) >= 2:
                u, v = map(int, parts[:2])
                if u == v:
                    continue  # skip self-loops
                adj_list[u].add(v)
                adj_list[v].add(u)

    edges_count = sum(len(neighbors) for neighbors in adj_list.values()) // 2

    with open(output_file, 'w') as f:
        f.write(f"{n} {edges_count}\n")
        for i in range(1, n + 1):
            neighbors = sorted(adj_list.get(i, []))
            f.write(" ".join(map(str, neighbors)) + "\n")

    print(f"Graph successfully converted and saved to {output_file}")

except FileNotFoundError:
    print(f"Error: The file {input_file} was not found.")
except Exception as e:
    print(f"An error occurred: {e}")
