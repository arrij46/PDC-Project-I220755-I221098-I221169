// Testing METIS partitioning and SSSP with edge weight updates

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <queue>
#include <limits>
#include <algorithm>
using namespace std;

// Structure to represent a weighted edge
struct Edge
{
    int target;
    int weight;

    Edge(int t, int w) : target(t), weight(w) {}
};

// Structure for priority queue
struct VertexDistance
{
    int vertex;
    int distance;

    VertexDistance(int v, int d) : vertex(v), distance(d) {}

    bool operator>(const VertexDistance &other) const
    {
        return distance > other.distance;
    }
};

// Reads the METIS format graph file with edge weights
vector<vector<Edge>> read_weighted_graph(const string &filename)
{
    ifstream infile(filename);
    if (!infile)
    {
        cerr << "Error opening graph file: " << filename << "\n";
        exit(1);
    }

    int num_vertices, num_edges;
    string line;

    // Read header
    getline(infile, line);
    stringstream header(line);
    header >> num_vertices >> num_edges;

    bool weighted = false;
    int fmt;
    if (header >> fmt && fmt == 1)
    {
        weighted = true;
    }

    vector<vector<Edge>> adjacency_list(num_vertices);

    for (int i = 0; i < num_vertices; ++i)
    {
        if (!getline(infile, line))
            break;
        stringstream ss(line);
        int neighbor, weight = 1; // Default weight is 1

        while (ss >> neighbor)
        {
            if (weighted && (ss >> weight))
            {
                // If weighted, read the weight
            }
            adjacency_list[i].push_back(Edge(neighbor - 1, weight)); // 1-based to 0-based
        }
    }

    return adjacency_list;
}

// Reads the METIS partition output file
vector<int> read_partition_file(const string &filename)
{
    ifstream infile(filename);
    if (!infile)
    {
        cerr << "Error opening partition file: " << filename << "\n";
        exit(1);
    }

    vector<int> partitions;
    string line;
    while (getline(infile, line))
    {
        partitions.push_back(stoi(line));
    }

    return partitions;
}

// Initialize single source shortest path distances
vector<int> initialize_sssp(const vector<vector<Edge>> &graph, int source)
{
    int n = graph.size();
    vector<int> dist(n, numeric_limits<int>::max());
    dist[source] = 0;

    priority_queue<VertexDistance, vector<VertexDistance>, greater<VertexDistance>> pq;
    pq.push(VertexDistance(source, 0));

    while (!pq.empty())
    {
        int u = pq.top().vertex;
        int d = pq.top().distance;
        pq.pop();

        if (d > dist[u])
            continue;

        for (const Edge &e : graph[u])
        {
            int v = e.target;
            int weight = e.weight;

            if (dist[u] + weight < dist[v])
            {
                dist[v] = dist[u] + weight;
                pq.push(VertexDistance(v, dist[v]));
            }
        }
    }

    return dist;
}

// Process an edge weight update as per the algorithm in the paper
void process_edge_update(vector<vector<Edge>> &graph, vector<int> &dist,
                         int u, int v, int new_weight, const vector<int> &partitions,
                         const unordered_map<int, vector<int>> &partition_map)
{

    // Find the edge and update its weight
    bool found = false;
    int old_weight = 0;

    for (auto &edge : graph[u])
    {
        if (edge.target == v)
        {
            old_weight = edge.weight;
            edge.weight = new_weight;
            found = true;
            break;
        }
    }

    if (!found)
    {
        cerr << "Edge (" << u << ", " << v << ") not found in the graph.\n";
        return;
    }

    int part_u = partitions[u];
    int part_v = partitions[v];

    // If the edge weight decreases
    if (new_weight < old_weight)
    {
        // Update according to decremental update procedure
        if (dist[u] + new_weight < dist[v])
        {
            dist[v] = dist[u] + new_weight;

            // Process affected vertices in each partition sequentially
            for (int p = 0; p < partition_map.size(); ++p)
            {
                priority_queue<VertexDistance, vector<VertexDistance>, greater<VertexDistance>> pq;

                // Start with vertices in the partition that need updating
                if (p == part_v)
                {
                    pq.push(VertexDistance(v, dist[v]));
                }

                // Process the partition
                while (!pq.empty())
                {
                    int x = pq.top().vertex;
                    int d = pq.top().distance;
                    pq.pop();

                    if (d > dist[x])
                        continue;

                    // Relax all neighbors of x
                    for (const Edge &e : graph[x])
                    {
                        int y = e.target;

                        if (partitions[y] == p && dist[x] + e.weight < dist[y])
                        {
                            dist[y] = dist[x] + e.weight;
                            pq.push(VertexDistance(y, dist[y]));
                        }
                        else if (partitions[y] != p && dist[x] + e.weight < dist[y])
                        {
                            // Cross-partition update
                            dist[y] = dist[x] + e.weight;
                        }
                    }
                }
            }
        }
    }
    // If the edge weight increases
    else if (new_weight > old_weight)
    {
        // Check if the edge is on the shortest path
        if (dist[v] == dist[u] + old_weight)
        {
            // Need to potentially update the distance of v and its dependents

            // Step 1: Reset potentially affected vertices
            vector<bool> affected(graph.size(), false);
            queue<int> q;
            q.push(v);
            affected[v] = true;

            while (!q.empty())
            {
                int x = q.front();
                q.pop();

                // Reset the distance to infinity (except for the source)
                if (dist[x] != 0)
                {
                    dist[x] = numeric_limits<int>::max();
                }

                // Find affected descendants
                for (const Edge &e : graph[x])
                {
                    int y = e.target;
                    if (!affected[y] && dist[y] == dist[x] + e.weight)
                    {
                        affected[y] = true;
                        q.push(y);
                    }
                }
            }

            // Step 2: Recompute distances for affected vertices in each partition
            for (int p = 0; p < partition_map.size(); ++p)
            {
                priority_queue<VertexDistance, vector<VertexDistance>, greater<VertexDistance>> pq;

                // Initialize with boundary vertices that have a valid distance
                for (int vertex : partition_map.at(p))
                {
                    if (affected[vertex])
                    {
                        for (const Edge &e : graph[vertex])
                        {
                            int y = e.target;
                            if (partitions[y] != p && dist[y] != numeric_limits<int>::max())
                            {
                                if (dist[y] + e.weight < dist[vertex] || dist[vertex] == numeric_limits<int>::max())
                                {
                                    dist[vertex] = dist[y] + e.weight;
                                }
                            }
                        }

                        if (dist[vertex] != numeric_limits<int>::max())
                        {
                            pq.push(VertexDistance(vertex, dist[vertex]));
                        }
                    }
                }

                // Process the partition
                while (!pq.empty())
                {
                    int x = pq.top().vertex;
                    int d = pq.top().distance;
                    pq.pop();

                    if (d > dist[x])
                        continue;

                    // Relax all neighbors of x in the same partition
                    for (const Edge &e : graph[x])
                    {
                        int y = e.target;

                        if (partitions[y] == p && affected[y] &&
                            (dist[x] + e.weight < dist[y] || dist[y] == numeric_limits<int>::max()))
                        {
                            dist[y] = dist[x] + e.weight;
                            pq.push(VertexDistance(y, dist[y]));
                        }
                    }
                }
            }
        }
    }
}

// Process each partition for the initial SSSP computation
void process_partition_sssp(int part_id, const vector<int> &vertices,
                            const vector<vector<Edge>> &graph, vector<int> &dist)
{
    cout << "\nProcessing Partition " << part_id << " for SSSP:\n";

    priority_queue<VertexDistance, vector<VertexDistance>, greater<VertexDistance>> pq;

    // Add all vertices in this partition that have a valid distance
    for (int v : vertices)
    {
        if (dist[v] != numeric_limits<int>::max())
        {
            pq.push(VertexDistance(v, dist[v]));
        }
    }

    // Process the partition using Dijkstra's algorithm
    while (!pq.empty())
    {
        int u = pq.top().vertex;
        int d = pq.top().distance;
        pq.pop();

        if (d > dist[u])
            continue;

        for (const Edge &e : graph[u])
        {
            int v = e.target;
            int weight = e.weight;

            if (dist[u] + weight < dist[v])
            {
                dist[v] = dist[u] + weight;
                // Only process vertices in this partition
                if (find(vertices.begin(), vertices.end(), v) != vertices.end())
                {
                    pq.push(VertexDistance(v, dist[v]));
                }
            }
        }
    }

    // Print the distances for vertices in this partition
    for (int v : vertices)
    {
        if (dist[v] != numeric_limits<int>::max())
        {
            cout << "Node " << v + 1 << " → Distance: ";
            /*if (dist[v] == numeric_limits<int>::max())
                cout << "INF" << endl;
            else
            */
            cout << dist[v] << endl;
        }
    }
}

int main(int argc, char *argv[])
{
    if (argc < 4)
    {
        cerr << "Usage: ./exe <path_to_graph_file> <num_partitions> <source_vertex>\n";
        return 1;
    }

    string graph_path = argv[1];
    int num_parts = stoi(argv[2]);
    int source = stoi(argv[3]) - 1; // Convert from 1-based to 0-based

    // Construct paths for graph and partition files
    string part_file = graph_path + ".part." + to_string(num_parts);

    // Read graph and partition files
    vector<vector<Edge>> graph = read_weighted_graph(graph_path);
    vector<int> partitions = read_partition_file(part_file);

    // Validate partition count
    if (partitions.size() != graph.size())
    {
        cerr << "Error: Partition file size does not match the number of vertices in the graph.\n";
        return 1;
    }

    // Group vertices by partition
    unordered_map<int, vector<int>> partition_map;
    for (int i = 0; i < partitions.size(); ++i)
    {
        partition_map[partitions[i]].push_back(i);
    }

    // Initialize SSSP from the source vertex
    cout << "Initializing SSSP from source vertex " << source + 1 << "...\n";
    vector<int> dist = initialize_sssp(graph, source);

    // Process each partition sequentially using the approach from the paper
    for (int i = 0; i < num_parts; ++i)
    {
        process_partition_sssp(i, partition_map[i], graph, dist);
    }

    // Simulate some edge updates
    cout << "\n=== Simulating Edge Updates ===\n";

    // Example: Decrease the weight of an edge
    int u = 0, v = 1, new_weight = 1;
    cout << "Decreasing weight of edge (" << u + 1 << ", " << v + 1 << ") to " << new_weight << endl;
    process_edge_update(graph, dist, u, v, new_weight, partitions, partition_map);

    // Print updated distances
    cout << "\nUpdated distances after edge weight decrease:\n";
    for (int i = 0; i < dist.size(); ++i)
    {
        cout << "Node " << i + 1 << " → Distance: ";
        if (dist[i] == numeric_limits<int>::max())
            cout << "INF" << endl;
        else
            cout << dist[i] << endl;
    }

    // Example: Increase the weight of an edge
    u = 0;
    v = 1;
    new_weight = 10;
    cout << "\nIncreasing weight of edge (" << u + 1 << ", " << v + 1 << ") to " << new_weight << endl;
    process_edge_update(graph, dist, u, v, new_weight, partitions, partition_map);

    // Print updated distances
    cout << "\nUpdated distances after edge weight increase:\n";
    for (int i = 0; i < dist.size(); ++i)
    {
        if (dist[i] != numeric_limits<int>::max())
        {
            cout << "Node " << i + 1 << " → Distance: ";
            if (dist[i] == numeric_limits<int>::max())
                cout << "INF" << endl;
            else
                cout << dist[i] << endl;
        }
    }

    return 0;
}