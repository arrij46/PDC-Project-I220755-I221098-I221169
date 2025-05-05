// To run from ../src
// g++ -o exe src/Serial.cpp
// ./exe
//#define graph_filename "Data/graph.txt"
#define graph_filename "Data/data.txt"
#include <iostream>
#include <vector>
#include <queue>
#include <limits>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <chrono>
#include <functional>  // Added this include for std::function

using namespace std;

// Structure to represent a weighted edge
struct Edge
{
    int source; // Source vertex
    int target; // Target vertex
    int weight; // Edge weight

    Edge(int s, int t, int w) : source(s), target(t), weight(w) {}
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

class DynamicSSSP
{
private:
    int num_vertices;
    vector<vector<Edge>> outgoing_edges; // Adjacency list for outgoing edges
    vector<vector<Edge>> incoming_edges; // Adjacency list for incoming edges
    vector<int> distances;               // Current shortest path distances
    vector<int> parents;                 // Parent pointers for shortest paths
    int source_vertex;                   // Source vertex for SSSP

public:
    // Constructor to initialize the graph
    DynamicSSSP(int n, int source) : num_vertices(n), source_vertex(source)
    {
        outgoing_edges.resize(n);
        incoming_edges.resize(n);
        distances.resize(n, numeric_limits<int>::max());
        parents.resize(n, -1);
    }

    // Add an edge to the graph
    void add_edge(int source, int target, int weight)
    {
        outgoing_edges[source].push_back(Edge(source, target, weight));
        incoming_edges[target].push_back(Edge(source, target, weight));
    }

    // Initialize graph from a file in simple format
    // Format: first line contains number of vertices and edges
    // Following lines contain source, target, weight for each edge
    bool load_graph(const string &filename)
    {
        ifstream infile(filename);
        if (!infile)
        {
            cerr << "Error opening graph file: " << filename << "\n";
            return false;
        }

        int n, m;
        infile >> n >> m;

        if (n <= 0 || m < 0)
        {
            cerr << "Invalid graph dimensions\n";
            return false;
        }

        cout << "Loaded graph with " << n << " vertices and " << m << " edges\n";

        // Reset the graph if needed
        if (n != num_vertices)
        {
            num_vertices = n;
            outgoing_edges.clear();
            incoming_edges.clear();
            outgoing_edges.resize(n);
            incoming_edges.resize(n);
            distances.resize(n, numeric_limits<int>::max());
            parents.resize(n, -1);
        }

        int u, v, w;
        for (int i = 0; i < m; i++)
        {
            if (!(infile >> u >> v >> w))
            {
                cerr << "Error reading edge " << i << endl;
                return false;
            }

            // Debug output
            cout << "Adding edge: " << u << " -> " << v << " with weight " << w << endl;

            // Ensure vertices are within valid range
            if (u < 0 || u >= n || v < 0 || v >= n)
            {
                cerr << "Edge vertices out of range: (" << u << "," << v << ")\n";
                continue;
            }

            add_edge(u, v, w);
        }

        return true;
    }

    // Compute initial SSSP using Dijkstra's algorithm
    void compute_initial_sssp()
    {
        // Initialize distances
        fill(distances.begin(), distances.end(), numeric_limits<int>::max());
        fill(parents.begin(), parents.end(), -1);

        // Set source distance to 0
        distances[source_vertex] = 0;

        cout << "Starting Dijkstra with source vertex " << source_vertex << endl;

        priority_queue<VertexDistance, vector<VertexDistance>, greater<VertexDistance>> pq;
        pq.push(VertexDistance(source_vertex, 0));

        while (!pq.empty())
        {
            int u = pq.top().vertex;
            int dist_u = pq.top().distance;
            pq.pop();

            // Debug output
            cout << "Processing vertex " << u << " with distance " << dist_u << endl;

            // Skip if we've found a better path already
            if (dist_u > distances[u])
                continue;

            // Relax all outgoing edges
            for (const Edge &edge : outgoing_edges[u])
            {
                int v = edge.target;
                int weight = edge.weight;

                cout << "  Checking edge to " << v << " with weight " << weight << endl;

                if (distances[u] + weight < distances[v])
                {
                    distances[v] = distances[u] + weight;
                    parents[v] = u;
                    pq.push(VertexDistance(v, distances[v]));

                    cout << "  Updated distance to " << v << " = " << distances[v] << endl;
                }
            }
        }
    }

    // Handle edge weight decrease (Part 1 of algorithm from paper)
    void handle_edge_decrease(int u, int v, int new_weight)
    {
        // Find and update the edge
        bool edge_found = false;
        for (Edge &edge : outgoing_edges[u])
        {
            if (edge.target == v)
            {
                // Removed unused variable old_weight
                edge.weight = new_weight;
                edge_found = true;

                // Also update the corresponding incoming edge
                for (Edge &in_edge : incoming_edges[v])
                {
                    if (in_edge.source == u)
                    {
                        in_edge.weight = new_weight;
                        break;
                    }
                }

                // Check if this creates a shorter path
                if (distances[u] + new_weight < distances[v])
                {
                    // Update distance and parent
                    distances[v] = distances[u] + new_weight;
                    parents[v] = u;

                    // Propagate the changes using Dijkstra
                    propagate_distance_decrease(v);
                }
                break;
            }
        }

        if (!edge_found)
        {
            cerr << "Edge (" << u << "," << v << ") not found in graph\n";
        }
    }

    // Propagate distance decreases to affected vertices
    void propagate_distance_decrease(int start)
    {
        priority_queue<VertexDistance, vector<VertexDistance>, greater<VertexDistance>> pq;
        pq.push(VertexDistance(start, distances[start]));

        while (!pq.empty())
        {
            int u = pq.top().vertex;
            int dist_u = pq.top().distance;
            pq.pop();

            // Skip if we've found a better path already
            if (dist_u > distances[u])
                continue;

            // Relax all outgoing edges
            for (const Edge &edge : outgoing_edges[u])
            {
                int v = edge.target;
                int weight = edge.weight;

                if (distances[u] + weight < distances[v])
                {
                    distances[v] = distances[u] + weight;
                    parents[v] = u;
                    pq.push(VertexDistance(v, distances[v]));
                }
            }
        }
    }

    // Handle edge weight increase (Part 2 of algorithm from paper)
    void handle_edge_increase(int u, int v, int new_weight)
    {
        // Find and update the edge
        bool edge_found = false;
        for (Edge &edge : outgoing_edges[u])
        {
            if (edge.target == v)
            {
                // Removed unused variable old_weight
                edge.weight = new_weight;
                edge_found = true;

                // Also update the corresponding incoming edge
                for (Edge &in_edge : incoming_edges[v])
                {
                    if (in_edge.source == u)
                    {
                        in_edge.weight = new_weight;
                        break;
                    }
                }

                // Check if this edge was part of the shortest path
                if (parents[v] == u)
                {
                    // Need to recompute paths for v and its descendants
                    vector<int> affected_vertices = identify_affected_vertices(v);
                    recompute_distances(affected_vertices);
                }
                break;
            }
        }

        if (!edge_found)
        {
            cerr << "Edge (" << u << "," << v << ") not found in graph\n";
        }
    }

    // Identify vertices affected by an edge weight increase
    vector<int> identify_affected_vertices(int start)
    {
        vector<int> affected;
        vector<bool> visited(num_vertices, false);
        queue<int> q;

        q.push(start);
        visited[start] = true;
        affected.push_back(start);

        while (!q.empty())
        {
            int u = q.front();
            q.pop();

            // Check all vertices that use u as parent in the shortest path
            for (const Edge &edge : outgoing_edges[u])
            {
                int v = edge.target;
                if (!visited[v] && parents[v] == u)
                {
                    visited[v] = true;
                    affected.push_back(v);
                    q.push(v);
                }
            }
        }

        return affected;
    }

    // Recompute distances for affected vertices
    void recompute_distances(const vector<int> &affected_vertices)
    {
        // Reset distances for affected vertices
        for (int v : affected_vertices)
        {
            if (v != source_vertex)
            { // Don't reset the source
                distances[v] = numeric_limits<int>::max();
                parents[v] = -1;
            }
        }

        // Find all non-affected vertices that point to affected ones
        priority_queue<VertexDistance, vector<VertexDistance>, greater<VertexDistance>> pq;

        for (int v : affected_vertices)
        {
            // Check all incoming edges
            for (const Edge &edge : incoming_edges[v])
            {
                int u = edge.source;
                bool is_affected = find(affected_vertices.begin(), affected_vertices.end(), u) != affected_vertices.end();

                // If the source vertex is not affected and has a valid distance
                if (!is_affected && distances[u] != numeric_limits<int>::max())
                {
                    int new_dist = distances[u] + edge.weight;
                    if (new_dist < distances[v])
                    {
                        distances[v] = new_dist;
                        parents[v] = u;
                    }
                }
            }

            // If we found a new valid distance, add to priority queue
            if (distances[v] != numeric_limits<int>::max())
            {
                pq.push(VertexDistance(v, distances[v]));
            }
        }

        // Run Dijkstra's algorithm starting from these affected vertices
        while (!pq.empty())
        {
            int u = pq.top().vertex;
            int dist_u = pq.top().distance;
            pq.pop();

            // Skip if we've found a better path already
            if (dist_u > distances[u])
                continue;

            // Only process edges to affected vertices
            for (const Edge &edge : outgoing_edges[u])
            {
                int v = edge.target;
                bool is_affected = find(affected_vertices.begin(), affected_vertices.end(), v) != affected_vertices.end();

                if (is_affected && distances[u] + edge.weight < distances[v])
                {
                    distances[v] = distances[u] + edge.weight;
                    parents[v] = u;
                    pq.push(VertexDistance(v, distances[v]));
                }
            }
        }
    }

    // Handle any edge weight update (convenience function)
    void update_edge_weight(int u, int v, int new_weight)
    {
        // Find the current weight
        int current_weight = -1;
        for (const Edge &edge : outgoing_edges[u])
        {
            if (edge.target == v)
            {
                current_weight = edge.weight;
                break;
            }
        }

        if (current_weight == -1)
        {
            cerr << "Edge (" << u << "," << v << ") not found in graph\n";
            return;
        }

        cout << "Updating edge (" << u << "," << v << ") from weight "
             << current_weight << " to " << new_weight << endl;

        if (new_weight < current_weight)
        {
            handle_edge_decrease(u, v, new_weight);
        }
        else if (new_weight > current_weight)
        {
            handle_edge_increase(u, v, new_weight);
        }
        else
        {
            cout << "No weight change, skipping update\n";
        }
    }

    // Print the current shortest path distances
    void print_distances()
    {
        cout << "\nCurrent shortest paths from source " << source_vertex << ":\n";
        cout << "----------------------------------------\n";
        cout << "Vertex | Distance | Parent\n";
        cout << "----------------------------------------\n";

        for (int i = 0; i < num_vertices; i++)
        {
            cout << setw(6) << i << " | ";
            if (distances[i] == numeric_limits<int>::max())
                cout << setw(8) << "INF" << " | ";
            else
                cout << setw(8) << distances[i] << " | ";

            if (parents[i] == -1)
                cout << setw(6) << "-" << "\n";
            else
                cout << setw(6) << parents[i] << "\n";
        }
        cout << "----------------------------------------\n";
    }

    void print_tree_structure()
    {
        cout << "\nShortest Path Tree from source " << source_vertex << ":\n";
        cout << "----------------------------------------\n";

        // Helper lambda to print indentation
        auto print_indent = [](int depth)
        {
            for (int i = 0; i < depth; i++)
                cout << "    ";
        };

        // Helper function to recursively print tree
        std::function<void(int, int)> print_subtree = [&](int vertex, int depth)
        {
            print_indent(depth);
            cout << "└── " << vertex;

            if (distances[vertex] == numeric_limits<int>::max())
                cout << " (INF)\n";
            else
                cout << " (" << distances[vertex] << ")\n";

            // Find all children (vertices that have this vertex as parent)
            for (int i = 0; i < num_vertices; i++)
            {
                if (parents[i] == vertex)
                {
                    print_subtree(i, depth + 1);
                }
            }
        };

        // Start printing from source vertex
        print_subtree(source_vertex, 0);
        cout << "----------------------------------------\n";
    }

    // Print the shortest path to a vertex
    void print_path(int v)
    {
        if (distances[v] == numeric_limits<int>::max())
        {
            cout << "No path exists to vertex " << v << "\n";
            return;
        }

        cout << "Shortest path to vertex " << v << " (distance = " << distances[v] << "): ";
        vector<int> path;

        // Reconstruct the path
        while (v != -1)
        {
            path.push_back(v);
            v = parents[v];
        }

        // Print in reverse order
        for (int i = path.size() - 1; i >= 0; i--)
        {
            cout << path[i];
            if (i > 0)
                cout << " -> ";
        }
        cout << "\n";
    }

    // Add getter for num_vertices
    int get_num_vertices() const
    {
        return num_vertices;
    }

    // Generate a small test graph (for debugging)
    void create_test_graph()
    {
        // Clear any existing graph
        outgoing_edges.clear();
        incoming_edges.clear();

        // Create a small test graph with 5 vertices
        num_vertices = 5;
        outgoing_edges.resize(num_vertices);
        incoming_edges.resize(num_vertices);
        distances.resize(num_vertices, numeric_limits<int>::max());
        parents.resize(num_vertices, -1);

        // Add edges
        add_edge(0, 1, 5); // Source (0) to vertex 1
        add_edge(0, 2, 3); // Source (0) to vertex 2
        add_edge(1, 3, 2); // Vertex 1 to vertex 3
        add_edge(2, 1, 1); // Vertex 2 to vertex 1
        add_edge(2, 3, 6); // Vertex 2 to vertex 3
        add_edge(2, 4, 8); // Vertex 2 to vertex 4
        add_edge(3, 4, 4); // Vertex 3 to vertex 4

        cout << "Created test graph with " << num_vertices << " vertices and 7 edges\n";
    }
};

int main() {

    cout << "\n\n----------------------------------------\n";
    cout << "----------------------------------------\n";
    cout << "Serial Implementation \n";
    cout << "----------------------------------------\n";
    cout << "----------------------------------------\n";
    int source_vertex = 0;

    cout << "Enter source node (0 to N-1): ";
    cin >> source_vertex;

    // You might want to add validation
    if (source_vertex < 0) {
        cerr << "Invalid source node." << endl;
        return 1;
    }

    cout << "Source node is: " << source_vertex << endl;

    // Create the dynamic SSSP object
    DynamicSSSP *sssp = new DynamicSSSP(0, source_vertex);

    // Load the initial graph
    //string graph_filename = "Data/graph.txt";
    cout << "Loading initial graph from: " << graph_filename << "\n";

    auto start = chrono::high_resolution_clock::now();

    bool loaded = sssp->load_graph(graph_filename);
    if (!loaded)
    {
        cerr << "Failed to load initial graph. Exiting.\n";
        delete sssp;
        return 1;
    }

    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end - start;
    cout << fixed << setprecision(6);
    cout << "Loading initial graph from: " << elapsed.count() << " seconds\n";

    start = chrono::high_resolution_clock::now();

    // Compute initial shortest paths
    cout << "\nComputing initial shortest paths from source vertex " << source_vertex << "...\n";
    sssp->compute_initial_sssp();

    end = chrono::high_resolution_clock::now();
    elapsed = end - start;
    cout << fixed << setprecision(6);
    cout << "Initial SSSP computation time: " << elapsed.count() << " seconds\n";

    sssp->print_distances();
    sssp->print_tree_structure();

    // Read and process changes from changes.txt
    start = chrono::high_resolution_clock::now();
    string changes_filename = "Data/changes.txt";
    ifstream changes_file(changes_filename);

    if (!changes_file)
    {
        cerr << "Error: Cannot open " << changes_filename << "\n";
        delete sssp;
        return 1;
    }

    cout << "\nProcessing edge updates from " << changes_filename << "...\n";

    int num_changes;
    changes_file >> num_changes; // Read number of changes

    int u, v, w;
    for (int i = 0; i < num_changes; i++)
    {
        changes_file >> u >> v >> w;
        cout << "\nProcessing change " << (i + 1) << "/" << num_changes
             << ": Edge (" << u << "," << v << ") -> weight " << w << "\n";
        sssp->update_edge_weight(u, v, w);
        sssp->print_distances();
        sssp->print_tree_structure();
    }

    changes_file.close();
    end = chrono::high_resolution_clock::now();
    elapsed = end - start;
    cout << fixed << setprecision(6);
    cout << "Changes incoporated in time: " << elapsed.count() << " seconds\n";

    delete sssp;
    return 0;
}