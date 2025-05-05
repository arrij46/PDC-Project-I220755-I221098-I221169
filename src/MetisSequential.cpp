#include <iostream>
#include <vector>
#include <queue>
#include <limits>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <chrono>
#include <string>
#include <unordered_map>
#include <functional>
#include <cstring>

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
    int num_partitions;                  // Number of partitions
    vector<int> vertex_to_partition;     // Maps each vertex to its partition
    bool adjust_indices;                 // Whether to adjust 1-indexed to 0-indexed

public:
    // Constructor to initialize the graph
    DynamicSSSP(int n, int source, bool adjust = false) : 
        num_vertices(n), source_vertex(source), num_partitions(0), adjust_indices(adjust)
    {
        outgoing_edges.resize(n);
        incoming_edges.resize(n);
        distances.resize(n, numeric_limits<int>::max());
        parents.resize(n, -1);
        vertex_to_partition.resize(n, -1);
    }

    // Add an edge to the graph
    void add_edge(int source, int target, int weight)
    {
        // Adjust indices if needed (convert from 1-indexed to 0-indexed)
        if (adjust_indices) {
            source--;
            target--;
        }

        if (source < 0 || source >= num_vertices || target < 0 || target >= num_vertices) {
            cerr << "Edge vertices out of range: (" << source << "," << target << ")\n";
            return;
        }

        outgoing_edges[source].push_back(Edge(source, target, weight));
        incoming_edges[target].push_back(Edge(source, target, weight));
    }

    // Load the METIS partition file
    bool load_partitions(const string &filename)
    {
        ifstream infile(filename);
        if (!infile)
        {
            cerr << "Error opening partition file: " << filename << "\n";
            return false;
        }

        cout << "Loading partition information from " << filename << endl;
        
        // Read partition assignment for each vertex
        int partition_id;
        int vertex_id = 0;
        
        while (infile >> partition_id && vertex_id < num_vertices)
        {
            vertex_to_partition[vertex_id] = partition_id;
            num_partitions = max(num_partitions, partition_id + 1);
            vertex_id++;
        }

        if (vertex_id < num_vertices)
        {
            cerr << "Warning: Partition file has fewer vertices than expected. Read " 
                 << vertex_id << " out of " << num_vertices << " vertices.\n";
            // Fill remaining vertices with default partition
            for (int i = vertex_id; i < num_vertices; i++)
            {
                vertex_to_partition[i] = 0;
            }
        }

        cout << "Loaded " << num_partitions << " partitions for " << vertex_id << " vertices\n";
        return true;
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
            vertex_to_partition.resize(n, -1);
        }

        int u, v, w;
        int edge_count = 0;
        int error_count = 0;
        
        while (edge_count < m && infile >> u >> v >> w)
        {
            // Validate edge before adding
            if (adjust_indices) {
                if (u <= 0 || u > n || v <= 0 || v > n) {
                    cerr << "Edge vertices out of range (before adjustment): (" << u << "," << v << ")\n";
                    error_count++;
                    continue;
                }
            } else {
                if (u < 0 || u >= n || v < 0 || v >= n) {
                    cerr << "Edge vertices out of range: (" << u << "," << v << ")\n";
                    error_count++;
                    continue;
                }
            }

            add_edge(u, v, w);
            edge_count++;
            
            // Print progress for large graphs
            if (edge_count % 100000 == 0) {
                cout << "Processed " << edge_count << " edges..." << endl;
            }
        }

        if (edge_count < m) {
            cerr << "Warning: Only read " << edge_count << " edges out of expected " << m << endl;
        }
        
        if (error_count > 0) {
            cerr << "Warning: Encountered " << error_count << " edges with invalid vertices." << endl;
        }

        cout << "Successfully added " << edge_count << " edges to the graph" << endl;
        return edge_count > 0; // Return true if at least some edges were loaded
    }

    // Process a single partition
    void process_partition(int partition_id)
    {
        cout << "\nProcessing partition " << partition_id << endl;
        
        // Find vertices in this partition
        vector<int> partition_vertices;
        for (int v = 0; v < num_vertices; v++)
        {
            if (vertex_to_partition[v] == partition_id)
            {
                partition_vertices.push_back(v);
            }
        }
        
        cout << "Partition " << partition_id << " contains " << partition_vertices.size() << " vertices" << endl;
        
        // Process each vertex in the partition
        for (int v : partition_vertices)
        {
            process_vertex(v);
        }
    }
    
    // Process a single vertex's outgoing edges
    void process_vertex(int u)
    {
        // Skip if vertex distance is not yet determined
        if (distances[u] == numeric_limits<int>::max())
            return;
            
        // Process all outgoing edges
        for (const Edge &edge : outgoing_edges[u])
        {
            int v = edge.target;
            int weight = edge.weight;
            
            if (distances[u] + weight < distances[v])
            {
                distances[v] = distances[u] + weight;
                parents[v] = u;
                
                // Note: We don't recursively process v here, as it will be processed
                // either in this partition or in a later partition depending on its assignment
            }
        }
    }

    // Compute initial SSSP using partitioned Dijkstra's algorithm
    void compute_partitioned_sssp()
    {
        // Initialize distances
        fill(distances.begin(), distances.end(), numeric_limits<int>::max());
        fill(parents.begin(), parents.end(), -1);

        // Set source distance to 0
        distances[source_vertex] = 0;

        cout << "Starting partitioned SSSP with source vertex " << source_vertex << endl;

        // Process each partition sequentially
        for (int p = 0; p < num_partitions; p++)
        {
            process_partition(p);
        }
        
        // We may need to run multiple passes because changes in one partition 
        // might affect vertices in earlier partitions
        bool changed = true;
        int pass = 1;
        const int MAX_PASSES = 3; // Limit number of passes to avoid infinite loops
        
        while (changed && pass < MAX_PASSES)
        {
            changed = false;
            cout << "\nRunning additional pass " << pass << endl;
            
            for (int p = 0; p < num_partitions; p++)
            {
                // Check if any distances improve in this partition
                vector<int> old_distances = distances;
                
                process_partition(p);
                
                // Check if any distances changed
                for (int v = 0; v < num_vertices; v++)
                {
                    if (old_distances[v] != distances[v])
                    {
                        changed = true;
                        break;
                    }
                }
            }
            
            pass++;
        }
    }

    // Handle edge weight decrease
    void handle_edge_decrease(int u, int v, int new_weight)
    {
        // Adjust indices if needed
        if (adjust_indices) {
            u--;
            v--;
        }
        
        // Find and update the edge
        bool edge_found = false;
        for (Edge &edge : outgoing_edges[u])
        {
            if (edge.target == v)
            {
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

                    // Need to process affected partitions
                    int v_partition = vertex_to_partition[v];
                    for (int p = v_partition; p < num_partitions; p++)
                    {
                        process_partition(p);
                    }
                }
                break;
            }
        }

        if (!edge_found)
        {
            cerr << "Edge (" << u << "," << v << ") not found in graph\n";
        }
    }

    // Handle edge weight increase
    void handle_edge_increase(int u, int v, int new_weight)
    {
        // Adjust indices if needed
        if (adjust_indices) {
            u--;
            v--;
        }
        
        // Find and update the edge
        bool edge_found = false;
        for (Edge &edge : outgoing_edges[u])
        {
            if (edge.target == v)
            {
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
                    // Recompute distances from scratch for simplicity
                    // A more optimized version could identify affected partitions
                    compute_partitioned_sssp();
                }
                break;
            }
        }

        if (!edge_found)
        {
            cerr << "Edge (" << u << "," << v << ") not found in graph\n";
        }
    }

    // Handle any edge weight update (convenience function)
    void update_edge_weight(int u, int v, int new_weight)
    {
        // Note: adjustment will happen in the specific handlers
        int u_idx = adjust_indices ? u - 1 : u;
        int v_idx = adjust_indices ? v - 1 : v;
        
        // Find the current weight
        int current_weight = -1;
        for (const Edge &edge : outgoing_edges[u_idx])
        {
            if (edge.target == v_idx)
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
    void print_distances(int limit = 20)
    {
        cout << "\nCurrent shortest paths from source " << source_vertex << ":\n";
        cout << "----------------------------------------\n";
        cout << "Vertex | Distance | Parent | Partition\n";
        cout << "----------------------------------------\n";

        // Print only up to limit vertices to avoid overwhelming output
        int count = 0;
        for (int i = 0; i < num_vertices && count < limit; i++)
        {
            // Only print vertices with finite distances or special cases
            if (distances[i] != numeric_limits<int>::max() || i == source_vertex || i < 10)
            {
                cout << setw(6) << i << " | ";
                if (distances[i] == numeric_limits<int>::max())
                    cout << setw(8) << "INF" << " | ";
                else
                    cout << setw(8) << distances[i] << " | ";

                if (parents[i] == -1)
                    cout << setw(6) << "-" << " | ";
                else
                    cout << setw(6) << parents[i] << " | ";

                cout << setw(8) << vertex_to_partition[i] << "\n";
                count++;
            }
        }
        
        if (count < num_vertices) {
            cout << "... (showing " << count << " out of " << num_vertices << " vertices)\n";
        }
        
        cout << "----------------------------------------\n";
    }

    // Print summary statistics by partition
    void print_partition_stats()
    {
        cout << "\nPartition Statistics:\n";
        cout << "----------------------------------------\n";
        cout << "Partition | Vertices | Avg Distance\n";
        cout << "----------------------------------------\n";
        
        for (int p = 0; p < num_partitions; p++)
        {
            int count = 0;
            long long sum = 0;
            
            for (int v = 0; v < num_vertices; v++)
            {
                if (vertex_to_partition[v] == p && distances[v] != numeric_limits<int>::max())
                {
                    sum += distances[v];
                    count++;
                }
            }
            
            cout << setw(9) << p << " | ";
            cout << setw(8) << count << " | ";
            
            if (count > 0)
                cout << setw(12) << fixed << setprecision(2) << (double)sum / count << "\n";
            else
                cout << setw(12) << "N/A" << "\n";
        }
        cout << "----------------------------------------\n";
    }

    // Print the shortest path tree structure
    void print_tree_structure(int max_depth = 3, int max_children = 5)
    {
        cout << "\nShortest Path Tree from source " << source_vertex << ":\n";
        cout << "----------------------------------------\n";

        // Helper lambda to print indentation
        auto print_indent = [](int depth)
        {
            for (int i = 0; i < depth; i++)
                cout << "    ";
        };

        // Use std::function for recursive lambda
        std::function<void(int, int)> print_subtree = [&](int vertex, int depth)
        {
            if (depth > max_depth) {
                print_indent(depth);
                cout << "... (max depth reached)\n";
                return;
            }
            
            print_indent(depth);
            cout << "└── " << vertex;

            if (distances[vertex] == numeric_limits<int>::max())
                cout << " (INF)";
            else
                cout << " (" << distances[vertex] << ")";

            cout << " [P" << vertex_to_partition[vertex] << "]\n";

            // Find all children (vertices that have this vertex as parent)
            vector<int> children;
            for (int i = 0; i < num_vertices; i++)
            {
                if (parents[i] == vertex)
                {
                    children.push_back(i);
                }
            }
            
            // Limit the number of children displayed
            int shown = 0;
            for (int child : children) {
                if (shown < max_children) {
                    print_subtree(child, depth + 1);
                    shown++;
                } else {
                    print_indent(depth + 1);
                    cout << "... (" << (children.size() - max_children) << " more children)\n";
                    break;
                }
            }
        };

        // Start printing from source vertex
        print_subtree(source_vertex, 0);
        cout << "----------------------------------------\n";
    }

    // Get the number of vertices
    int get_num_vertices() const
    {
        return num_vertices;
    }
    
    // Get the number of partitions
    int get_num_partitions() const
    {
        return num_partitions;
    }
};

int main(int argc, char* argv[]) {
    cout << "\n\n----------------------------------------\n";
    cout << "----------------------------------------\n";
    cout << "Metis + Serial Implementation \n";
    cout << "----------------------------------------\n";
    cout << "----------------------------------------\n";
    
    // Default values
    int source_vertex = 0;
    string graph_filename = "MetisPartition/data.graph";
    string partition_filename = "";
    bool adjust_indices = false; // Default to 0-indexed
    
    // Parse command line arguments
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--adjust-indices") == 0 || strcmp(argv[i], "-a") == 0) {
            adjust_indices = true;
            cout << "Using index adjustment (1-indexed to 0-indexed conversion)" << endl;
        } else if (i == 1 && argv[i][0] != '-') {
            // First non-flag argument is the partition file
            partition_filename = argv[i];
            cout << "Using partition file: " << partition_filename << endl;
        } else if (i == 2 && argv[i][0] != '-') {
            // Second non-flag argument is the graph file
            graph_filename = argv[i];
            cout << "Using graph file: " << graph_filename << endl;
        }
    }
    
    if (partition_filename.empty()) {
        // Auto-generate partition filename if not provided
        partition_filename = graph_filename + ".part.3";
        cout << "Using default partition file: " << partition_filename << endl;
    }

    cout << "Enter source node (0 to N-1): ";
    cin >> source_vertex;

    // You might want to add validation
    if (source_vertex < 0) {
        cerr << "Invalid source node." << endl;
        return 1;
    }

    cout << "Source node is: " << source_vertex << endl;
    
    // Adjust source vertex if needed
    int internal_source = source_vertex;
    if (adjust_indices) {
        internal_source--;
        cout << "Internal source node (after adjustment): " << internal_source << endl;
    }

    // Create the dynamic SSSP object
    DynamicSSSP *sssp = new DynamicSSSP(0, internal_source, adjust_indices);

    // Load the initial graph
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
    cout << "Graph loading time: " << elapsed.count() << " seconds\n";

    // Load partition information
    start = chrono::high_resolution_clock::now();
    
    bool partitions_loaded = sssp->load_partitions(partition_filename);
    if (!partitions_loaded)
    {
        cerr << "Failed to load partition information. Exiting.\n";
        delete sssp;
        return 1;
    }
    
    end = chrono::high_resolution_clock::now();
    elapsed = end - start;
    cout << "Partition loading time: " << elapsed.count() << " seconds\n";

    // Compute initial shortest paths using partitioned approach
    start = chrono::high_resolution_clock::now();
    cout << "\nComputing initial shortest paths from source vertex " << source_vertex 
         << " using partitioned approach...\n";
    sssp->compute_partitioned_sssp();

    end = chrono::high_resolution_clock::now();
    elapsed = end - start;
    cout << "Initial partitioned SSSP computation time: " << elapsed.count() << " seconds\n";

    sssp->print_distances();
    sssp->print_partition_stats();
    sssp->print_tree_structure();

    // Read and process changes from changes.txt
    start = chrono::high_resolution_clock::now();
    string changes_filename = "Data/changes.txt";
    ifstream changes_file(changes_filename);

    if (!changes_file)
    {
        cout << "No changes file found at " << changes_filename << ". Skipping edge updates.\n";
    }
    else
    {
        cout << "\nProcessing edge updates from " << changes_filename << "...\n";

        int num_changes;
        changes_file >> num_changes; // Read number of changes

        int u, v, w;
        for (int i = 0; i < num_changes; i++)
        {
            if (!(changes_file >> u >> v >> w)) {
                cerr << "Error reading change " << (i+1) << ". Stopping edge updates.\n";
                break;
            }
            
            cout << "\nProcessing change " << (i + 1) << "/" << num_changes
                 << ": Edge (" << u << "," << v << ") -> weight " << w << "\n";
                 
            sssp->update_edge_weight(u, v, w);
            sssp->print_distances();
            sssp->print_partition_stats();
        }

        changes_file.close();
    }
    
    end = chrono::high_resolution_clock::now();
    elapsed = end - start;
    cout << "Changes incorporated in time: " << elapsed.count() << " seconds\n";

    delete sssp;
    return 0;
}