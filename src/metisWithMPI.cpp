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
#include <mpi.h>

using namespace std;

// Structure to represent a weighted edge
struct Edge {
    int source;
    int target;
    int weight;

    Edge(int s, int t, int w) : source(s), target(t), weight(w) {}
};

struct VertexDistance {
    int vertex;
    int distance;

    VertexDistance(int v, int d) : vertex(v), distance(d) {}

    bool operator>(const VertexDistance& other) const {
        return distance > other.distance;
    }
};

class DynamicSSSP {
private:
    int num_vertices;
    vector<vector<Edge>> outgoing_edges;
    vector<vector<Edge>> incoming_edges;
    vector<int> distances;
    vector<int> parents;
    int source_vertex;
    int num_partitions;
    vector<int> vertex_to_partition;
    bool adjust_indices;
    
    // MPI-related members
    int mpi_rank;
    int mpi_size;
    vector<int> local_partitions; // Partitions this process is responsible for
    unordered_map<int, vector<int>> boundary_vertices; // Key: partition, Value: boundary vertices

public:
    DynamicSSSP(int n, int source, bool adjust = false) : 
        num_vertices(n), source_vertex(source), num_partitions(0), adjust_indices(adjust) {
        outgoing_edges.resize(n);
        incoming_edges.resize(n);
        distances.resize(n, numeric_limits<int>::max());
        parents.resize(n, -1);
        vertex_to_partition.resize(n, -1);
        
        MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    }

    // Add edge (unchanged from original)
    void add_edge(int source, int target, int weight) {
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

    // Load partitions (modified for MPI)
    bool load_partitions(const string &filename) {
        ifstream infile(filename);
        if (!infile) {
            cerr << "Error opening partition file: " << filename << "\n";
            return false;
        }

        if (mpi_rank == 0) {
            cout << "Loading partition information from " << filename << endl;
        }
        
        int partition_id;
        int vertex_id = 0;
        
        while (infile >> partition_id && vertex_id < num_vertices) {
            vertex_to_partition[vertex_id] = partition_id;
            num_partitions = max(num_partitions, partition_id + 1);
            vertex_id++;
        }

        if (vertex_id < num_vertices) {
            if (mpi_rank == 0) {
                cerr << "Warning: Partition file has fewer vertices than expected. Read " 
                     << vertex_id << " out of " << num_vertices << " vertices.\n";
            }
            for (int i = vertex_id; i < num_vertices; i++) {
                vertex_to_partition[i] = 0;
            }
        }

        // Broadcast the number of partitions to all processes
        MPI_Bcast(&num_partitions, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (mpi_rank == 0) {
            cout << "Loaded " << num_partitions << " partitions for " << vertex_id << " vertices\n";
        }

        // Distribute partitions among processes
        distribute_partitions();

        // Identify boundary vertices (vertices with edges crossing partitions)
        identify_boundary_vertices();

        return true;
    }

    // Distribute partitions among MPI processes
    void distribute_partitions() {
        // Simple round-robin distribution
        for (int p = 0; p < num_partitions; p++) {
            if (p % mpi_size == mpi_rank) {
                local_partitions.push_back(p);
            }
        }

        // Print local partition assignment
        cout << "Process " << mpi_rank << " responsible for partitions: ";
        for (int p : local_partitions) {
            cout << p << " ";
        }
        cout << endl;
    }

    // Identify boundary vertices (vertices with edges crossing partitions)
    void identify_boundary_vertices() {
        for (int u = 0; u < num_vertices; u++) {
            int u_part = vertex_to_partition[u];
            
            // Check if this vertex is in one of our local partitions
            if (find(local_partitions.begin(), local_partitions.end(), u_part) != local_partitions.end()) {
                for (const Edge& edge : outgoing_edges[u]) {
                    int v = edge.target;
                    int v_part = vertex_to_partition[v];
                    
                    // If target is in a different partition, mark as boundary
                    if (v_part != u_part) {
                        boundary_vertices[u_part].push_back(u);
                        break; // No need to check other edges
                    }
                }
            }
        }
    }

    // Load graph (unchanged from original)
    bool load_graph(const string &filename) {
        ifstream infile(filename);
        if (!infile) {
            cerr << "Error opening graph file: " << filename << "\n";
            return false;
        }

        int n, m;
        infile >> n >> m;

        if (n <= 0 || m < 0) {
            cerr << "Invalid graph dimensions\n";
            return false;
        }

        if (mpi_rank == 0) {
            cout << "Loaded graph with " << n << " vertices and " << m << " edges\n";
        }

        if (n != num_vertices) {
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
        
        while (edge_count < m && infile >> u >> v >> w) {
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
            
            if (mpi_rank == 0 && edge_count % 100000 == 0) {
                cout << "Processed " << edge_count << " edges..." << endl;
            }
        }

        if (edge_count < m) {
            cerr << "Warning: Only read " << edge_count << " edges out of expected " << m << endl;
        }
        
        if (error_count > 0) {
            cerr << "Warning: Encountered " << error_count << " edges with invalid vertices." << endl;
        }

        if (mpi_rank == 0) {
            cout << "Successfully added " << edge_count << " edges to the graph" << endl;
        }
        return edge_count > 0;
    }

    // Process a single partition (modified for MPI)
    void process_partition(int partition_id) {
        // Only process if this is our local partition
        if (find(local_partitions.begin(), local_partitions.end(), partition_id) == local_partitions.end()) {
            return;
        }

        cout << "Process " << mpi_rank << " processing partition " << partition_id << endl;
        
        vector<int> partition_vertices;
        for (int v = 0; v < num_vertices; v++) {
            if (vertex_to_partition[v] == partition_id) {
                partition_vertices.push_back(v);
            }
        }
        
        cout << "Partition " << partition_id << " contains " << partition_vertices.size() 
             << " vertices (processed by rank " << mpi_rank << ")" << endl;
        
        for (int v : partition_vertices) {
            process_vertex(v);
        }
    }
    
    // Process a single vertex's outgoing edges
    void process_vertex(int u) {
        if (distances[u] == numeric_limits<int>::max()) {
            return;
        }
            
        for (const Edge &edge : outgoing_edges[u]) {
            int v = edge.target;
            int weight = edge.weight;
            
            if (distances[u] + weight < distances[v]) {
                distances[v] = distances[u] + weight;
                parents[v] = u;
            }
        }
    }

    // Exchange boundary vertex distances with other processes
    void exchange_boundary_updates() {
        // Create a vector to store all boundary updates
        vector<int> updates_to_send;
        
        // Prepare updates for each partition we own
        for (int p : local_partitions) {
            for (int u : boundary_vertices[p]) {
                updates_to_send.push_back(u); // Vertex ID
                updates_to_send.push_back(distances[u]); // Distance
            }
        }
        
        // First, gather the sizes of updates from all processes
        int local_update_size = updates_to_send.size();
        vector<int> all_update_sizes(mpi_size);
        MPI_Allgather(&local_update_size, 1, MPI_INT, 
                     all_update_sizes.data(), 1, MPI_INT, MPI_COMM_WORLD);
        
        // Calculate displacements for Allgatherv
        vector<int> displacements(mpi_size, 0);
        int total_updates = 0;
        for (int i = 0; i < mpi_size; ++i) {
            displacements[i] = total_updates;
            total_updates += all_update_sizes[i];
        }
        
        // Gather all updates from all processes
        vector<int> all_updates(total_updates);
        MPI_Allgatherv(updates_to_send.data(), local_update_size, MPI_INT,
                      all_updates.data(), all_update_sizes.data(), 
                      displacements.data(), MPI_INT, MPI_COMM_WORLD);
        
        // Process received updates
        for (size_t i = 0; i < all_updates.size(); i += 2) {
            int u = all_updates[i];
            int new_dist = all_updates[i+1];
            
            if (new_dist < distances[u]) {
                distances[u] = new_dist;
                // Note: We don't update parent here as it might not be consistent across partitions
            }
        }
    }

    // Compute SSSP using parallel partitioned approach
    void compute_partitioned_sssp() {
        // Initialize distances
        fill(distances.begin(), distances.end(), numeric_limits<int>::max());
        fill(parents.begin(), parents.end(), -1);

        // Set source distance to 0
        distances[source_vertex] = 0;

        if (mpi_rank == 0) {
            cout << "Starting partitioned SSSP with source vertex " << source_vertex << endl;
        }

        // Process partitions in parallel
        bool changed = true;
        int pass = 1;
        const int MAX_PASSES = 3;
        
        while (changed && pass < MAX_PASSES) {
            changed = false;
            
            if (mpi_rank == 0) {
                cout << "\nRunning additional pass " << pass << endl;
            }
            
            // Process all partitions (each process handles its own partitions)
            for (int p = 0; p < num_partitions; p++) {
                process_partition(p);
            }
            
            // Exchange boundary updates
            exchange_boundary_updates();
            
            // Check if any distances changed (local check only)
            // In a real implementation, we'd need a global reduction here
            for (int v = 0; v < num_vertices; v++) {
                if (vertex_to_partition[v] == -1) continue;
                
                for (const Edge& edge : outgoing_edges[v]) {
                    int u = edge.target;
                    if (distances[v] + edge.weight < distances[u]) {
                        changed = true;
                        break;
                    }
                }
                if (changed) break;
            }
            
            // Global check if any process had changes
            int global_changed = changed ? 1 : 0;
            MPI_Allreduce(MPI_IN_PLACE, &global_changed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
            changed = (global_changed == 1);
            
            pass++;
        }
    }

    // Handle edge weight decrease (modified for MPI)
    void handle_edge_decrease(int u, int v, int new_weight) {
        if (adjust_indices) {
            u--;
            v--;
        }
        
        bool edge_found = false;
        for (Edge &edge : outgoing_edges[u]) {
            if (edge.target == v) {
                edge.weight = new_weight;
                edge_found = true;

                for (Edge &in_edge : incoming_edges[v]) {
                    if (in_edge.source == u) {
                        in_edge.weight = new_weight;
                        break;
                    }
                }

                if (distances[u] + new_weight < distances[v]) {
                    distances[v] = distances[u] + new_weight;
                    parents[v] = u;

                    // Need to process affected partitions
                    int v_partition = vertex_to_partition[v];
                    for (int p = v_partition; p < num_partitions; p++) {
                        process_partition(p);
                    }
                    
                    // Exchange updates after processing
                    exchange_boundary_updates();
                }
                break;
            }
        }

        if (!edge_found && mpi_rank == 0) {
            cerr << "Edge (" << u << "," << v << ") not found in graph\n";
        }
    }

    // Handle edge weight increase (modified for MPI)
    void handle_edge_increase(int u, int v, int new_weight) {
        if (adjust_indices) {
            u--;
            v--;
        }
        
        bool edge_found = false;
        for (Edge &edge : outgoing_edges[u]) {
            if (edge.target == v) {
                edge.weight = new_weight;
                edge_found = true;

                for (Edge &in_edge : incoming_edges[v]) {
                    if (in_edge.source == u) {
                        in_edge.weight = new_weight;
                        break;
                    }
                }

                if (parents[v] == u) {
                    // Recompute in parallel
                    compute_partitioned_sssp();
                }
                break;
            }
        }

        if (!edge_found && mpi_rank == 0) {
            cerr << "Edge (" << u << "," << v << ") not found in graph\n";
        }
    }

    // Print distances (only from rank 0)
    void print_distances(int limit = 20) {
        if (mpi_rank != 0) return;

        cout << "\nCurrent shortest paths from source " << source_vertex << ":\n";
        cout << "----------------------------------------\n";
        cout << "Vertex | Distance | Parent | Partition\n";
        cout << "----------------------------------------\n";

        int count = 0;
        for (int i = 0; i < num_vertices && count < limit; i++) {
            if (distances[i] != numeric_limits<int>::max() || i == source_vertex || i < 10) {
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

    // Print partition stats (only from rank 0)
    void print_partition_stats() {
        if (mpi_rank != 0) return;

        cout << "\nPartition Statistics:\n";
        cout << "----------------------------------------\n";
        cout << "Partition | Vertices | Avg Distance\n";
        cout << "----------------------------------------\n";
        
        for (int p = 0; p < num_partitions; p++) {
            int count = 0;
            long long sum = 0;
            
            for (int v = 0; v < num_vertices; v++) {
                if (vertex_to_partition[v] == p && distances[v] != numeric_limits<int>::max()) {
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
    MPI_Init(&argc, &argv);

    int mpi_rank, mpi_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    if (mpi_rank == 0) {
        cout << "\n\n----------------------------------------\n";
        cout << "----------------------------------------\n";
        cout << "Metis + MPI Parallel Implementation \n";
        cout << "----------------------------------------\n";
        cout << "----------------------------------------\n";
        cout << "Running with " << mpi_size << " MPI processes\n";
    }

    // Default values
    int source_vertex = 0;
    string graph_filename = "MetisPartition/data.graph";
    string partition_filename = "";
    bool adjust_indices = false;
    
    // Parse command line arguments (only rank 0)
    if (mpi_rank == 0) {
        for (int i = 1; i < argc; i++) {
            if (strcmp(argv[i], "--adjust-indices") == 0 || strcmp(argv[i], "-a") == 0) {
                adjust_indices = true;
                cout << "Using index adjustment (1-indexed to 0-indexed conversion)" << endl;
            } else if (i == 1 && argv[i][0] != '-') {
                partition_filename = argv[i];
                cout << "Using partition file: " << partition_filename << endl;
            } else if (i == 2 && argv[i][0] != '-') {
                graph_filename = argv[i];
                cout << "Using graph file: " << graph_filename << endl;
            }
        }
        
        if (partition_filename.empty()) {
            partition_filename = graph_filename + ".part.3";
            cout << "Using default partition file: " << partition_filename << endl;
        }

        cout << "Enter source node (0 to N-1): ";
        cin >> source_vertex;

        if (source_vertex < 0) {
            cerr << "Invalid source node." << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
            return 1;
        }

        cout << "Source node is: " << source_vertex << endl;
    }
    
    // Broadcast source vertex and adjustment flag to all processes
    MPI_Bcast(&source_vertex, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&adjust_indices, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
    
    // Adjust source vertex if needed
    int internal_source = source_vertex;
    if (adjust_indices) {
        internal_source--;
        if (mpi_rank == 0) {
            cout << "Internal source node (after adjustment): " << internal_source << endl;
        }
    }

    // Create the dynamic SSSP object
    DynamicSSSP *sssp = new DynamicSSSP(0, internal_source, adjust_indices);

    // Load the initial graph (all processes load the same graph)
    if (mpi_rank == 0) {
        cout << "Loading initial graph from: " << graph_filename << "\n";
    }

    auto start = chrono::high_resolution_clock::now();
    bool loaded = sssp->load_graph(graph_filename);
    
    if (!loaded) {
        cerr << "Process " << mpi_rank << ": Failed to load initial graph.\n";
        delete sssp;
        MPI_Finalize();
        return 1;
    }

    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end - start;
    if (mpi_rank == 0) {
        cout << fixed << setprecision(6);
        cout << "Graph loading time: " << elapsed.count() << " seconds\n";
    }

    // Load partition information
    start = chrono::high_resolution_clock::now();
    
    bool partitions_loaded = sssp->load_partitions(partition_filename);
    if (!partitions_loaded) {
        cerr << "Process " << mpi_rank << ": Failed to load partition information.\n";
        delete sssp;
        MPI_Finalize();
        return 1;
    }
    
    end = chrono::high_resolution_clock::now();
    elapsed = end - start;
    if (mpi_rank == 0) {
        cout << "Partition loading time: " << elapsed.count() << " seconds\n";
    }

    // Compute initial shortest paths using partitioned approach
    start = chrono::high_resolution_clock::now();
    if (mpi_rank == 0) {
        cout << "\nComputing initial shortest paths from source vertex " << source_vertex 
             << " using partitioned approach...\n";
    }
    
    sssp->compute_partitioned_sssp();
    MPI_Barrier(MPI_COMM_WORLD); // Ensure all processes complete before timing
    
    end = chrono::high_resolution_clock::now();
    elapsed = end - start;
    if (mpi_rank == 0) {
        cout << "Initial partitioned SSSP computation time: " << elapsed.count() << " seconds\n";
    }

    sssp->print_distances();
    sssp->print_partition_stats();

    // Read and process changes from changes.txt (only rank 0 reads file)
    start = chrono::high_resolution_clock::now();
    string changes_filename = "Data/changes.txt";
    int num_changes = 0;
    vector<int> changes_data;

    if (mpi_rank == 0) {
        ifstream changes_file(changes_filename);

        if (!changes_file) {
            cout << "No changes file found at " << changes_filename << ". Skipping edge updates.\n";
        } else {
            cout << "\nProcessing edge updates from " << changes_filename << "...\n";
            changes_file >> num_changes;
            changes_data.resize(num_changes * 3); // Each change has u, v, w

            for (int i = 0; i < num_changes; i++) {
                changes_file >> changes_data[i*3] >> changes_data[i*3+1] >> changes_data[i*3+2];
            }
            changes_file.close();
        }
    }

    // Broadcast number of changes and changes data
    MPI_Bcast(&num_changes, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    if (num_changes > 0) {
        if (mpi_rank != 0) {
            changes_data.resize(num_changes * 3);
        }
        MPI_Bcast(changes_data.data(), num_changes * 3, MPI_INT, 0, MPI_COMM_WORLD);

        for (int i = 0; i < num_changes; i++) {
            int u = changes_data[i*3];
            int v = changes_data[i*3+1];
            int w = changes_data[i*3+2];
            
            if (mpi_rank == 0) {
                cout << "\nProcessing change " << (i + 1) << "/" << num_changes
                     << ": Edge (" << u << "," << v << ") -> weight " << w << "\n";
            }
            
            sssp->update_edge_weight(u, v, w);
            
            if (mpi_rank == 0) {
                sssp->print_distances();
                sssp->print_partition_stats();
            }
        }
    }
    
    end = chrono::high_resolution_clock::now();
    elapsed = end - start;
    if (mpi_rank == 0) {
        cout << "Changes incorporated in time: " << elapsed.count() << " seconds\n";
    }

    delete sssp;
    MPI_Finalize();
    return 0;
}