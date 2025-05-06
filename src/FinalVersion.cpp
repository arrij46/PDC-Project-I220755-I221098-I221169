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

#include <mpi.h> // Include MPI header

#include <omp.h> // Include OpenMP header



using namespace std;



// Structure to represent a weighted edge

struct Edge

{

    int source; // Source vertex

    int target; // Target vertex

    int weight; // Edge weight



    Edge(int s, int t, int w) : source(s), target(t), weight(w) {}



    // Default constructor for MPI communication

    Edge() : source(-1), target(-1), weight(-1) {}

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



    // MPI related variables

    int mpi_rank;                     // Current process rank

    int mpi_size;                     // Total number of processes

    vector<int> partition_to_process; // Maps each partition to a process



    // OpenMP related variables

    int omp_num_threads; // Number of OpenMP threads per process



public:

    // Constructor to initialize the graph

    DynamicSSSP(int n, int source, bool adjust = false, int rank = 0, int size = 1) : num_vertices(n), source_vertex(source), num_partitions(0), adjust_indices(adjust),

                                                                                      mpi_rank(rank), mpi_size(size)

    {

        outgoing_edges.resize(n);

        incoming_edges.resize(n);

        distances.resize(n, numeric_limits<int>::max());

        parents.resize(n, -1);

        vertex_to_partition.resize(n, -1);



// Get the number of OpenMP threads available

#pragma omp parallel

        {

#pragma omp master

            {

                omp_num_threads = omp_get_num_threads();

            }

        }

    }



    // Add an edge to the graph

    void add_edge(int source, int target, int weight)

    {

        // Adjust indices if needed (convert from 1-indexed to 0-indexed)

        if (adjust_indices)

        {

            source--;

            target--;

        }



        if (source < 0 || source >= num_vertices || target < 0 || target >= num_vertices)

        {

            if (mpi_rank == 0)

            {

                cerr << "Edge vertices out of range: (" << source << "," << target << ")\n";

            }

            return;

        }



#pragma omp critical

        {

            outgoing_edges[source].push_back(Edge(source, target, weight));

            incoming_edges[target].push_back(Edge(source, target, weight));

        }

    }



    // Load the METIS partition file

    bool load_partitions(const string &filename)

    {

        // Only rank 0 reads the file and broadcasts the data

        vector<int> partition_data;



        if (mpi_rank == 0)

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



            partition_data = vertex_to_partition;

            cout << "Loaded " << num_partitions << " partitions for " << vertex_id << " vertices\n";

        }



        // Broadcast number of partitions

        MPI_Bcast(&num_partitions, 1, MPI_INT, 0, MPI_COMM_WORLD);



        // Broadcast partition assignments

        if (mpi_rank != 0)

        {

            partition_data.resize(num_vertices);

        }



        MPI_Bcast(partition_data.data(), num_vertices, MPI_INT, 0, MPI_COMM_WORLD);



        if (mpi_rank != 0)

        {

            vertex_to_partition = partition_data;

        }



        // Assign partitions to processes (round-robin)

        partition_to_process.resize(num_partitions);

        for (int p = 0; p < num_partitions; p++)

        {

            partition_to_process[p] = p % mpi_size;

        }



        // Print partition assignment

        if (mpi_rank == 0)

        {

            cout << "Partition to process assignment:" << endl;

            for (int p = 0; p < num_partitions; p++)

            {

                cout << "Partition " << p << " -> Process " << partition_to_process[p] << endl;

            }

            cout << "Each process will use " << omp_num_threads << " OpenMP threads" << endl;

        }



        return true;

    }



    // Initialize graph from a file in simple format

    // Format: first line contains number of vertices and edges

    // Following lines contain source, target, weight for each edge

    bool load_graph(const string &filename)

    {

        int n = 0;

        int m = 0;



        if (mpi_rank == 0)

        {

            ifstream infile(filename);

            if (!infile)

            {

                cerr << "Error opening graph file: " << filename << "\n";

                return false;

            }



            infile >> n >> m;



            if (n <= 0 || m < 0)

            {

                cerr << "Invalid graph dimensions\n";

                return false;

            }



            cout << "Loaded graph with " << n << " vertices and " << m << " edges\n";

        }



        // Broadcast vertices and edges count

        MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

        MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);



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



        // Each process reads and processes the whole graph for now

        // A more optimized approach would be to distribute the graph loading

        if (mpi_rank == 0)

        {

            ifstream infile(filename);

            infile >> n >> m; // Skip first line (already read)



            int u, v, w;

            int edge_count = 0;

            int error_count = 0;



            // Read edges in chunks for parallel processing

            vector<tuple<int, int, int>> edge_data;

            edge_data.reserve(m);



            while (edge_count < m && infile >> u >> v >> w)

            {

                // Validate edge before adding

                if (adjust_indices)

                {

                    if (u <= 0 || u > n || v <= 0 || v > n)

                    {

                        cerr << "Edge vertices out of range (before adjustment): (" << u << "," << v << ")\n";

                        error_count++;

                        continue;

                    }

                }

                else

                {

                    if (u < 0 || u >= n || v < 0 || v >= n)

                    {

                        cerr << "Edge vertices out of range: (" << u << "," << v << ")\n";

                        error_count++;

                        continue;

                    }

                }



                edge_data.push_back(make_tuple(u, v, w));

                edge_count++;



                // Print progress for large graphs

                if (edge_count % 100000 == 0)

                {

                    cout << "Read " << edge_count << " edges..." << endl;

                }

            }



            cout << "Processing " << edge_count << " edges with OpenMP parallelism..." << endl;



// Process edges in parallel with OpenMP

#pragma omp parallel

            {

#pragma omp single

                {

                    cout << "Using " << omp_get_num_threads() << " threads for edge processing" << endl;

                }



#pragma omp for schedule(dynamic, 1000)

                for (int i = 0; i < edge_data.size(); i++)

                {

                    int u = get<0>(edge_data[i]);

                    int v = get<1>(edge_data[i]);

                    int w = get<2>(edge_data[i]);



                    // Add edge locally - no need for broadcast yet

                    add_edge(u, v, w);

                }

            }



            // After all edges are added locally, broadcast the edge data to other ranks

            int edge_count_bcast = edge_data.size();

            MPI_Bcast(&edge_count_bcast, 1, MPI_INT, 0, MPI_COMM_WORLD);



            // Send each edge individually - could be optimized with custom MPI type

            for (auto &edge : edge_data)

            {

                int edge_info[3] = {get<0>(edge), get<1>(edge), get<2>(edge)};

                MPI_Bcast(edge_info, 3, MPI_INT, 0, MPI_COMM_WORLD);

            }



            if (edge_count < m)

            {

                cerr << "Warning: Only read " << edge_count << " edges out of expected " << m << endl;

            }



            if (error_count > 0)

            {

                cerr << "Warning: Encountered " << error_count << " edges with invalid vertices." << endl;

            }



            cout << "Successfully added " << edge_count << " edges to the graph" << endl;

            return edge_count > 0; // Return true if at least some edges were loaded

        }

        else

        {

            // Other processes receive edge information and add edges

            int edge_count_bcast;

            MPI_Bcast(&edge_count_bcast, 1, MPI_INT, 0, MPI_COMM_WORLD);



            vector<tuple<int, int, int>> edge_data;

            edge_data.reserve(edge_count_bcast);



            for (int i = 0; i < edge_count_bcast; i++)

            {

                int edge_info[3];

                MPI_Bcast(edge_info, 3, MPI_INT, 0, MPI_COMM_WORLD);

                edge_data.push_back(make_tuple(edge_info[0], edge_info[1], edge_info[2]));

            }



// Process received edges in parallel

#pragma omp parallel for schedule(dynamic, 1000)

            for (int i = 0; i < edge_data.size(); i++)

            {

                int u = get<0>(edge_data[i]);

                int v = get<1>(edge_data[i]);

                int w = get<2>(edge_data[i]);

                add_edge(u, v, w);

            }



            return edge_count_bcast > 0;

        }

    }



    // Process a single partition

    void process_partition(int partition_id)

    {

        // Only process if this partition is assigned to this process

        if (partition_to_process[partition_id] != mpi_rank)

        {

            return;

        }



        if (mpi_rank == 0)

        {

            cout << "\nProcess " << mpi_rank << " processing partition " << partition_id << endl;

        }



        // Find vertices in this partition

        vector<int> partition_vertices;

        for (int v = 0; v < num_vertices; v++)

        {

            if (vertex_to_partition[v] == partition_id)

            {

                partition_vertices.push_back(v);

            }

        }



        if (mpi_rank == 0)

        {

            cout << "Partition " << partition_id << " contains " << partition_vertices.size() << " vertices" << endl;

        }



// Process vertices in parallel with OpenMP

#pragma omp parallel for schedule(dynamic)

        for (int i = 0; i < partition_vertices.size(); i++)

        {

            int v = partition_vertices[i];

            process_vertex(v);

        }

    }



    // Process a single vertex's outgoing edges

    void process_vertex(int u)

    {

        // Skip if vertex distance is not yet determined

        if (distances[u] == numeric_limits<int>::max())

            return;



        // Calculate relaxation for all neighbors

        // This is thread-safe as we use atomic operations for distance updates

        for (const Edge &edge : outgoing_edges[u])

        {

            int v = edge.target;

            int weight = edge.weight;

            int new_distance = distances[u] + weight;



            if (new_distance < distances[v])

            {

// We need atomic operation here to avoid race conditions

#pragma omp critical

                {

                    if (new_distance < distances[v])

                    {

                        distances[v] = new_distance;

                        parents[v] = u;

                    }

                }

            }

        }

    }



    // Compute initial SSSP using partitioned Dijkstra's algorithm with MPI and OpenMP

    void compute_partitioned_sssp()

    {

        // Initialize distances

        fill(distances.begin(), distances.end(), numeric_limits<int>::max());

        fill(parents.begin(), parents.end(), -1);



        // Set source distance to 0

        distances[source_vertex] = 0;



        if (mpi_rank == 0)

        {

            cout << "Starting partitioned SSSP with source vertex " << source_vertex

                 << " using " << mpi_size << " MPI processes and "

                 << omp_num_threads << " OpenMP threads per process" << endl;

        }



        // Each process processes its assigned partitions

        for (int p = 0; p < num_partitions; p++)

        {

            // Each process works on partitions assigned to it

            process_partition(p);



            // Synchronize distances and parents after each partition

            sync_distances_and_parents();

        }



        // We may need to run multiple passes because changes in one partition

        // might affect vertices in earlier partitions

        bool changed = true;

        int pass = 1;

        const int MAX_PASSES = 3; // Limit number of passes to avoid infinite loops



        while (changed && pass < MAX_PASSES)

        {

            changed = false;



            if (mpi_rank == 0)

            {

                cout << "\nRunning additional pass " << pass << endl;

            }



            // Each process keeps its own copy of previous distances

            vector<int> old_distances = distances;



            for (int p = 0; p < num_partitions; p++)

            {

                // Each process works only on its assigned partitions

                process_partition(p);



                // Synchronize after each partition

                sync_distances_and_parents();

            }



            // Check if any distances changed globally - can use OpenMP to parallelize

            int local_changed = 0;



#pragma omp parallel

            {

                int thread_changed = 0;



#pragma omp for schedule(static)

                for (int v = 0; v < num_vertices; v++)

                {

                    if (old_distances[v] != distances[v])

                    {

                        thread_changed = 1;

                    }

                }



                // Combine results from all threads

                if (thread_changed)

                {

#pragma omp critical

                    local_changed = 1;

                }

            }

            // Reduce to determine if any process had changes

            int global_changed = 0;

            MPI_Allreduce(&local_changed, &global_changed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

            changed = (global_changed > 0);



            pass++;

        }

    }



    // Synchronize distances and parents across all processes

    void sync_distances_and_parents()

    {

        // Create temporary buffers for the reduction

        vector<int> min_distances(num_vertices, numeric_limits<int>::max());

        vector<int> best_parents(num_vertices, -1);



        // Find minimum distance for each vertex across all processes

        MPI_Allreduce(distances.data(), min_distances.data(), num_vertices, MPI_INT, MPI_MIN, MPI_COMM_WORLD);



// For each vertex, find the process with the minimum distance

// and broadcast the parent - can be optimized with custom MPI operation

#pragma omp parallel for schedule(static)

        for (int v = 0; v < num_vertices; v++)

        {

            // Prepare values for reduction

            int local_value = (distances[v] == min_distances[v] &&

                               distances[v] != numeric_limits<int>::max())

                                  ? (mpi_rank + 1)

                                  : 0;



            int owner_process = 0;

            MPI_Allreduce(&local_value, &owner_process, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);



            // If this process has the min distance, prepare to broadcast its parent value

            int parent_value = (owner_process == mpi_rank + 1) ? parents[v] : -1;



            // Collect the parent from the process with min distance

            MPI_Allreduce(&parent_value, &best_parents[v], 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

        }



// Update local copies with synchronized values

#pragma omp parallel for schedule(static)

        for (int v = 0; v < num_vertices; v++)

        {

            distances[v] = min_distances[v];

            if (best_parents[v] != -1)

            {

                parents[v] = best_parents[v];

            }

        }

    }



    // Handle edge weight decrease

    void handle_edge_decrease(int u, int v, int new_weight)

    {

        // Adjust indices if needed

        if (adjust_indices)

        {

            u--;

            v--;

        }



        // Find and update the edge

        bool edge_found = false;



#pragma omp parallel for

        for (int i = 0; i < outgoing_edges[u].size(); i++)

        {

            Edge &edge = outgoing_edges[u][i];

            if (edge.target == v)

            {

#pragma omp critical

                {

                    edge.weight = new_weight;

                    edge_found = true;

                }

            }

        }



        // Update corresponding incoming edge

        if (edge_found)

        {

#pragma omp parallel for

            for (int i = 0; i < incoming_edges[v].size(); i++)

            {

                Edge &in_edge = incoming_edges[v][i];

                if (in_edge.source == u)

                {

#pragma omp critical

                    {

                        in_edge.weight = new_weight;

                    }

                }

            }



            // Synchronize edge update across all processes

            MPI_Barrier(MPI_COMM_WORLD);



            // Check if this creates a shorter path

            if (distances[u] != numeric_limits<int>::max() &&

                distances[u] + new_weight < distances[v])

            {

                // Update distance and parent

                distances[v] = distances[u] + new_weight;

                parents[v] = u;



                // Process all partitions starting from the affected vertex's partition

                int v_partition = vertex_to_partition[v];

                for (int p = v_partition; p < num_partitions; p++)

                {

                    process_partition(p);

                    sync_distances_and_parents();

                }

            }

        }



        if (!edge_found && mpi_rank == 0)

        {

            cerr << "Edge (" << u << "," << v << ") not found in graph\n";

        }

    }



    // Handle edge weight increase

    void handle_edge_increase(int u, int v, int new_weight)

    {

        // Adjust indices if needed

        if (adjust_indices)

        {

            u--;

            v--;

        }



        // Find and update the edge

        bool edge_found = false;

        bool requires_recompute = false;



#pragma omp parallel for

        for (int i = 0; i < outgoing_edges[u].size(); i++)

        {

            Edge &edge = outgoing_edges[u][i];

            if (edge.target == v)

            {

#pragma omp critical

                {

                    edge.weight = new_weight;

                    edge_found = true;

                    // Check if this edge was part of the shortest path

                    if (parents[v] == u)

                    {

                        requires_recompute = true;

                    }

                }

            }

        }



        // Update corresponding incoming edge

        if (edge_found)

        {

#pragma omp parallel for

            for (int i = 0; i < incoming_edges[v].size(); i++)

            {

                Edge &in_edge = incoming_edges[v][i];

                if (in_edge.source == u)

                {

#pragma omp critical

                    {

                        in_edge.weight = new_weight;

                    }

                }

            }



            // Synchronize edge update across all processes

            MPI_Barrier(MPI_COMM_WORLD);



            // If the edge was part of the shortest path, recompute

            // Reduce to make sure all processes agree on the need to recompute

            int local_recompute = requires_recompute ? 1 : 0;

            int global_recompute = 0;



            MPI_Allreduce(&local_recompute, &global_recompute, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);



            if (global_recompute)

            {

                if (mpi_rank == 0)

                {

                    cout << "Edge weight increase requires recomputation of SSSP" << endl;

                }

                compute_partitioned_sssp();

            }

        }



        if (!edge_found && mpi_rank == 0)

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



// OpenMP parallel search for the edge

#pragma omp parallel

        {

#pragma omp for nowait

            for (int i = 0; i < outgoing_edges[u_idx].size(); i++)

            {

                if (outgoing_edges[u_idx][i].target == v_idx)

                {

#pragma omp critical

                    {

                        current_weight = outgoing_edges[u_idx][i].weight;

                    }

                }

            }

        }



        // Gather current weight from all processes

        int global_weight = -1;

        MPI_Allreduce(&current_weight, &global_weight, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

        current_weight = global_weight;



        if (current_weight == -1)

        {

            if (mpi_rank == 0)

            {

                cerr << "Edge (" << u << "," << v << ") not found in graph\n";

            }

            return;

        }



        if (mpi_rank == 0)

        {

            cout << "Updating edge (" << u << "," << v << ") from weight "

                 << current_weight << " to " << new_weight << endl;

        }



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

            if (mpi_rank == 0)

            {

                cout << "No weight change, skipping update\n";

            }

        }

    }



    // Print the current shortest path distances (only rank 0)

    void print_distances(int limit = 20)

    {

        if (mpi_rank != 0)

            return;



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



        if (count < num_vertices)

        {

            cout << "... (showing " << count << " out of " << num_vertices << " vertices)\n";

        }



        cout << "----------------------------------------\n";

    }



    // Print summary statistics by partition (only rank 0)

    void print_partition_stats()

    {

        if (mpi_rank != 0)

            return;



        cout << "\nPartition Statistics:\n";

        cout << "----------------------------------------\n";

        cout << "Partition | Vertices | Avg Distance | Process\n";

        cout << "----------------------------------------\n";



        // Compute statistics in parallel

        vector<int> partition_counts(num_partitions, 0);

        vector<long long> partition_sums(num_partitions, 0);



#pragma omp parallel for schedule(dynamic)

        for (int p = 0; p < num_partitions; p++)

        {

            for (int v = 0; v < num_vertices; v++)

            {

                if (vertex_to_partition[v] == p && distances[v] != numeric_limits<int>::max())

                {

#pragma omp atomic

                    partition_counts[p]++;



#pragma omp atomic

                    partition_sums[p] += distances[v];

                }

            }

        }



        // Print the results

        for (int p = 0; p < num_partitions; p++)

        {

            cout << setw(9) << p << " | ";

            cout << setw(8) << partition_counts[p] << " | ";



            if (partition_counts[p] > 0)

                cout << setw(12) << fixed << setprecision(2) << (double)partition_sums[p] / partition_counts[p] << " | ";

            else

                cout << setw(12) << "N/A" << " | ";



            cout << setw(7) << partition_to_process[p] << "\n";

        }

        cout << "----------------------------------------\n";

    }



    // Print the shortest path tree structure (only rank 0)

    void print_tree_structure(int max_depth = 3, int max_children = 5)

    {

        if (mpi_rank != 0)

            return;



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

            if (depth > max_depth)

            {

                print_indent(depth);

                cout << "... (max depth reached)\n";

                return;

            }



            print_indent(depth);

            cout << "└── " << vertex;



            if (distances[vertex] == numeric_limits<int>::max())

                cout << " (INF)" << endl;

            else

                cout << " (dist=" << distances[vertex] << ")" << endl;



            // Find children in shortest path tree

            vector<int> children;

            for (int v = 0; v < num_vertices; v++)

            {

                if (parents[v] == vertex && v != vertex)

                {

                    children.push_back(v);

                }

            }

            sort(children.begin(), children.end(), [this](int a, int b)

                 { return distances[a] < distances[b]; });



            int child_count = 0;

            for (int child : children)

            {

                if (child_count < max_children)

                {

                    print_subtree(child, depth + 1);

                    child_count++;

                }

                else

                {

                    print_indent(depth + 1);

                    cout << "... (" << (children.size() - max_children) << " more children)\n";

                    break;

                }

            }

        };



        print_subtree(source_vertex, 0);

        cout << "----------------------------------------\n";

    }



    // Add getter for num_vertices

    int get_num_vertices() const

    {

        return num_vertices;

    }



    // Add getter for mpi_rank

    int get_mpi_rank() const

    {

        return mpi_rank;

    }

};

int main()

{

    // Initialize MPI

    int argc = 0;

    char **argv = nullptr;

    MPI_Init(nullptr, nullptr);



    int mpi_rank, mpi_size;

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);



    if (mpi_rank == 0)

    {

        cout << "\n\n----------------------------------------\n";

        cout << "----------------------------------------\n";

        cout << "MPI + OpenMP Implementation \n";

        cout << "----------------------------------------\n";

        cout << "----------------------------------------\n";

    }



    // Static input files and source vertex

    string graph_filename = "data.graph";

    string partition_filename = "data.graph.part.3";

    int source_vertex = 0; // You can change this to any valid source vertex



    if (mpi_rank == 0)

    {

        cout << "Graph file: " << graph_filename << "\n";

        cout << "Partition file: " << partition_filename << "\n";

        cout << "Source vertex: " << source_vertex << "\n";

    }



    // Create the dynamic SSSP object

    DynamicSSSP *sssp = new DynamicSSSP(0, source_vertex, true, mpi_rank, mpi_size);



    // Load the graph

    auto start = chrono::high_resolution_clock::now();

    bool graph_loaded = sssp->load_graph(graph_filename);

    auto end = chrono::high_resolution_clock::now();

    chrono::duration<double> elapsed = end - start;



    if (!graph_loaded)

    {

        if (mpi_rank == 0)

        {

            cerr << "Failed to load graph. Exiting.\n";

        }

        delete sssp;

        MPI_Finalize();

        return 1;

    }



    if (mpi_rank == 0)

    {

        cout << fixed << setprecision(6);

        cout << "Loading graph from file: " << elapsed.count() << " seconds\n";

    }



    // Load the partitions

    start = chrono::high_resolution_clock::now();

    bool partitions_loaded = sssp->load_partitions(partition_filename);

    end = chrono::high_resolution_clock::now();

    elapsed = end - start;



    if (!partitions_loaded)

    {

        if (mpi_rank == 0)

        {

            cerr << "Failed to load partitions. Exiting.\n";

        }

        delete sssp;

        MPI_Finalize();

        return 1;

    }



    if (mpi_rank == 0)

    {

        cout << fixed << setprecision(6);

        cout << "Loading partitions from file: " << elapsed.count() << " seconds\n";

    }



    // Compute initial shortest paths

    start = chrono::high_resolution_clock::now();

    sssp->compute_partitioned_sssp();

    end = chrono::high_resolution_clock::now();

    elapsed = end - start;



    if (mpi_rank == 0)

    {

        cout << fixed << setprecision(6);

        cout << "Initial SSSP computation time: " << elapsed.count() << " seconds\n";

    }



    sssp->print_distances();

    sssp->print_partition_stats();

    sssp->print_tree_structure();



    // Simulate some edge updates

    if (mpi_rank == 0)

    {

        cout << "\n=== Simulating Edge Updates ===\n";

    }



    // Decrease weight example

    if (sssp->get_num_vertices() > 1)

    {

        if (mpi_rank == 0)

        {

            cout << "\n[1] Decreasing edge weight:\n";

        }

        sssp->update_edge_weight(1, 2, 1);

        sssp->print_distances();

        sssp->print_tree_structure();

    }



    // Increase weight example

    if (sssp->get_num_vertices() > 1)

    {

        if (mpi_rank == 0)

        {

            cout << "\n[2] Increasing edge weight:\n";

        }

        sssp->update_edge_weight(1, 2, 10);

        sssp->print_distances();

        sssp->print_tree_structure();

    }



    delete sssp;



    // Finalize MPI

    MPI_Finalize();



    return 0;

}