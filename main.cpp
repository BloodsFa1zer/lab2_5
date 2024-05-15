// Block 0: 1, 3
// Block 1: 5, 6
// Block 2: 11
// Block 3: 14
// Block 4: 18
// Block 5: 19
// Block 6: 21


#include <iostream>
#include <vector>
#include <bitset>
#include <cstdlib> // For rand() function
#include <ctime>   // For seeding srand()
#include <algorithm>
#include <queue>
#include <set>
#include <stack>
#include <thread>

using namespace std;
using namespace chrono;

#define INF 1e9

const int MAX_VERTICES = 32; // Maximum number of vertices for the graph


void loadingAnimation() {
    for (int i = 0; i < 3; ++i) {
        cout << "Loading.";
        cout.flush();
        this_thread::sleep_for(milliseconds(250));
        cout << "\r";

        cout << "Loading..";
        cout.flush();
        this_thread::sleep_for(milliseconds(250));
        cout << "\r";

        cout << "Loading...";
        cout.flush();
        this_thread::sleep_for(milliseconds(250));
        cout << "\r";
    }
}

int getIntInput(const string& prompt, int minValue, int maxValue) {
    int value;
    bool inputIsValid;
    do {
        cout << prompt;
        if (cin >> value && value >= minValue && value <= maxValue) {
            inputIsValid = true;
        } else {
            cerr << "\n### Invalid input. Please enter a valid integer between "
                 << minValue << " and " << maxValue << ". ###\n\n";
            cin.clear();
            cin.ignore(numeric_limits<streamsize>::max(), '\n');
            inputIsValid = false;
        }
    } while (!inputIsValid);
    return value;
}

// Define a structure to represent an edge
struct Edge {
    int src, dest; // Source and destination vertices of the edge
    int weight; // Weight of the edge (not used in bit vector graph)

    Edge(int src, int dest, int weight) : src(src), dest(dest), weight(weight) {}
};

// Function to perform find operation in disjoint set
int find(vector<int>& parent, int vertex) {
    if (parent[vertex] != vertex)
        parent[vertex] = find(parent, parent[vertex]);
    return parent[vertex];
}

// Function to perform union operation in disjoint set
void unionSets(vector<int>& parent, int x, int y) {
    int rootX = find(parent, x);
    int rootY = find(parent, y);
    parent[rootX] = rootY;
}

// Graph data structure based on a bit vector adjacency list
struct BitVectorGraph {
    bitset<MAX_VERTICES> adjList[MAX_VERTICES]; // Adjacency list using bit vectors
    int V; // Number of vertices

    // Constructor to initialize the graph with V vertices
    BitVectorGraph(int V) : V(V) {}

// Function to add an edge between vertices u and v
    // Function to add an edge between vertices u and v
    void addEdge(int u, int v) {
        // Check if the vertices are valid
        if (u < 0 || u >= V || v < 0 || v >= V) {
            cout << "Error: Invalid vertices." << endl;
            return;
        }

        // Add the edge between u and v
        adjList[u].set(v);
        adjList[v].set(u); // For undirected graphs, set the bit for the reverse direction as well
    }


    // Function to print the adjacency list
    void printGraph() {
        cout << "Bit Vector Adjacency List:" << endl;
        for (int i = 0; i < V; ++i) {
            cout << "Vertex " << i << " is adjacent to: ";
            for (int j = 0; j < V; ++j) {
                if (adjList[i].test(j)) {
                    cout << j << " ";
                }
            }
            cout << endl;
        }
    }

    // Function to check if the graph is directed
    bool isDirectedBitVector() {
        for (int i = 0; i < V; ++i) {
            for (int j = i + 1; j < V; ++j) {
                if (adjList[i].test(j) != adjList[j].test(i)) {
                    return true;
                }
            }
        }
        return false;
    }


    // Function to create a random graph with specified number of edges
    void createRandomGraph(int E, double maxWeight) {
        srand(time(NULL)); // Seed the random number generator

        int edgesAdded = 0;
        while (edgesAdded < E) {
            int u = rand() % V; // Randomly select a source vertex
            int v = rand() % V; // Randomly select a destination vertex
            if (u != v && !adjList[u].test(v)) { // Ensure the edge does not already exist and is not a self-loop
                double rawWeight = static_cast<double>(rand()) / (static_cast<double>(RAND_MAX / maxWeight));
                double weight = round(rawWeight * 10) / 10.0; // Round to one decimal place
                addEdge(u, v);
                ++edgesAdded;
            }
        }
    }


// Depth-First Search (DFS) function to explore connected components
    void DFS(int v, vector<bool>& visited) {
        visited[v] = true; // Mark the current vertex as visited
        cout << v << " "; // Print the current vertex

        // Traverse all adjacent vertices of v
        for (int i = 0; i < V; ++i) {
            // If an edge exists from v to i and i is not visited, recursively call DFS
            if (adjList[v][i] && !visited[i]) {
                DFS(i, visited);
            }
        }
    }

    // Function to find and print connected components
    void findConnectedComponents() {
        vector<bool> visited(V, false); // Initialize all vertices as not visited

        cout << "Connected Components:" << endl;
        for (int v = 0; v < V; ++v) {
            if (!visited[v]) {
                DFS(v, visited); // Call DFS for each unvisited vertex
                cout << endl; // Separate components by newline
            }
        }
    }

    // Function to perform DFS traversal on the graph with neighbors traversed in order of increasing edge weights
    void dfsIncreasingEdgeWeights(int v, vector<bool>& visited) {
        visited[v] = true; // Mark current vertex as visited
        cout << v << " "; // Print the current vertex

        // Create a set to store adjacent vertices in increasing order of edge weights
        set<int> adjVertices;

        // Insert adjacent vertices of vertex v into the set in increasing order of edge weights
        for (int i = 0; i < V; ++i) {
            if (adjList[v].test(i) && !visited[i]) {
                adjVertices.insert(i);
            }
        }

        // Traverse the adjacent vertices in increasing order of edge weights
        for (int vertex : adjVertices) {
            if (!visited[vertex]) {
                dfsIncreasingEdgeWeights(vertex, visited); // Recursively visit unvisited adjacent vertices
            }
        }
    }


    // Dijkstra's algorithm to find the shortest distances from a source vertex to all other vertices
    void dijkstra(int source) {
        vector<double> dist(V, INF); // Initialize distance array with infinity
        queue<int> Q;

        dist[source] = 0;
        Q.push(source);

        while (!Q.empty()) {
            int u = Q.front();
            Q.pop();

            for (int v = 0; v < V; ++v) {
                if (adjList[u].test(v)) {
                    if (dist[u] + 1 < dist[v]) { // Assuming unweighted graph, hence weight is 1
                        dist[v] = dist[u] + 1;
                        Q.push(v);
                    }
                }
            }
        }

        // Print the shortest distances from the source vertex to all other vertices
        cout << "Vertex \t Distance from Source\n";
        for (int i = 0; i < V; ++i) {
            cout << i << " \t ";
            (dist[i] != INF ? cout << dist[i] : cout << "Path doesn't exist");
            cout << endl;
        }
    }

    void DFS(int v, vector<bool>& visited, stack<int>& stk) {
        visited[v] = true;

        for (int i = 0; i < V; ++i) {
            if (adjList[v].test(i) && !visited[i]) {
                DFS(i, visited, stk);
            }
        }

        stk.push(v);
    }

    // Function to perform topological sorting
    void topologicalSort() {
        stack<int> stk;
        vector<bool> visited(V, false);

        for (int i = 0; i < V; ++i) {
            if (!visited[i]) {
                DFS(i, visited, stk);
            }
        }

        cout << "Topological Sorting: ";
        while (!stk.empty()) {
            cout << stk.top() << " ";
            stk.pop();
        }
        cout << endl;
    }

};


// DFS traversal function for BitVectorGraph
void dfsBit(const BitVectorGraph& graph, int v, vector<bool>& visited) {
    visited[v] = true; // Mark current vertex as visited
    // Traverse all adjacent vertices of vertex v
    for (int i = 0; i < graph.V; ++i) {
        if (graph.adjList[v].test(i) && !visited[i]) {
            dfsBit(graph, i, visited); // Recursively visit unvisited adjacent vertices
        }
    }
}

// Function to check if the graph is connected using Depth-First Search (DFS)
bool isGraphConnectedBit(const BitVectorGraph& graph) {
    vector<bool> visited(graph.V, false); // Mark all vertices as not visited
    dfsBit(graph, 0, visited); // Start DFS traversal from vertex 0

    // Check if all vertices were visited during the DFS traversal
    for (bool v : visited) {
        if (!v) {
            return false; // If any vertex was not visited, the graph is disconnected
        }
    }
    return true; // All vertices were visited, graph is connected
}

// Function to perform DFS and construct a spanning tree based on depth-first search
void dfsSpanningTreeBit(int u, const BitVectorGraph& graph, vector<bool>& visited, vector<int>& parent) {
    visited[u] = true;
    for (int v = 0; v < graph.V; ++v) {
        if (graph.adjList[u].test(v) && !visited[v]) {
            parent[v] = u; // Record the parent of vertex v
            dfsSpanningTreeBit(v, graph, visited, parent);
        }
    }
}

// Function to construct a spanning tree based on depth-first search
void buildSpanningTreeDFSBit(const BitVectorGraph& graph) {
    int V = graph.V;
    vector<bool> visited(V, false);
    vector<int> parent(V, -1); // Initialize parent array to store the spanning tree edges
    int totalWeight = 0; // Total weight of the spanning tree

    for (int u = 0; u < V; ++u) {
        if (!visited[u]) {
            dfsSpanningTreeBit(u, graph, visited, parent);
        }
    }

    // Calculate total weight of the spanning tree
    for (int v = 0; v < V; ++v) {
        if (parent[v] != -1) {
            // Since the bit vector graph does not have weights, we cannot calculate the total weight
            // We'll simply count the number of edges in the spanning tree
            totalWeight += 1;
        }
    }

    // Output the total weight of the spanning tree
    cout << "Total number of edges in the spanning tree: " << totalWeight << endl;
}

// Function to perform Kruskal's algorithm for the bit vector graph
void kruskalBit(const BitVectorGraph& graph, bool benchmark = false) {
    int V = graph.V;
    vector<int> parent(V);
    for (int i = 0; i < V; ++i)
        parent[i] = i;

    vector<Edge> edges;
    for (int u = 0; u < V; ++u) {
        for (int v = 0; v < V; ++v) {
            if (graph.adjList[u].test(v)) {
                edges.emplace_back(u, v, 0); // Weight is not used for bit vector graph
            }
        }
    }

    double totalWeight = 0; // Total weight of the minimum spanning tree

    // Sort edges by their source and destination vertices
    sort(edges.begin(), edges.end(), [](const Edge& a, const Edge& b) {
        if (a.src != b.src) return a.src < b.src;
        return a.dest < b.dest;
    });

    for (const Edge& edge : edges) {
        int rootX = find(parent, edge.src);
        int rootY = find(parent, edge.dest);
        if (rootX != rootY) {
            if (!benchmark) cout << edge.src << " -> " << edge.dest << endl;
            totalWeight += 1; // Increase the total weight by 1 for each edge
            unionSets(parent, edge.src, edge.dest);
        }
    }

    if (!benchmark) cout << "Total weight of MST: " << totalWeight << endl;
}




// Graph data structure based on adjacency matrix
struct AdjMatrixGraph {
    int V; // Number of vertices
    vector<vector<int>> adjMatrix; // Adjacency matrix

    // Constructor to initialize the graph with V vertices
    AdjMatrixGraph(int V) : V(V) {
        // Initialize the adjacency matrix with all zeros
        adjMatrix.resize(V, vector<int>(V, 0));
    }


    void addEdge(int u, int v, int weight) {
        // Check if the vertices are valid
        if (u < 0 || u >= V || v < 0 || v >= V) {
            cout << "Error: Invalid vertices." << endl;
            return;
        }

        // Check if the edge already exists
        if (adjMatrix[u][v] != 0 || adjMatrix[v][u] != 0) {
            cout << "Error: Edge already exists." << endl;
            return;
        }

        // Check for self-loop
        if (u == v) {
            cout << "Error: Self-loops are not allowed." << endl;
            return;
        }

        // Add the edge between u and v with weight
        adjMatrix[u][v] = weight;
        adjMatrix[v][u] = weight; // For undirected graphs, add the edge in reverse direction too
    }

    // Function to check if the graph is directed
    bool isDirectedAdjMatrix() {
        for (int i = 0; i < V; ++i) {
            for (int j = i + 1; j < V; ++j) {
                if (adjMatrix[i][j] != adjMatrix[j][i]) {
                    return true;
                }
            }
        }
        return false;
    }


    // Function to print the adjacency matrix
    void printGraph() {
        cout << "Adjacency Matrix:" << endl;
        for (int i = 0; i < V; ++i) {
            for (int j = 0; j < V; ++j) {
                cout << adjMatrix[i][j] << " ";
            }
            cout << endl;
        }
    }

    void createRandomGraph(int E, double maxWeight) {
        srand(time(NULL));

        int edgesAdded = 0;
        while (edgesAdded < E) {
            int u = rand() % V; // Randomly select a source vertex
            int v = rand() % V; // Randomly select a destination vertex
            if (u != v && !adjMatrix[u][v]) { // Ensure the edge does not already exist and is not a self-loop
                double rawWeight = static_cast<double>(rand()) / (static_cast<double>(RAND_MAX / maxWeight));
                double weight = round(rawWeight * 10) / 10.0; // Round to one decimal place
                addEdge(u, v, weight);
                ++edgesAdded;
            }
        }
    }


    // Depth-First Search (DFS) function to explore connected components
    void DFS(int v, vector<bool>& visited) {
        visited[v] = true; // Mark the current vertex as visited
        cout << v << " "; // Print the current vertex

        // Traverse all adjacent vertices of v
        for (int i = 0; i < V; ++i) {
            // If an edge exists from v to i and i is not visited, recursively call DFS
            if (adjMatrix[v][i] && !visited[i]) {
                DFS(i, visited);
            }
        }
    }

    // Function to find and print connected components
    void findConnectedComponents() {
        vector<bool> visited(V, false); // Initialize all vertices as not visited

        cout << "Connected Components:" << endl;
        for (int v = 0; v < V; ++v) {
            if (!visited[v]) {
                DFS(v, visited); // Call DFS for each unvisited vertex
                cout << endl; // Separate components by newline
            }
        }
    }





    // Function to perform DFS traversal on the graph with neighbors traversed in order of increasing edge weights
    void dfsIncreasingEdgeWeights(int v, vector<bool>& visited) {
        visited[v] = true; // Mark current vertex as visited
        cout << v << " "; // Print the current vertex

        // Create a vector to store adjacent vertices in increasing order of edge weights
        vector<int> adjVertices;

        // Insert adjacent vertices of vertex v into the vector in increasing order of edge weights
        for (int i = 0; i < V; ++i) {
            if (adjMatrix[v][i] && !visited[i]) {
                adjVertices.push_back(i);
            }
        }

        // Sort the adjacent vertices based on edge weights
        sort(adjVertices.begin(), adjVertices.end(), [&](int a, int b) {
            return adjMatrix[v][a] < adjMatrix[v][b];
        });

        // Traverse the adjacent vertices in increasing order of edge weights
        for (int vertex : adjVertices) {
            if (!visited[vertex]) {
                dfsIncreasingEdgeWeights(vertex, visited); // Recursively visit unvisited adjacent vertices
            }
        }
    }


    // Dijkstra's algorithm to find the shortest distances from a source vertex to all other vertices
    void dijkstra(int source) {
        vector<double> dist(V, INF); // Initialize distance array with infinity
        queue<int> Q;

        dist[source] = 0;
        Q.push(source);

        while (!Q.empty()) {
            int u = Q.front();
            Q.pop();

            for (int v = 0; v < V; ++v) {
                if (adjMatrix[u][v] && dist[u] != INF && dist[u] + adjMatrix[u][v] < dist[v]) {
                    dist[v] = dist[u] + adjMatrix[u][v];
                    Q.push(v);
                }
            }
        }

        // Print the shortest distances from the source vertex to all other vertices
        cout << "Vertex \t Distance from Source\n";
        for (int i = 0; i < V; ++i) {
            cout << i << " \t ";
            (dist[i] != INF ? cout << dist[i] : cout << "Path doesn't exist");
            cout << endl;
        }
    }

    // Utility function for topological sorting
    void topologicalSortUtil(int v, vector<bool>& visited, stack<int>& stack) {
        visited[v] = true;

        // Recur for all adjacent vertices
        for (int i = 0; i < V; ++i) {
            if (adjMatrix[v][i] && !visited[i]) {
                topologicalSortUtil(i, visited, stack);
            }
        }

        // Push the current vertex to the stack
        stack.push(v);
    }

    // Function to perform topological sorting
    void topologicalSort() {
        // Mark all the vertices as not visited
        vector<bool> visited(V, false);
        stack<int> stack;

        // Call the utility function for each vertex
        for (int i = 0; i < V; ++i) {
            if (!visited[i]) {
                topologicalSortUtil(i, visited, stack);
            }
        }

        // Print the contents of the stack
        cout << "Topological Sort: ";
        while (!stack.empty()) {
            cout << stack.top() << " ";
            stack.pop();
        }
        cout << endl;
    }
};

// DFS traversal function for AdjMatrixGraph
void dfsAdj(const AdjMatrixGraph& graph, int v, vector<bool>& visited) {
    visited[v] = true; // Mark current vertex as visited
    // Traverse all adjacent vertices of vertex v
    for (int i = 0; i < graph.V; ++i) {
        if (graph.adjMatrix[v][i] && !visited[i]) {
            dfsAdj(graph, i, visited); // Recursively visit unvisited adjacent vertices
        }
    }
}

// Function to check if the graph is connected using Depth-First Search (DFS)
bool isGraphConnectedAdj(const AdjMatrixGraph& graph) {
    vector<bool> visited(graph.V, false); // Mark all vertices as not visited
    dfsAdj(graph, 0, visited); // Start DFS traversal from vertex 0

    // Check if all vertices were visited during the DFS traversal
    for (bool v : visited) {
        if (!v) {
            return false; // If any vertex was not visited, the graph is disconnected
        }
    }
    return true; // All vertices were visited, graph is connected
}

// Function to perform DFS and construct a spanning tree based on depth-first search
void dfsSpanningTreeAdj(int u, vector<vector<int>>& adjMatrix, vector<bool>& visited, vector<int>& parent) {
    visited[u] = true;
    for (int v = 0; v < adjMatrix[u].size(); ++v) {
        if (adjMatrix[u][v] != 0 && !visited[v]) {
            parent[v] = u; // Record the parent of vertex v
            dfsSpanningTreeAdj(v, adjMatrix, visited, parent);
        }
    }
}

// Function to construct a spanning tree based on depth-first search
void buildSpanningTreeDFSAdj(AdjMatrixGraph& graph) {
    int V = graph.V;
    vector<bool> visited(V, false);
    vector<int> parent(V, -1); // Initialize parent array to store the spanning tree edges
    int totalWeight = 0; // Total weight of the spanning tree

    for (int u = 0; u < V; ++u) {
        if (!visited[u]) {
            dfsSpanningTreeAdj(u, graph.adjMatrix, visited, parent);
        }
    }

    // Calculate total weight of the spanning tree
    for (int v = 0; v < V; ++v) {
        if (parent[v] != -1) {
            totalWeight += graph.adjMatrix[parent[v]][v];
        }
    }

    // Output the total weight of the spanning tree
    cout << "Total weight of the spanning tree: " << totalWeight << endl;
}

// Function to perform Kruskal's algorithm for the adjacency matrix graph
void kruskalAdj(const AdjMatrixGraph& graph, bool benchmark = false) {
    int V = graph.V;
    vector<int> parent(V);
    for (int i = 0; i < V; ++i)
        parent[i] = i;

    vector<Edge> edges;
    for (int u = 0; u < V; ++u) {
        for (int v = u + 1; v < V; ++v) { // Only traverse the upper triangular part of the matrix
            if (graph.adjMatrix[u][v] != 0) {
                edges.emplace_back(u, v, graph.adjMatrix[u][v]);
            }
        }
    }

    double totalWeight = 0; // Total weight of the minimum spanning tree

    // Sort edges by their weights
    sort(edges.begin(), edges.end(), [](const Edge& a, const Edge& b) {
        return a.weight < b.weight;
    });

    for (const Edge& edge : edges) {
        int rootX = find(parent, edge.src);
        int rootY = find(parent, edge.dest);
        if (rootX != rootY) {
            if (!benchmark) cout << edge.src << " -> " << edge.dest << " (" << edge.weight << ")" << endl;
            totalWeight += edge.weight;
            unionSets(parent, edge.src, edge.dest);
        }
    }

    if (!benchmark) cout << "Total weight of MST: " << totalWeight << endl;
}


void userManagementMode();
void userManagementModeBitGraph();
void userManagementModeAdjGraph();
void demonstrationMode();
void benchmarkMode();

int main() {
    srand(time(nullptr));
    int mode;
    mode = getIntInput("Select mode:\n"
                       "1. User management mode\n"
                       "2. Demonstration mode \n"
                       "3. Benchmark mode\n", 1, 3);

    switch (mode) {
        case 1:
            userManagementMode();
            break;
        case 2:
            demonstrationMode();
            break;
        case 3:
            benchmarkMode();
    }




//    int V = 5; // Number of vertices in the graph
//    int E = 7; // Number of edges in the graph (for random graph creation)
//
//
//    BitVectorGraph bitVectorGraph(V);
//    bitVectorGraph.createRandomGraph(7, 10);
//    bitVectorGraph.printGraph();
//
//// Check if the bit vector graph is directed
//    cout << "Bit Vector Graph: ";
//    if (bitVectorGraph.isDirectedBitVector()) {
//        // Do something if the graph is directed
//    } else {
//        // Do something if the graph is undirected
//    }
//
//
//    // Check if the BitVectorGraph is connected
//    cout << "BitVectorGraph is " << (isGraphConnectedBit(bitVectorGraph) ? "connected" : "disconnected") << endl;
//
//    vector<bool> visitedBit(bitVectorGraph.V, false);
//    cout << "DFS Arbitrary Order (Bit Vector Graph): ";
//    bitVectorGraph.dfsArbitraryOrder(0, visitedBit);
//    cout << endl;
//
//    fill(visitedBit.begin(), visitedBit.end(), false); // Reset visited array
//    cout << "DFS Increasing Edge Weights (Bit Vector Graph): ";
//    bitVectorGraph.dfsIncreasingEdgeWeights(0, visitedBit);
//    cout << endl;
//
//    // Find the shortest paths from vertex 0
//    cout << "Shortest paths using BitVectorGraph:" << endl;
//    bitVectorGraph.dijkstra(0);
//    cout << endl;
//
//
//    cout << "Topological Sorting:" << endl;
//    bitVectorGraph.topologicalSort();
//
//    // Build and display the spanning tree using DFS
//    buildSpanningTreeDFSBit(bitVectorGraph);
//
//
//
//
//    // Create and manipulate the adjacency matrix graph
//    AdjMatrixGraph adjMatrixGraph(V);
//    adjMatrixGraph.createRandomGraph(7, 10);
//    adjMatrixGraph.printGraph();
//
//    // Check if the adjacency matrix graph is directed
//    cout << "Adjacency Matrix Graph: ";
//    if (adjMatrixGraph.isDirectedAdjMatrix()) {
//        // Do something if the graph is directed
//    } else {
//        // Do something if the graph is undirected
//    }
//
//    // Check if the AdjMatrixGraph is connected
//    cout << "AdjMatrixGraph is " << (isGraphConnectedAdj(adjMatrixGraph) ? "connected" : "disconnected") << endl;
//
//    vector<bool> visitedMatrix(adjMatrixGraph.V, false);
//    cout << "DFS Arbitrary Order (Adjacency Matrix Graph): ";
//    adjMatrixGraph.dfsArbitraryOrder(0, visitedMatrix);
//    cout << endl;
//
//    fill(visitedMatrix.begin(), visitedMatrix.end(), false); // Reset visited array
//    cout << "DFS Increasing Edge Weights (Adjacency Matrix Graph): ";
//    adjMatrixGraph.dfsIncreasingEdgeWeights(0, visitedMatrix);
//    cout << endl;
//
//    // Find the shortest paths from vertex 0
//    cout << "Shortest paths using AdjMatrixGraph:" << endl;
//    adjMatrixGraph.dijkstra(0);
//
//
//    cout << "Topological Sorting:" << endl;
//    adjMatrixGraph.topologicalSort();
//
//    // Build and display the spanning tree using DFS
//    buildSpanningTreeDFSAdj(adjMatrixGraph);
//
//    // Apply Kruskal's algorithm to find the minimum spanning tree
//    kruskalAdj(adjMatrixGraph);

    return 0;
}

//void printMenuUserManagement() {
//    cout << "\n====== Graph Operations Menu ======" << endl;
//    cout << "1. Create a graph" << endl;
//    cout << "2. Add an edge" << endl;
//    cout << "3. Perform algorithm" << endl;
//    cout << "4. Print graph details" << endl;
//    cout << "5. Exit" << endl;
//    cout << "Enter your choice: ";
//}

const string menu = "Select mode:\n"
                    "1. Create a graph\n"
                    "2. Create a random graph\n"
                    "3. Add an edge\n"
                    "4. Check if graph is directed\n"
                    "5. Print graph\n"
                    "6. Check if graph is connected\n"
                    "7. Depth-First Search (DFS) - Find Connected Components\n"
                    "8. Depth-First Search (DFS) - Increasing Edge Weights\n"
                    "9. Find the shortest paths from vertex(Dijkstra algorithm)\n"
                    "10. Topological Sorting\n"
                    "11. Build and display the spanning tree using DFS\n"
                    "12. Apply Kruskal's algorithm\n"
                    "0. Exit";


void userManagementMode() {

    int optionForGraphs = 0;
    optionForGraphs = getIntInput("Select graph you want to interact with:"
                        "1. Bit Vector Graph"
                       "2. Adjacency Matrix Graph"
                       "0. Exit", 0, 2);
    switch (optionForGraphs) {
        case 0: {
            cout << "Thanks for using system" << endl;
            return;
        }
        case 1: {
            userManagementModeBitGraph();
            break;
        }
        case 2: {
            userManagementModeAdjGraph();
            break;
        }
    }
}

void userManagementModeBitGraph() {
    BitVectorGraph bitVectorGraph(0);
    int mode = 0;
    do {
        mode = getIntInput(menu, 0, 12);
        switch (mode) {
            case 0: {
                cout << "Thanks for using system" << endl;
                return;
            }
            case 1: {
                int graphV = getIntInput("Enter number of vertices in the graph (from 0 to 1000)", 0, 1000);
                bitVectorGraph = BitVectorGraph(graphV);
                break;
            }
            case 2: {
                int graphV = getIntInput("Enter number of vertices in the graph (from 0 to 32)", 0, MAX_VERTICES);
                bitVectorGraph = BitVectorGraph(graphV);
                int graphE = getIntInput("Enter number of edges in the graph", 0, 1000);
                int graphMaxWeight = getIntInput("Enter max weight of the graph", 0, 1000);
                bitVectorGraph.createRandomGraph(graphE, graphMaxWeight);
                break;
            }
            case 3: {
                int u, v;
                cout << "Enter the source vertex (u): ";
                cin >> u;
                v = getIntInput("Enter the destination vertex(v) (from 0 to 32):", 0, MAX_VERTICES);

                bitVectorGraph.addEdge(u, v);
                break;
            }
            case 4: {
                if (bitVectorGraph.isDirectedBitVector()) {
                    cout << "The graph is directed." << endl;
                } else {
                    cout << "The graph is undirected." << endl;
                }
                break;
            }
            case 5: {
                bitVectorGraph.printGraph();
                break;
            }
            case 6: {
                cout << "BitVectorGraph is " << (isGraphConnectedBit(bitVectorGraph) ? "connected" : "disconnected") << endl;
                break;
            }
            case 7:
            case 8: {
                bool increasingWeights = (mode == 7); // Check if option 7 is selected
                vector<bool> visited(bitVectorGraph.V, false);
                if (increasingWeights) {
                    cout << "DFS Increasing Edge Weights:" << endl;
                    bitVectorGraph.dfsIncreasingEdgeWeights(0, visited);
                    cout << endl;
                }
                else {
                    bitVectorGraph.findConnectedComponents();
                    cout << endl;
                }
                break;
            }
            case 9: {
                bitVectorGraph.dijkstra(0);
                break;
            }
            case 10: {
                cout << "Topological Sorting:" << endl;
                bitVectorGraph.topologicalSort();
                break;
            }
            case 11: {
                // Build and display the spanning tree using DFS
                buildSpanningTreeDFSBit(bitVectorGraph);
                break;
            }
            case 12: {
                // Apply Kruskal's algorithm to find the minimum spanning tree
                kruskalBit(bitVectorGraph);
                break;
            }

        }
    } while (true);
}

void userManagementModeAdjGraph() {
    AdjMatrixGraph adjVectorGraph(0);
    int mode = 0;
    do {
        mode = getIntInput(menu, 0, 12);
        switch (mode) {
            case 0: {
                cout << "Thanks for using system" << endl;
                return;
            }
            case 1: {
                int graphV = getIntInput("Enter number of vertices in the graph (from 0 to 1000)", 0, 1000);
                adjVectorGraph = AdjMatrixGraph(graphV);
                break;
            }
            case 2: {
                int graphV = getIntInput("Enter number of vertices in the graph (from 0 to 32)", 0, MAX_VERTICES);
                adjVectorGraph = AdjMatrixGraph(graphV);
                int graphE = getIntInput("Enter number of edges in the graph", 0, 1000);
                int graphMaxWeight = getIntInput("Enter max weight of the graph", 0, 1000);
                adjVectorGraph.createRandomGraph(graphE, graphMaxWeight);
                break;
            }
            case 3: {
                int u, v;
                double weight;
                cout << "Enter the source vertex (u): ";
                cin >> u;
                v = getIntInput("Enter the destination vertex(v) (from 0 to 32):", 0, MAX_VERTICES);
                cout << "Enter the weight of the edge: ";
                cin >> weight;

                adjVectorGraph.addEdge(u, v, weight);
                break;
            }
            case 4: {
                if (adjVectorGraph.isDirectedAdjMatrix()) {
                    cout << "The graph is directed." << endl;
                } else {
                    cout << "The graph is undirected." << endl;
                }
                break;
            }
            case 5: {
                adjVectorGraph.printGraph();
                break;
            }
            case 6: {
                cout << "BitVectorGraph is " << (isGraphConnectedAdj(adjVectorGraph) ? "connected" : "disconnected") << endl;
                break;
            }
            case 7:
            case 8: {
                bool increasingWeights = (mode == 7); // Check if option 7 is selected
                vector<bool> visited(adjVectorGraph.V, false);
                if (increasingWeights) {
                    cout << "DFS Increasing Edge Weights:" << endl;
                    adjVectorGraph.dfsIncreasingEdgeWeights(0, visited);
                    cout << endl;
                }
                else {
                    adjVectorGraph.findConnectedComponents();
                    cout << endl;
                }
                break;
            }
            case 9: {
                adjVectorGraph.dijkstra(0);
                break;
            }
            case 10: {
                cout << "Topological Sorting:" << endl;
                adjVectorGraph.topologicalSort();
                break;
            }
            case 11: {
                // Build and display the spanning tree using DFS
                buildSpanningTreeDFSAdj(adjVectorGraph);
                break;
            }
            case 12: {
                // Apply Kruskal's algorithm to find the minimum spanning tree
                kruskalAdj(adjVectorGraph);
                break;
            }

        }
    } while (true);
}


void demonstrationMode() {
    cout << "\n**** DEMONSTRATION MODE *****\n";


    cout << "\n\n@@@@@ Creating a random graph with random vertices and edges(BitVectorGraph). @@@@@\n\n";
    loadingAnimation();
    loadingAnimation();
    int V = rand() % MAX_VERTICES; // Generate a random number for V
    int E = rand() % (V * (V - 1) / 2 + 1); // Generate a random number between 0 and maximum possible edges for V

    cout << "Randomly generated values:" << endl;
    cout << "Number of vertices (V): " << V << endl;
    cout << "Number of edges (E): " << E << endl;
    double weight = rand() % 100 / 10.0; // Generate a random
    BitVectorGraph bitVectorGraph(V);
    bitVectorGraph.createRandomGraph(E, weight);



    int u = rand() % MAX_VERTICES; // Generate a random vertex index
    int v = rand() % MAX_VERTICES; //  Generate a random vertex index
//    weight = rand() % 100 / 10.0;
    cout << "\n\n@@@@@  Adding an edge with random source vertex, destination vertex and the weight of the edge(BitVectorGraph).  @@@@@\n\n" << endl;
    loadingAnimation();
    loadingAnimation();
    cout << "Randomly generated values:" << endl;
    cout << "Source vertex (u): " << u << endl;
    cout << "Destination vertex (v): " << v << endl;
    cout << "Weight of the edge: " << weight << endl;
    bitVectorGraph.addEdge(u, v);




    cout << "\n\n@@@@@  Checking if the vector is Directed or not(BitVectorGraph). @@@@@\n\n" << endl;
    loadingAnimation();
    loadingAnimation();
    if (bitVectorGraph.isDirectedBitVector()) {
        cout << "The graph is directed." << endl;
    } else {
        cout << "The graph is undirected." << endl;
    }


    cout << "\n\n@@@@@  Printing the Graph(BitVectorGraph). @@@@@\n\n" << endl;
    loadingAnimation();
    bitVectorGraph.printGraph();


    cout << "\n\n@@@@@  Checking whether graph is connected or not(BitVectorGraph). @@@@@\n\n" << endl;
    loadingAnimation();
    loadingAnimation();
    cout << "BitVectorGraph is " << (isGraphConnectedBit(bitVectorGraph) ? "connected" : "disconnected") << endl;


    cout << "\n\n@@@@@  Applying DFS Algorithms for the graph(BitVectorGraph). @@@@@\n\n" << endl;
    loadingAnimation();
    loadingAnimation();
    vector<bool> visitedBit(bitVectorGraph.V, false);
    cout << "DFS Increasing Edge Weights:" << endl;
    bitVectorGraph.dfsIncreasingEdgeWeights(0, visitedBit);
    cout << endl;

    loadingAnimation();
    loadingAnimation();
    cout << "Looking for a connected components of the graph(BitVectorGraph)." << endl;
    bitVectorGraph.findConnectedComponents();
    cout << endl;


    cout << "\n\n@@@@@  Applying Dijkstra Algorithms for the graph(BitVectorGraph). @@@@@\n\n" << endl;
    loadingAnimation();
    loadingAnimation();
    bitVectorGraph.dijkstra(0);


    cout << "\n\n@@@@@ Topological Sorting for the graph(BitVectorGraph). @@@@@\n\n" << endl;
    loadingAnimation();
    loadingAnimation();
    bitVectorGraph.topologicalSort();


    cout << "\n\n@@@@@ Building Spanning Tree for the graph(BitVectorGraph). @@@@@\n\n" << endl;
    loadingAnimation();
    loadingAnimation();
    buildSpanningTreeDFSBit(bitVectorGraph);


    cout << "\n\n@@@@@ Finding the minimum spanning tree(Kruskal Algorithm) for the graph(BitVectorGraph). @@@@@\n\n" << endl;
    loadingAnimation();
    loadingAnimation();
    kruskalBit(bitVectorGraph);




    cout << "\n\n@@@@@ Creating a random graph with random vertices and edges(AdjMatrixGraph). @@@@@\n\n";
    loadingAnimation();
    loadingAnimation();

    cout << "Randomly generated values:" << endl;
    cout << "Number of vertices (V): " << V << endl;
    cout << "Number of edges (E): " << E << endl;
    AdjMatrixGraph adjMatrixGraph(V);
    adjMatrixGraph.createRandomGraph(E, weight);


    cout << "\n\n@@@@@  Adding an edge with random source vertex, destination vertex and the weight of the edge(AdjMatrixGraph).  @@@@@\n\n" << endl;
    loadingAnimation();
    loadingAnimation();
    cout << "Randomly generated values:" << endl;
    cout << "Source vertex (u): " << u << endl;
    cout << "Destination vertex (v): " << v << endl;
    cout << "Weight of the edge: " << weight << endl;
    adjMatrixGraph.addEdge(u, v, weight);



    cout << "\n\n@@@@@  Checking if the vector is Directed or not(AdjMatrixGraph). @@@@@\n\n" << endl;
    loadingAnimation();
    loadingAnimation();
    if (adjMatrixGraph.isDirectedAdjMatrix()) {
        cout << "The graph is directed." << endl;
    } else {
        cout << "The graph is undirected." << endl;
    }


    cout << "\n\n@@@@@  Printing the Graph(AdjMatrixGraph). @@@@@\n\n" << endl;
    loadingAnimation();
    adjMatrixGraph.printGraph();


    cout << "\n\n@@@@@  Checking whether graph is connected or not(AdjMatrixGraph). @@@@@\n\n" << endl;
    loadingAnimation();
    loadingAnimation();
    cout << "BitVectorGraph is " << (isGraphConnectedAdj(adjMatrixGraph) ? "connected" : "disconnected") << endl;


    cout << "\n\n@@@@@  Applying DFS Algorithms for the graph(AdjMatrixGraph). @@@@@\n\n" << endl;
    loadingAnimation();
    loadingAnimation();
    vector<bool> visitedAdj(adjMatrixGraph.V, false);
    cout << "DFS Increasing Edge Weights:" << endl;
    adjMatrixGraph.dfsIncreasingEdgeWeights(0, visitedAdj);
    cout << endl;

    loadingAnimation();
    loadingAnimation();

    cout << "Looking for a connected components of the graph(AdjMatrixGraph)." << endl;
    adjMatrixGraph.findConnectedComponents();
    cout << endl;


    cout << "\n\n@@@@@  Applying Dijkstra Algorithms for the graph(AdjMatrixGraph). @@@@@\n\n" << endl;
    loadingAnimation();
    loadingAnimation();
    adjMatrixGraph.dijkstra(0);


    cout << "\n\n@@@@@ Topological Sorting for the graph(AdjMatrixGraph). @@@@@\n\n" << endl;
    loadingAnimation();
    loadingAnimation();
    adjMatrixGraph.topologicalSort();


    cout << "\n\n@@@@@ Building Spanning Tree for the graph(AdjMatrixGraph). @@@@@\n\n" << endl;
    loadingAnimation();
    loadingAnimation();
    buildSpanningTreeDFSAdj(adjMatrixGraph);


    cout << "\n\n@@@@@ Finding the minimum spanning tree(Kruskal Algorithm) for the graph(AdjMatrixGraph). @@@@@\n\n" << endl;
    loadingAnimation();
    loadingAnimation();
    kruskalAdj(adjMatrixGraph);

    cout << "\n\n\n\n\n@@@@@ SUMMARY: @@@@@\n\n\n" << endl;
}


void benchmarkMode() {


}