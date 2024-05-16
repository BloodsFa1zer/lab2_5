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
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <queue>
#include <set>
#include <stack>
#include <thread>
#include "benchmarkOS_X.h"
#include "benchmark.h"

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
    int src, dest;
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
    bitset<MAX_VERTICES> adjList[MAX_VERTICES];
    int V; // Number of vertices

    BitVectorGraph(int V) : V(V) {}


    // Function to add an edge between vertices u and v
    void addEdge(int u, int v) {
        // Check if the vertices are valid
        if (u < 0 || u >= V || v < 0 || v >= V) {
            cout << "Error: Invalid vertices." << endl;
            return;
        }

        // Add the edge between u and v
        adjList[u].set(v);
        adjList[v].set(u);
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
    void createRandomGraph(int E, double maxWeight, bool directed) {
        srand(time(NULL));

        int edgesAdded = 0;
        while (edgesAdded < E) {
            int u = rand() % V;
            int v = rand() % V;
            if (u != v) { // Ensure the edge is not a self-loop

                // Check if the edge already exists
                if (!directed && !adjList[u].test(v) && !adjList[v].test(u)) {
                    // For undirected graph, add edge in both directions
                    addEdge(u, v);
                    addEdge(v, u);
                    ++edgesAdded;
                } else if (directed && !adjList[u].test(v)) {
                    // For directed graph, add edge only if it doesn't already exist
                    addEdge(u, v);
                    ++edgesAdded;
                }
            }
        }
    }



// Depth-First Search (DFS) function to explore connected components
    void DFS(int v, vector<bool>& visited) {
        visited[v] = true;
        cout << v << " ";

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
                DFS(v, visited);
                cout << endl;
            }
        }
    }

    // Function to perform DFS traversal on the graph with neighbors traversed in order of increasing edge weights
    void dfsIncreasingEdgeWeights(int v, vector<bool>& visited) {
        visited[v] = true;
        cout << v << " ";

        set<int> adjVertices;

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
        vector<double> dist(V, INF);
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
    visited[v] = true;

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
            parent[v] = u;
            dfsSpanningTreeBit(v, graph, visited, parent);
        }
    }
}

// Function to construct a spanning tree based on depth-first search
void buildSpanningTreeDFSBit(const BitVectorGraph& graph) {
    int V = graph.V;
    vector<bool> visited(V, false);
    vector<int> parent(V, -1);
    int totalWeight = 0;

    for (int u = 0; u < V; ++u) {
        if (!visited[u]) {
            dfsSpanningTreeBit(u, graph, visited, parent);
        }
    }

    // Calculate total weight of the spanning tree
    for (int v = 0; v < V; ++v) {
        if (parent[v] != -1) {
            // Since the bit vector graph does not have weights, we cannot calculate the total weight
            // we'll simply count the number of edges in the spanning tree
            totalWeight += 1;
        }
    }

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

    double totalWeight = 0;

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
    int V;
    vector<vector<int>> adjMatrix;


    AdjMatrixGraph(int V) : V(V) {
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
            return;
        }

        // Check for self-loop
        if (u == v) {
            cout << "Error: Self-loops are not allowed." << endl;
            return;
        }

        // Add the edge between u and v with weight
        adjMatrix[u][v] = weight;
        adjMatrix[v][u] = weight;
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

    // Function to create a random graph with specified number of edges
    void createRandomGraph(int E, double maxWeight, bool directed) {
        srand(time(NULL));

        int edgesAdded = 0;
        while (edgesAdded < E) {
            int u = rand() % V;
            int v = rand() % V;
            if (u != v) { // Ensure the edge is not a self-loop

                // Check if the edge already exists
                if (!directed && adjMatrix[u][v] == 0 && adjMatrix[v][u] == 0) {
                    // For undirected graph, add edge in both directions
                    addEdge(u, v, maxWeight);
                    addEdge(v, u, maxWeight);
                    ++edgesAdded;
                } else if (directed && adjMatrix[u][v] == 0) {
                    // For directed graph, add edge only if it doesn't already exist
                    addEdge(u, v, maxWeight);
                    ++edgesAdded;
                }
            }
        }
    }



    // Depth-First Search (DFS) function to explore connected components
    void DFS(int v, vector<bool>& visited) {
        visited[v] = true;
        cout << v << " ";

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
        cout << endl;
        cout << "Connected Components:" << endl;
        for (int v = 0; v < V; ++v) {
            if (!visited[v]) {
                DFS(v, visited); // Call DFS for each unvisited vertex
                cout << endl;
            }
        }
    }


    // Function to perform DFS traversal on the graph with neighbors traversed in order of increasing edge weights
    void dfsIncreasingEdgeWeights(int v, vector<bool>& visited) {
        visited[v] = true;
        cout << v << " ";

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

        for (int i = 0; i < V; ++i) {
            if (adjMatrix[v][i] && !visited[i]) {
                topologicalSortUtil(i, visited, stack);
            }
        }

        stack.push(v);
    }

    // Function to perform topological sorting
    void topologicalSort() {
        vector<bool> visited(V, false);
        stack<int> stack;

        for (int i = 0; i < V; ++i) {
            if (!visited[i]) {
                topologicalSortUtil(i, visited, stack);
            }
        }

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

    for (int i = 0; i < graph.V; ++i) {
        if (graph.adjMatrix[v][i] && !visited[i]) {
            dfsAdj(graph, i, visited); // Recursively visit unvisited adjacent vertices
        }
    }
}

// Function to check if the graph is connected using Depth-First Search (DFS)
bool isGraphConnectedAdj(const AdjMatrixGraph& graph) {
    vector<bool> visited(graph.V, false);
    dfsAdj(graph, 0, visited);

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
            parent[v] = u;
            dfsSpanningTreeAdj(v, adjMatrix, visited, parent);
        }
    }
}

// Function to construct a spanning tree based on depth-first search
void buildSpanningTreeDFSAdj(AdjMatrixGraph& graph) {
    int V = graph.V;
    vector<bool> visited(V, false);
    vector<int> parent(V, -1);
    int totalWeight = 0;

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
        for (int v = u + 1; v < V; ++v) {
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

    return 0;
}

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
                bitVectorGraph.createRandomGraph(graphE, graphMaxWeight, false);
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
                adjVectorGraph.createRandomGraph(graphE, graphMaxWeight, false);
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
    bitVectorGraph.createRandomGraph(E, weight, false);



    int u = rand() % MAX_VERTICES; // Generate a random vertex index
    int v = rand() % MAX_VERTICES; //  Generate a random vertex index
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
    adjMatrixGraph.createRandomGraph(E, weight, false);


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
}


void benchmarkMode() {
    cout << "\nYou are now in the benchmark mode, time of each operation will be written afterwards\n";

    cout << "Bit Vector Graph, directed" << endl;
    int V = rand() % MAX_VERTICES;
    int E = rand() % (V * (V - 1) / 2 + 1);

    double weight = rand() % 100 / 10.0;
    auto startTimeBitDirected = high_resolution_clock::now();
    BitVectorGraph bitVectorGraphDirected(V);
    auto startTimeBitDirectedCreating = high_resolution_clock::now();
     BenchmarkMax memoryBenchmarkDirBit("Creating a random graph(Directed Bit Graph)");
    bitVectorGraphDirected.createRandomGraph(E, weight, true);
    auto endTimeBitDirectedCreating = high_resolution_clock::now();


    auto startTimeBitDirectedPrinting = high_resolution_clock::now();
     BenchmarkMax printBenchmarkDirBit("Printing graph(Directed Bit Graph)");
    bitVectorGraphDirected.printGraph();
    auto endTimeBitDirectedPrinting = high_resolution_clock::now();

    auto startTimeBitDirectedIsConnected = high_resolution_clock::now();
     BenchmarkMax connectionBenchmarkDirBit("Checking whether graph connected or not(Directed Bit Graph)");
    isGraphConnectedBit(bitVectorGraphDirected);
    auto endTimeBitDirectedIsConnected = high_resolution_clock::now();


    auto startTimeBitDirectedIncreasingEdgeWeights = high_resolution_clock::now();
     BenchmarkMax dfsIncreasingEdgeWeightsBenchmarkDirBit("Applying DFS Algorithms for the graph - IncreasingEdgeWeights(Directed Bit Graph)");
    vector<bool> visitedBitDirected(bitVectorGraphDirected.V, false);
    bitVectorGraphDirected.dfsIncreasingEdgeWeights(0, visitedBitDirected);
    auto endTimeBitDirectedIncreasingEdgeWeights = high_resolution_clock::now();


    auto startTimeBitDirectedFindingConnectedComponents = high_resolution_clock::now();
     BenchmarkMax findConnectedComponentsBenchmarkDirBit("Applying DFS Algorithms for the graph - findConnectedComponents(Directed Bit Graph)");
    bitVectorGraphDirected.findConnectedComponents();
    auto endTimeBitDirectedFindingConnectedComponents = high_resolution_clock::now();



    auto startTimeBitDirectedDijkstra = high_resolution_clock::now();
     BenchmarkMax dijkstraBenchmarkDirBit(" Applying Dijkstra Algorithms for the graph(Directed Bit Graph)");
    bitVectorGraphDirected.dijkstra(0);
    auto endTimeBitDirectedDijkstra = high_resolution_clock::now();


    auto startTimeBitDirectedTopologicalSort = high_resolution_clock::now();
     BenchmarkMax topologicalSortBenchmarkDirBit("Applying topological sorting for the graph(Directed Bit Graph)");
    bitVectorGraphDirected.topologicalSort();
    auto endTimeBitDirectedTopologicalSort= high_resolution_clock::now();


    auto startTimeBitDirectedBuildSpanningTreeDFSBit = high_resolution_clock::now();
     BenchmarkMax buildSpanningTreeDFSBitBenchmarkDirBit("Building Spanning Tree for the graph(Directed Bit Graph)");
    buildSpanningTreeDFSBit(bitVectorGraphDirected);
    auto endTimeBitDirectedBuildSpanningTreeDFSBit = high_resolution_clock::now();


    auto startTimeBitDirectedKruskalBit = high_resolution_clock::now();
     BenchmarkMax kruskalBitBenchmarkDirBit("Finding the minimum spanning tree(Kruskal Algorithm) for the graph(Directed Bit Graph).");
    kruskalBit(bitVectorGraphDirected);
    auto endTimeBitDirectedKruskalBit = high_resolution_clock::now();

    auto endTimeBitDirected = high_resolution_clock::now();





    cout << "Bit Vector Graph, Undirected" << endl;
    auto startTimeBitUndirected = high_resolution_clock::now();

    BitVectorGraph bitVectorGraphUndirected(V);
    auto startTimeBitUndirectedCreating = high_resolution_clock::now();
     BenchmarkMax memoryBenchmarkUndirectedBit("Creating a random graph(Undirected Bit Graph).");
    bitVectorGraphUndirected.createRandomGraph(E, weight, false);
    auto endTimeBitUndirectedCreating = high_resolution_clock::now();


    auto startTimeBitUndirectedPrinting = high_resolution_clock::now();
     BenchmarkMax printBenchmarkUndirectedBit("Printing graph(Undirected Bit Graph).");
    bitVectorGraphUndirected.printGraph();
    auto endTimeBitUndirectedPrinting = high_resolution_clock::now();

    auto startTimeBitUndirectedIsConnected = high_resolution_clock::now();
     BenchmarkMax connectionBenchmarkUndirectedBit("Checking whether graph connected or not(Undirected Bit Graph).");
    isGraphConnectedBit(bitVectorGraphUndirected);
    auto endTimeBitUndirectedIsConnected = high_resolution_clock::now();


    auto startTimeBitUndirectedIncreasingEdgeWeights = high_resolution_clock::now();
     BenchmarkMax dfsIncreasingEdgeWeightsBenchmarkUndirectedBit("Applying DFS Algorithms for the graph - IncreasingEdgeWeights(Undirected Bit Graph).");
    vector<bool> visitedBit(bitVectorGraphUndirected.V, false);
    bitVectorGraphUndirected.dfsIncreasingEdgeWeights(0, visitedBit);
    auto endTimeBitUndirectedIncreasingEdgeWeights = high_resolution_clock::now();


    auto startTimeBitUndirectedFindingConnectedComponents = high_resolution_clock::now();
     BenchmarkMax findConnectedComponentsBenchmarkUndirectedBit("Applying DFS Algorithms for the graph - findConnectedComponents(Undirected Bit Graph).");
    bitVectorGraphUndirected.findConnectedComponents();
    auto endTimeBitUndirectedFindingConnectedComponents = high_resolution_clock::now();



    auto startTimeBitUndirectedDijkstra = high_resolution_clock::now();
     BenchmarkMax dijkstraBenchmarkUndirectedBit(" Applying Dijkstra Algorithms for the graph(Undirected Bit Graph).");
    bitVectorGraphUndirected.dijkstra(0);
    auto endTimeBitUndirectedDijkstra = high_resolution_clock::now();


    auto startTimeBitUndirectedTopologicalSort = high_resolution_clock::now();
     BenchmarkMax topologicalSortBenchmarkUndirectedBit("Applying topological sorting for the graph(Undirected Bit Graph).");
    bitVectorGraphUndirected.topologicalSort();
    auto endTimeBitUndirectedTopologicalSort= high_resolution_clock::now();


    auto startTimeBitUndirectedBuildSpanningTreeDFSBit = high_resolution_clock::now();
     BenchmarkMax buildSpanningTreeDFSBitBenchmarkUndirectedBit("Building Spanning Tree for the graph(Undirected Bit Graph)");
    buildSpanningTreeDFSBit(bitVectorGraphUndirected);
    auto endTimeBitUndirectedBuildSpanningTreeDFSBit = high_resolution_clock::now();


    auto startTimeBitUndirectedKruskalBit = high_resolution_clock::now();
     BenchmarkMax kruskalBitBenchmarkUndirectedBit("Finding the minimum spanning tree(Kruskal Algorithm) for the graph(Undirected Bit Graph).");
    kruskalBit(bitVectorGraphUndirected);
    auto endTimeBitUndirectedKruskalBit = high_resolution_clock::now();


    auto endTimeBitUndirected = high_resolution_clock::now();




    cout << "Adjusted Matrix Graph, directed" << endl;


    auto startTimeAdjDirected = high_resolution_clock::now();
    AdjMatrixGraph adjMatrixGraphDirected(V);
    auto startTimeAdjDirectedCreating = high_resolution_clock::now();
     BenchmarkMax memoryBenchmarkDirAdj("Creating a random graph(Directed Adj Graph)");
    adjMatrixGraphDirected.createRandomGraph(E, weight, true);
    auto endTimeAdjDirectedCreating = high_resolution_clock::now();


    auto startTimeAdjDirectedPrinting = high_resolution_clock::now();
     BenchmarkMax printBenchmarkDirAdj("Printing graph(Directed Adj Graph)");
    adjMatrixGraphDirected.printGraph();
    auto endTimeAdjDirectedPrinting = high_resolution_clock::now();

    auto startTimeAdjDirectedIsConnected = high_resolution_clock::now();
     BenchmarkMax connectionBenchmarkDirAdj("Checking whether graph connected or not(Directed Adj Graph)");
    isGraphConnectedAdj(adjMatrixGraphDirected);
    auto endTimeAdjDirectedIsConnected = high_resolution_clock::now();


    auto startTimeAdjDirectedIncreasingEdgeWeights = high_resolution_clock::now();
     BenchmarkMax dfsIncreasingEdgeWeightsBenchmarkDirAdj("Applying DFS Algorithms for the graph - IncreasingEdgeWeights(Directed Adj Graph)");
    vector<bool> visitedAdjDirected(adjMatrixGraphDirected.V, false);
    adjMatrixGraphDirected.dfsIncreasingEdgeWeights(0, visitedAdjDirected);
    auto endTimeAdjDirectedIncreasingEdgeWeights = high_resolution_clock::now();


    auto startTimeAdjDirectedFindingConnectedComponents = high_resolution_clock::now();
     BenchmarkMax findConnectedComponentsBenchmarkDirAdj("Applying DFS Algorithms for the graph - findConnectedComponents(Directed Adj Graph)");
    adjMatrixGraphDirected.findConnectedComponents();
    auto endTimeAdjDirectedFindingConnectedComponents = high_resolution_clock::now();



    auto startTimeAdjDirectedDijkstra = high_resolution_clock::now();
     BenchmarkMax dijkstraBenchmarkDirAdj("Applying Dijkstra Algorithms for the graph(Directed Adj Graph)");
    adjMatrixGraphDirected.dijkstra(0);
    auto endTimeAdjDirectedDijkstra = high_resolution_clock::now();


    auto startTimeAdjDirectedTopologicalSort = high_resolution_clock::now();
     BenchmarkMax topologicalSortBenchmarkDirAdj("Applying topological sorting for the graph(Directed Adj Graph)");
    adjMatrixGraphDirected.topologicalSort();
    auto endTimeAdjDirectedTopologicalSort= high_resolution_clock::now();


    auto startTimeAdjDirectedBuildSpanningTreeDFS = high_resolution_clock::now();
     BenchmarkMax buildSpanningTreeDFSBenchmarkDirAdj("Building Spanning Tree for the graph(Directed Adj Graph)");
    buildSpanningTreeDFSAdj(adjMatrixGraphDirected);
    auto endTimeAdjDirectedBuildSpanningTreeDFS = high_resolution_clock::now();


    auto startTimeAdjDirectedKruskal = high_resolution_clock::now();
     BenchmarkMax kruskalBitBenchmarkDirAdj("Finding the minimum spanning tree(Kruskal Algorithm) for the graph(Directed Adj Graph).");
    kruskalAdj(adjMatrixGraphDirected);
    auto endTimeAdjDirectedKruskal = high_resolution_clock::now();

    auto endTimeAdjDirected = high_resolution_clock::now();




    cout << "Adjusted Matrix Graph, undirected" << endl;


    auto startTimeAdjUndirected = high_resolution_clock::now();
    AdjMatrixGraph adjMatrixGraphUndirected(V);
    auto startTimeAdjUndirectedCreating = high_resolution_clock::now();
     BenchmarkMax memoryBenchmarkUndirectedAdj("Creating a random graph(Undirected Adj Graph)");
    adjMatrixGraphUndirected.createRandomGraph(E, weight, false);
    auto endTimeAdjUndirectedCreating = high_resolution_clock::now();


    auto startTimeAdjUndirectedPrinting = high_resolution_clock::now();
     BenchmarkMax printBenchmarkUndirectedAdj("Printing graph(Undirected Adj Graph)");
    adjMatrixGraphUndirected.printGraph();
    auto endTimeAdjUndirectedPrinting = high_resolution_clock::now();

    auto startTimeAdjUndirectedIsConnected = high_resolution_clock::now();
     BenchmarkMax connectionBenchmarkUndirectedAdj("Checking whether graph connected or not(Undirected Adj Graph)");
    isGraphConnectedAdj(adjMatrixGraphUndirected);
    auto endTimeAdjUndirectedIsConnected = high_resolution_clock::now();


    auto startTimeAdjUndirectedIncreasingEdgeWeights = high_resolution_clock::now();
     BenchmarkMax dfsIncreasingEdgeWeightsBenchmarkUndirectedAdj("Applying DFS Algorithms for the graph - IncreasingEdgeWeights(Undirected Adj Graph)");
    vector<bool> visitedAdjUndirected(adjMatrixGraphUndirected.V, false);
    adjMatrixGraphUndirected.dfsIncreasingEdgeWeights(0, visitedAdjUndirected);
    auto endTimeAdjUndirectedIncreasingEdgeWeights = high_resolution_clock::now();


    auto startTimeAdjUndirectedFindingConnectedComponents = high_resolution_clock::now();
     BenchmarkMax findConnectedComponentsBenchmarkUndirectedAdj("Applying DFS Algorithms for the graph - findConnectedComponents(Undirected Adj Graph)");
    adjMatrixGraphUndirected.findConnectedComponents();
    auto endTimeAdjUndirectedFindingConnectedComponents = high_resolution_clock::now();



    auto startTimeAdjUndirectedDijkstra = high_resolution_clock::now();
     BenchmarkMax dijkstraBenchmarkUndirectedAdj("Applying Dijkstra Algorithms for the graph(Undirected Adj Graph)");
    adjMatrixGraphUndirected.dijkstra(0);
    auto endTimeAdjUndirectedDijkstra = high_resolution_clock::now();



    auto startTimeAdjUndirectedTopologicalSort = high_resolution_clock::now();
     BenchmarkMax topologicalSortBenchmarkUndirectedAdj("Applying topological sorting for the graph(Undirected Adj Graph)");
    adjMatrixGraphUndirected.topologicalSort();
    auto endTimeAdjUndirectedTopologicalSort= high_resolution_clock::now();


    auto startTimeAdjUndirectedBuildSpanningTreeDFS = high_resolution_clock::now();
     BenchmarkMax buildSpanningTreeDFSBenchmarkUndirectedAdj("Building Spanning Tree for the graph(Undirected Adj Graph)");
    buildSpanningTreeDFSAdj(adjMatrixGraphUndirected);
    auto endTimeAdjUndirectedBuildSpanningTreeDFS = high_resolution_clock::now();


    auto startTimeAdjUndirectedKruskal = high_resolution_clock::now();
     BenchmarkMax kruskalAdjBenchmarkUndirectedAdj("Finding the minimum spanning tree(Kruskal Algorithm) for the graph(Undirected Adj Graph).");
    kruskalAdj(adjMatrixGraphUndirected);
    auto endTimeAdjUndirectedKruskal = high_resolution_clock::now();

    auto endTimeAdjUndirected = high_resolution_clock::now();




    auto durationBitDirected = duration_cast<duration<double, milli>>(endTimeBitDirected - startTimeBitDirected).count();
    auto durationBitDirectedCreating = duration_cast<duration<double, nano>>(endTimeBitDirectedCreating - startTimeBitDirectedCreating).count();
    auto durationBitDirectedPrinting = duration_cast<duration<double, nano>>(endTimeBitDirectedPrinting - startTimeBitDirectedPrinting).count();
    auto durationBitDirectedConnection = duration_cast<duration<double, nano>>(endTimeBitDirectedIsConnected - startTimeBitDirectedIsConnected).count();
    auto durationBitDirectedDfsIncreasingEdgeWeights = duration_cast<duration<double, nano>>(endTimeBitDirectedBuildSpanningTreeDFSBit - startTimeBitDirectedBuildSpanningTreeDFSBit).count();
    auto durationBitDirectedFindConnectedComponents = duration_cast<duration<double, nano>>(endTimeBitDirectedFindingConnectedComponents - startTimeBitDirectedFindingConnectedComponents).count();
    auto durationBitDirectedDijkstra = duration_cast<duration<double, nano>>(endTimeBitDirectedDijkstra - startTimeBitDirectedDijkstra).count();
    auto durationBitDirectedTopologicalSort = duration_cast<duration<double, nano>>(endTimeBitDirectedTopologicalSort - startTimeBitDirectedTopologicalSort).count();
    auto durationBitDirectedBuildSpanningTreeDFS = duration_cast<duration<double, nano>>(endTimeBitDirectedBuildSpanningTreeDFSBit - startTimeBitDirectedBuildSpanningTreeDFSBit).count();
    auto durationBitDirectedKruskal = duration_cast<duration<double, nano>>(endTimeBitDirectedKruskalBit - startTimeBitDirectedKruskalBit).count();



    auto durationBitUndirected = duration_cast<duration<double, milli>>(endTimeBitUndirected - startTimeBitUndirected).count();
    auto durationBitUndirectedCreating = duration_cast<duration<double, nano>>(endTimeBitUndirectedCreating - startTimeBitUndirectedCreating).count();
    auto durationBitUndirectedPrinting = duration_cast<duration<double, nano>>(endTimeBitUndirectedPrinting - startTimeBitUndirectedPrinting).count();
    auto durationBitUndirectedConnection = duration_cast<duration<double, nano>>(endTimeBitUndirectedIsConnected - startTimeBitUndirectedIsConnected).count();
    auto durationBitUndirectedDfsIncreasingEdgeWeights = duration_cast<duration<double, nano>>(endTimeBitUndirectedBuildSpanningTreeDFSBit - startTimeBitUndirectedBuildSpanningTreeDFSBit).count();
    auto durationBitUndirectedFindConnectedComponents = duration_cast<duration<double, nano>>(endTimeBitUndirectedFindingConnectedComponents - startTimeBitUndirectedFindingConnectedComponents).count();
    auto durationBitUndirectedDijkstra = duration_cast<duration<double, nano>>(endTimeBitUndirectedDijkstra - startTimeBitUndirectedDijkstra).count();
    auto durationBitUndirectedTopologicalSort = duration_cast<duration<double, nano>>(endTimeBitUndirectedTopologicalSort - startTimeBitUndirectedTopologicalSort).count();
    auto durationBitUndirectedBuildSpanningTreeDFS = duration_cast<duration<double, nano>>(endTimeBitUndirectedBuildSpanningTreeDFSBit - startTimeBitUndirectedBuildSpanningTreeDFSBit).count();
    auto durationBitUndirectedKruskal = duration_cast<duration<double, nano>>(endTimeBitUndirectedKruskalBit - startTimeBitUndirectedKruskalBit).count();



    auto durationAdjDirected = duration_cast<duration<double, milli>>(endTimeAdjDirected - startTimeAdjDirected).count();
    auto durationAdjDirectedCreating = duration_cast<duration<double, nano>>(endTimeAdjDirectedCreating - startTimeAdjDirectedCreating).count();
    auto durationAdjDirectedPrinting = duration_cast<duration<double, nano>>(endTimeAdjDirectedPrinting - startTimeAdjDirectedPrinting).count();
    auto durationAdjDirectedConnection = duration_cast<duration<double, nano>>(endTimeAdjDirectedIsConnected - startTimeAdjDirectedIsConnected).count();
    auto durationAdjDirectedDfsIncreasingEdgeWeights = duration_cast<duration<double, nano>>(endTimeAdjDirectedBuildSpanningTreeDFS - startTimeAdjDirectedBuildSpanningTreeDFS).count();
    auto durationAdjDirectedFindConnectedComponents = duration_cast<duration<double, nano>>(endTimeAdjDirectedFindingConnectedComponents - startTimeAdjDirectedFindingConnectedComponents).count();
    auto durationAdjDirectedDijkstra = duration_cast<duration<double, nano>>(endTimeAdjDirectedDijkstra - startTimeAdjDirectedDijkstra).count();
    auto durationAdjDirectedTopologicalSort = duration_cast<duration<double, nano>>(endTimeAdjDirectedTopologicalSort - startTimeAdjDirectedTopologicalSort).count();
    auto durationAdjDirectedBuildSpanningTreeDFS = duration_cast<duration<double, nano>>(endTimeAdjDirectedBuildSpanningTreeDFS - startTimeAdjDirectedBuildSpanningTreeDFS).count();
    auto durationAdjDirectedKruskal = duration_cast<duration<double, nano>>(endTimeAdjDirectedKruskal - startTimeAdjDirectedKruskal).count();



    auto durationAdjUndirected = duration_cast<duration<double, milli>>(endTimeAdjUndirected - startTimeAdjUndirected).count();
    auto durationAdjUndirectedCreating = duration_cast<duration<double, nano>>(endTimeAdjUndirectedCreating - startTimeBitUndirectedCreating).count();
    auto durationAdjUndirectedPrinting = duration_cast<duration<double, nano>>(endTimeAdjUndirectedPrinting - startTimeBitUndirectedPrinting).count();
    auto durationAdjUndirectedConnection = duration_cast<duration<double, nano>>(endTimeAdjUndirectedIsConnected - startTimeBitUndirectedIsConnected).count();
    auto durationAdjUndirectedDfsIncreasingEdgeWeights = duration_cast<duration<double, nano>>(endTimeAdjUndirectedBuildSpanningTreeDFS - startTimeAdjUndirectedBuildSpanningTreeDFS).count();
    auto durationAdjUndirectedFindConnectedComponents = duration_cast<duration<double, nano>>(endTimeAdjUndirectedFindingConnectedComponents - startTimeAdjUndirectedFindingConnectedComponents).count();
    auto durationAdjUndirectedDijkstra = duration_cast<duration<double, nano>>(endTimeAdjUndirectedDijkstra - startTimeAdjUndirectedDijkstra).count();
    auto durationAdjUndirectedTopologicalSort = duration_cast<duration<double, nano>>(endTimeAdjUndirectedTopologicalSort - startTimeAdjUndirectedTopologicalSort).count();
    auto durationAdjUndirectedBuildSpanningTreeDFS = duration_cast<duration<double, nano>>(endTimeAdjUndirectedBuildSpanningTreeDFS - startTimeAdjUndirectedBuildSpanningTreeDFS).count();
    auto durationAdjUndirectedKruskal = duration_cast<duration<double, nano>>(endTimeAdjUndirectedKruskal - startTimeAdjUndirectedKruskal).count();


    cout << "\n\n\n\n---------------------------------------------------" << endl;
    cout << "\nTime All Bit Graph Directed: " << durationBitDirected << " milliseconds\n" <<
         "\nTime Creating Bit Graph Directed: " << durationBitDirectedCreating << " nanoseconds\n" <<
         "\nTime Printing Bit Graph Directed: " << durationBitDirectedPrinting << " nanoseconds\n" <<
         "\nTime Bit Directed Connection Directed:" << durationBitDirectedConnection << " nanoseconds\n" <<
         "\nTime Dfs Increasing Edge Weights Bit Graph Directed: " << durationBitDirectedDfsIncreasingEdgeWeights << " nanoseconds\n" <<
         "\nTime Find Connected Components Bit Graph Directed: " << durationBitDirectedFindConnectedComponents << " nanoseconds\n" <<
         "\nTime Dijkstra Bit Graph Directed: " << durationBitDirectedDijkstra << " nanoseconds\n" <<
         "\nTime Topological Sort Bit Graph Directed: " << durationBitDirectedTopologicalSort << " nanoseconds\n" <<
         "\nTime Build Spanning Tree DFS Bit Graph Directed: " << durationBitDirectedBuildSpanningTreeDFS << " nanoseconds\n" <<
         "\nTime Kruskal Bit Graph Directed: " << durationBitDirectedKruskal << " nanoseconds\n";


    cout << "\n\n\n\n---------------------------------------------------" << endl;
    cout << "\nTime All Bit Graph Undirected: " << durationBitUndirected << " milliseconds\n" <<
         "\nTime Creating Bit Graph Undirected: " << durationBitUndirectedCreating << " nanoseconds\n" <<
         "\nTime Printing Bit Graph Undirected: " << durationBitUndirectedPrinting << " nanoseconds\n" <<
         "\nTime Bit Directed Connection Undirected:" << durationBitUndirectedConnection << " nanoseconds\n" <<
         "\nTime Dfs Increasing Edge Weights Bit Graph Undirected: " << durationBitUndirectedDfsIncreasingEdgeWeights << " nanoseconds\n" <<
         "\nTime Find Connected Components Bit Graph Undirected: " << durationBitUndirectedFindConnectedComponents << " nanoseconds\n" <<
         "\nTime Dijkstra Bit Graph Undirected: " << durationBitUndirectedDijkstra << " nanoseconds\n" <<
         "\nTime Topological Sort Bit Graph Undirected: " << durationBitUndirectedTopologicalSort << " nanoseconds\n" <<
         "\nTime Build Spanning Tree DFS Bit Graph Undirected: " << durationBitUndirectedBuildSpanningTreeDFS << " nanoseconds\n" <<
         "\nTime Kruskal Bit Graph Undirected: " << durationBitUndirectedKruskal << " nanoseconds\n";



    cout << "\n\n\n\n---------------------------------------------------" << endl;
    cout << "\nTime All Adj Graph Directed: " << durationAdjDirected << " milliseconds\n" <<
         "\nTime Creating Adj Graph Directed: " << durationAdjDirectedCreating << " nanoseconds\n" <<
         "\nTime Printing Adj Graph Directed: " << durationAdjDirectedPrinting << " nanoseconds\n" <<
         "\nTime Adj Directed Connection Directed:" << durationAdjDirectedConnection << " nanoseconds\n" <<
         "\nTime Dfs Increasing Edge Weights Adj Graph Directed: " << durationAdjDirectedDfsIncreasingEdgeWeights << " nanoseconds\n" <<
         "\nTime Find Connected Components Adj Graph Directed: " << durationAdjDirectedFindConnectedComponents << " nanoseconds\n" <<
         "\nTime Dijkstra Adj Graph Directed: " << durationAdjDirectedDijkstra << " nanoseconds\n" <<
         "\nTime Topological Sort Adj Graph Directed: " << durationAdjDirectedTopologicalSort << " nanoseconds\n" <<
         "\nTime Build Spanning Tree DFS Adj Graph Directed: " << durationAdjDirectedBuildSpanningTreeDFS << " nanoseconds\n" <<
         "\nTime Kruskal Adj Graph Directed: " << durationAdjDirectedKruskal << " nanoseconds\n";



    cout << "\n\n\n\n---------------------------------------------------" << endl;
    cout << "\nTime All Adj Graph Undirected: " << durationAdjUndirected << " milliseconds\n" <<
         "\nTime Creating Adj Graph Undirected: " << durationAdjUndirectedCreating << " nanoseconds\n" <<
         "\nTime Printing Adj Graph Undirected: " << durationAdjUndirectedPrinting << " nanoseconds\n" <<
         "\nTime Adj Directed Connection Undirected:" << durationAdjUndirectedConnection << " nanoseconds\n" <<
         "\nTime Dfs Increasing Edge Weights Adj Graph Undirected: " << durationAdjUndirectedDfsIncreasingEdgeWeights << " nanoseconds\n" <<
         "\nTime Find Connected Components Adj Graph Undirected: " << durationAdjUndirectedFindConnectedComponents << " nanoseconds\n" <<
         "\nTime Dijkstra Adj Graph Undirected: " << durationAdjUndirectedDijkstra << " nanoseconds\n" <<
         "\nTime Topological Sort Adj Graph Undirected: " << durationAdjUndirectedTopologicalSort << " nanoseconds\n" <<
         "\nTime Build Spanning Tree DFS Adj Graph Undirected: " << durationAdjUndirectedBuildSpanningTreeDFS << " nanoseconds\n" <<
         "\nTime Kruskal Adj Graph Undirected: " << durationAdjUndirectedKruskal << " nanoseconds\n";


    cout << "\n\n    BENCHMARK IS ENDED UP WITH THE ABOVE-MENTIONED RESULTS    \n";
}