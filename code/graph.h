#pragma once
#include <vector>
#include <unordered_map>
#include <string>
#include <utility>

// Adjacency list: adj[u][v] = weight of edge (u,v)
struct Graph {
    int n;  // number of nodes (0-indexed)
    std::vector<std::unordered_map<int,int>> adj;
    std::vector<bool> label;  // label[v] = true means v is in s-core (active)
    // Tracks the order in which nodes first appeared when edges were added.
    // Used to replicate Python's G.nodes() insertion-order CC building.
    std::vector<int> node_insertion_order;

    Graph() : n(0) {}
    explicit Graph(int num_nodes) : n(num_nodes), adj(num_nodes), label(num_nodes, false) {
        // Default: nodes in 0..n-1 order
        node_insertion_order.resize(num_nodes);
        for (int i = 0; i < num_nodes; i++) node_insertion_order[i] = i;
    }

    bool has_edge(int u, int v) const {
        if (u < 0 || u >= n) return false;
        return adj[u].count(v) > 0;
    }

    int get_weight(int u, int v) const {
        auto it = adj[u].find(v);
        if (it == adj[u].end()) return 0;
        return it->second;
    }

    void add_edge(int u, int v, int w) {
        adj[u][v] = w;
        adj[v][u] = w;
    }

    void increase_weight(int u, int v, int dw) {
        adj[u][v] += dw;
        adj[v][u] += dw;
    }

    void remove_edge(int u, int v) {
        adj[u].erase(v);
        adj[v].erase(u);
    }

    // Deep copy
    Graph copy() const {
        Graph g(n);
        g.adj = adj;
        g.label = label;
        g.node_insertion_order = node_insertion_order;
        return g;
    }
};

// Load graph from weighted edgelist file (space/tab separated: u v w per line)
// Applies 0-indexing relabel and ceil on weights.
Graph load_graph(const std::string& path);
