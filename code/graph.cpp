#include "graph.h"
#include <fstream>
#include <sstream>
#include <cmath>
#include <map>
#include <set>
#include <algorithm>
#include <stdexcept>

Graph load_graph(const std::string& path) {
    std::ifstream fin(path);
    if (!fin.is_open()) {
        throw std::runtime_error("Cannot open file: " + path);
    }

    // Read all edges; nodes may not be contiguous or 0-indexed
    std::vector<std::tuple<int,int,double>> raw_edges;
    std::set<int> node_set;
    std::string line;
    while (std::getline(fin, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        int u, v;
        double w;
        if (!(ss >> u >> v >> w)) continue;
        node_set.insert(u);
        node_set.insert(v);
        raw_edges.emplace_back(u, v, w);
    }
    fin.close();

    // Build relabel map: sorted original ids -> 0-indexed
    std::vector<int> sorted_nodes(node_set.begin(), node_set.end());
    std::sort(sorted_nodes.begin(), sorted_nodes.end());
    std::map<int,int> relabel;
    for (int i = 0; i < (int)sorted_nodes.size(); i++) {
        relabel[sorted_nodes[i]] = i;
    }

    int n = (int)sorted_nodes.size();
    Graph G(n);

    // Track node insertion order (first appearance order, matching Python's G.nodes() order)
    std::vector<bool> seen(n, false);
    std::vector<int> insertion_order;
    insertion_order.reserve(n);

    for (auto& [ou, ov, ow] : raw_edges) {
        int u = relabel[ou];
        int v = relabel[ov];
        int w = (int)std::ceil(ow);
        if (u == v) continue;  // skip self loops in edge list

        // Track first appearance
        if (!seen[u]) { seen[u] = true; insertion_order.push_back(u); }
        if (!seen[v]) { seen[v] = true; insertion_order.push_back(v); }

        if (!G.has_edge(u, v)) {
            G.add_edge(u, v, w);
        } else {
            G.adj[u][v] = w;
            G.adj[v][u] = w;
        }
    }
    G.node_insertion_order = insertion_order;

    return G;
}
