#include "functions.h"
#include <queue>
#include <algorithm>
#include <deque>
#include <limits>
#include <cmath>

// -------------------------------------------------------------------------
// calculate_s_core
// -------------------------------------------------------------------------
int calculate_s_core(Graph& G, const std::vector<int>& nodes, int s, Coreness& coreness) {
    // Initialize all nodes in the set to label=true
    // Build a set for O(1) membership test — critical for subset (T3 reuse) calls
    std::unordered_set<int> nodes_set(nodes.begin(), nodes.end());
    for (int node : nodes) {
        G.label[node] = true;
    }
    int s_core_num = (int)nodes.size();

    // weight_sum[node] = sum of ALL neighbor weights (includes s-core neighbors
    // outside the subset — they contribute to the initial weight but are never
    // removed, so their weight stays accounted for implicitly)
    std::unordered_map<int,int> weight_sum;
    weight_sum.reserve(nodes.size());
    for (int node : nodes) {
        int ws = 0;
        for (auto& [nb, w] : G.adj[node]) {
            ws += w;
        }
        weight_sum[node] = ws;
    }

    // Min-heap: (weight, node) — lazy deletion via label check
    using PII = std::pair<int,int>;
    std::priority_queue<PII, std::vector<PII>, std::greater<PII>> heap;
    for (int node : nodes) {
        heap.push({weight_sum[node], node});
    }

    // Core decomposition with layer tracking
    while (!heap.empty()) {
        auto [current_core, top_node] = heap.top();
        if (current_core >= s) break;

        int layer = 0;
        while (!heap.empty() && heap.top().first <= current_core) {
            layer++;
            std::unordered_map<int,int> temp;

            while (!heap.empty() && heap.top().first <= current_core) {
                auto [weight, node] = heap.top();
                heap.pop();

                // lazy deletion: skip if already processed
                if (!G.label[node]) continue;

                coreness[node] = {current_core, layer};

                // Update neighbors — ONLY those in the nodes_set.
                // Python: "if neighbor not in nodes: continue"
                // Using label==false as "already removed" is NOT sufficient here:
                // s-core nodes outside the subset also have label==true.
                for (auto& [neighbor, ew] : G.adj[node]) {
                    if (!nodes_set.count(neighbor)) continue;  // not in our subset
                    if (!G.label[neighbor]) continue;           // already removed
                    weight_sum[neighbor] -= ew;
                    temp[neighbor] = weight_sum[neighbor];
                }

                G.label[node] = false;
                s_core_num--;
            }

            for (auto& [nd, w] : temp) {
                heap.push({w, nd});
            }
        }
    }

    return s_core_num;
}

// -------------------------------------------------------------------------
// computeDelta
// -------------------------------------------------------------------------
int computeDelta(Graph& G, int s, int u, int v, const Coreness& coreness) {
    int c_u = G.label[u] ? s : coreness[u].first;
    int c_v = G.label[v] ? s : coreness[v].first;
    return s - std::min(c_u, c_v);
}

// -------------------------------------------------------------------------
// FindFollowers
// -------------------------------------------------------------------------
std::unordered_set<int> FindFollowers(int u, int v, int delta_e, Graph& G, int s, Coreness& coreness) {
    std::unordered_set<int> F;

    // Add or increase edge
    bool edge_added = false;
    if (G.has_edge(u, v)) {
        G.increase_weight(u, v, delta_e);
        edge_added = false;
    } else {
        G.add_edge(u, v, delta_e);
        edge_added = true;
    }

    // Build priority queue with non-s-core endpoints
    // PQ element: (coreness, index, node) - min-heap by coreness then index
    using Elem = std::tuple<std::pair<int,int>, int, int>;
    std::priority_queue<Elem, std::vector<Elem>, std::greater<Elem>> PQ;
    // Use a counter (not a set) to track PQ membership.
    // Python's heap may hold duplicate entries for the same node; the "in PQ" check
    // in sigma_plus looks at ALL heap entries. So we count how many copies each node
    // has in the heap: the node is "in PQ" as long as count > 0.
    std::unordered_map<int,int> pq_count;
    int index_PQ = 0;

    auto push_pq = [&](int nd) {
        PQ.push({coreness[nd], index_PQ++, nd});
        pq_count[nd]++;
    };

    if (u == v) {
        if (!G.label[u]) {
            push_pq(u);
        }
    } else {
        for (int nd : {u, v}) {
            if (!G.label[nd]) {
                push_pq(nd);
            }
        }
    }

    std::unordered_map<int,int> sigma_plus;

    while (!PQ.empty()) {
        auto [cor_x, idx_x, x] = PQ.top();
        PQ.pop();
        // Decrement count; remove from map when count reaches 0
        if (--pq_count[x] == 0) pq_count.erase(x);

        // σ⁺(x) = sum of weights to neighbors that are:
        // label==true (s-core) OR coreness > coreness[x] OR in F OR in PQ
        int sp = 0;
        for (auto& [nb, w] : G.adj[x]) {
            bool qualifies = false;
            if (G.label[nb]) {
                qualifies = true;
            } else if (coreness[nb] > coreness[x]) {
                qualifies = true;
            } else if (F.count(nb)) {
                qualifies = true;
            } else if (pq_count.count(nb)) {
                qualifies = true;
            }
            if (qualifies) sp += w;
        }
        sigma_plus[x] = sp;

        if (sp >= s) {
            F.insert(x);
            // Push neighbors that are non-s-core, not in F, coreness > coreness[x]
            for (auto& [y, w] : G.adj[x]) {
                if (!G.label[y] && !F.count(y) && coreness[y] > coreness[x]) {
                    push_pq(y);
                }
            }
        } else {
            // BFS cascade removal.
            // IMPORTANT: the neighbor propagation runs for EVERY node popped from Q,
            // not just those in F. This matches Python: the "for z in neighbors(y)"
            // loop is outside the "if y in F" block.
            std::deque<int> Q;
            Q.push_back(x);
            bool early_exit = false;
            while (!Q.empty() && !early_exit) {
                int y = Q.front();
                Q.pop_front();
                if (F.count(y)) {
                    F.erase(y);
                    if (y == u || y == v) {
                        early_exit = true;
                        break;
                    }
                }
                // Always propagate to F-neighbors (even if y was not in F)
                for (auto& [z, w] : G.adj[y]) {
                    if (F.count(z)) {
                        sigma_plus[z] -= w;
                        if (sigma_plus[z] < s) {
                            Q.push_back(z);
                        }
                    }
                }
            }
            if (early_exit) {
                if (edge_added) {
                    G.remove_edge(u, v);
                } else {
                    G.increase_weight(u, v, -delta_e);
                }
                return {};
            }
        }
    }

    // Roll back edge change
    if (edge_added) {
        G.remove_edge(u, v);
    } else {
        G.increase_weight(u, v, -delta_e);
    }

    return F;
}

// -------------------------------------------------------------------------
// Upperbound
// -------------------------------------------------------------------------
double Upperbound(Graph& G, int u, const Coreness& coreness, int s, const std::string& delta_tactic) {
    if (G.label[u]) return 0.0;

    int n = G.n;
    std::vector<bool> visited(n, false);
    std::deque<int> Q;
    int count = 1;
    Q.push_back(u);
    visited[u] = true;

    while (!Q.empty()) {
        int v = Q.front();
        Q.pop_front();
        for (auto& [w, ew] : G.adj[v]) {
            if (!G.label[w] && coreness[v] < coreness[w]) {
                if (!visited[w]) {
                    count++;
                    Q.push_back(w);
                    visited[w] = true;
                }
            }
        }
    }

    double denominator;
    if (delta_tactic == "compute") {
        denominator = (double)(s - coreness[u].first);
    } else {
        int max_weight = 0;
        for (auto& [nb, w] : G.adj[u]) {
            max_weight = std::max(max_weight, w);
        }
        denominator = (double)(s - max_weight);
    }

    if (denominator <= 0) {
        // Avoid division by zero; return large value
        return (double)count;
    }
    return (double)count / denominator;
}

// -------------------------------------------------------------------------
// U_single
// -------------------------------------------------------------------------
double U_single(int u, const std::unordered_map<int,double>& upperbound) {
    auto it = upperbound.find(u);
    if (it == upperbound.end()) return 0.0;
    return it->second;
}

// -------------------------------------------------------------------------
// U_double
// -------------------------------------------------------------------------
double U_double(int u, int v, const std::unordered_map<int,double>& upperbound,
                const Coreness& coreness, const Graph& G, int s) {
    if (G.label[v]) return U_single(u, upperbound);
    if (G.label[u]) return U_single(v, upperbound);

    if ((G.has_edge(u, v) && coreness[u] < coreness[v]) || u == v) {
        return U_single(u, upperbound);
    } else if (G.has_edge(u, v) && coreness[u] > coreness[v]) {
        return U_single(v, upperbound);
    } else {
        return U_single(u, upperbound) + U_single(v, upperbound);
    }
}
