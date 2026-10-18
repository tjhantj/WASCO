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
    std::unordered_set<int> nodes_set(nodes.begin(), nodes.end());
    for (int node : nodes) {
        G.label[node] = true;
    }
    int s_core_num = (int)nodes.size();

    std::unordered_map<int,int> weight_sum;
    weight_sum.reserve(nodes.size());
    for (int node : nodes) {
        int ws = 0;
        for (auto& [nb, w] : G.adj[node]) {
            ws += w;
        }
        weight_sum[node] = ws;
    }

    using PII = std::pair<int,int>;
    std::priority_queue<PII, std::vector<PII>, std::greater<PII>> heap;
    for (int node : nodes) {
        heap.push({weight_sum[node], node});
    }

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

                if (!G.label[node]) continue;

                coreness[node] = {current_core, layer};

                for (auto& [neighbor, ew] : G.adj[node]) {
                    if (!nodes_set.count(neighbor)) continue;
                    if (!G.label[neighbor]) continue;
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
    bool edge_added = !G.has_edge(u, v);
    if (edge_added) {
        G.add_edge(u, v, delta_e);
    } else {
        G.increase_weight(u, v, delta_e);
    }
    auto rollback = [&]() {
        if (edge_added) {
            G.remove_edge(u, v);
        } else {
            G.increase_weight(u, v, -delta_e);
        }
    };

    enum : char { IN_PQ = 1, IN_F = 2, DEAD = 3 };
    std::unordered_map<int,char> state;
    std::unordered_map<int,long long> sigma_plus;

    using Elem = std::tuple<std::pair<int,int>, int, int>;
    std::priority_queue<Elem, std::vector<Elem>, std::greater<Elem>> PQ;
    int index_PQ = 0;

    auto push_pq = [&](int nd) {
        state[nd] = IN_PQ;
        PQ.push({coreness[nd], index_PQ++, nd});
    };

    for (int nd : {u, v}) {
        if (!G.label[nd] && !state.count(nd)) {
            push_pq(nd);
        }
    }

    while (!PQ.empty()) {
        auto [cor_x, idx_x, x] = PQ.top();
        PQ.pop();

        long long sp = 0;
        for (auto& [nb, w] : G.adj[x]) {
            if (G.label[nb]) {
                sp += w;
                continue;
            }
            auto it = state.find(nb);
            if (it == state.end()) {
                if (coreness[nb] > coreness[x]) sp += w;
            } else if (it->second != DEAD) {
                sp += w;
            }
        }

        if (sp >= s) {
            sigma_plus[x] = sp;
            state[x] = IN_F;
            for (auto& [y, w] : G.adj[x]) {
                if (!G.label[y] && !state.count(y) && coreness[y] > coreness[x]) {
                    push_pq(y);
                }
            }
            continue;
        }

        state[x] = DEAD;
        if (x == u || x == v) {
            rollback();
            return {};
        }

        std::deque<int> Q;
        Q.push_back(x);
        while (!Q.empty()) {
            int y = Q.front();
            Q.pop_front();
            for (auto& [z, w] : G.adj[y]) {
                auto it = state.find(z);
                if (it == state.end() || it->second != IN_F) continue;
                sigma_plus[z] -= w;
                if (sigma_plus[z] < s) {
                    it->second = DEAD;
                    if (z == u || z == v) {
                        rollback();
                        return {};
                    }
                    Q.push_back(z);
                }
            }
        }
    }

    rollback();

    std::unordered_set<int> F;
    for (auto& [nd, st] : state) {
        if (st == IN_F) F.insert(nd);
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
    if (delta_tactic == "compute" || delta_tactic == "mincost") {
        denominator = (double)(s - coreness[u].first);
    } else {
        int max_weight = 0;
        for (auto& [nb, w] : G.adj[u]) {
            max_weight = std::max(max_weight, w);
        }
        denominator = (double)(s - max_weight);
    }

    if (denominator <= 0) {
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
