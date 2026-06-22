#include "exp_func.h"
#include "functions.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

// -------------------------------------------------------------------------
// Helper: compute delta_e based on tactic
// -------------------------------------------------------------------------
int compute_delta_e(Graph& G, int s, int u, int v,
                    const std::string& delta_tactic,
                    const Coreness& coreness,
                    int orig_u, int orig_v) {
    // For naive/percentage has_edge check, use original (u,v) before self-edge transform
    // If not provided (-1), fall back to the transformed (u,v)
    int hu = (orig_u >= 0) ? orig_u : u;
    int hv = (orig_v >= 0) ? orig_v : v;

    if (delta_tactic == "compute") {
        return computeDelta(G, s, u, v, coreness);
    }

    // Check if delta_tactic is a digit string (e.g. "50")
    bool is_digit = !delta_tactic.empty();
    for (char c : delta_tactic) {
        if (!std::isdigit(c)) { is_digit = false; break; }
    }

    if (is_digit) {
        double ratio = std::stoi(delta_tactic) / 100.0;
        int mee = computeDelta(G, s, u, v, coreness);
        int naive;
        if (G.has_edge(hu, hv)) {
            naive = s - G.get_weight(hu, hv);
        } else {
            naive = s;
        }
        int diff = naive - mee;
        return mee + (int)std::ceil(diff * ratio);
    }

    // naive tactic
    if (G.has_edge(hu, hv)) {
        return s - G.get_weight(hu, hv);
    } else {
        return s;
    }
}

// -------------------------------------------------------------------------
// is_better: check if candidate replaces current best
// -------------------------------------------------------------------------
bool is_better(const BestEdge& candidate, const BestEdge& current) {
    if (!current.valid) return candidate.valid;
    if (!candidate.valid) return false;
    if (candidate.FR > current.FR) return true;
    if (candidate.FR == current.FR && candidate.edge < current.edge) return true;
    return false;
}

// -------------------------------------------------------------------------
// self_edge_pruning
// -------------------------------------------------------------------------
std::pair<std::vector<int>, int> self_edge_pruning(const Graph& G) {
    std::vector<int> non_s_core;
    int s_cand = -1;
    for (int n = 0; n < G.n; n++) {
        if (!G.label[n]) {
            non_s_core.push_back(n);
        } else {
            if (s_cand == -1 || n < s_cand) {
                s_cand = n;
            }
        }
    }
    if (s_cand == -1) {
        throw std::runtime_error("No node in s-core. Change s value");
    }
    return {non_s_core, s_cand};
}

// -------------------------------------------------------------------------
// make_candidate_nodes (T1 tactic)
// -------------------------------------------------------------------------
std::vector<int> make_candidate_nodes(Graph& G, const std::vector<int>& nodes, int s,
                                       int budget_left, const Coreness& coreness,
                                       std::unordered_map<int,double>& upperbound,
                                       bool T2_upperbound,
                                       const std::string& delta_tactic) {
    std::vector<int> candidate_nodes;
    for (int u : nodes) {
        if (!G.label[u] && s - coreness[u].first <= budget_left) {
            candidate_nodes.push_back(u);
            if (T2_upperbound) {
                upperbound[u] = Upperbound(G, u, coreness, s, delta_tactic);
            }
        }
    }
    if (T2_upperbound) {
        std::sort(candidate_nodes.begin(), candidate_nodes.end(),
                  [&](int a, int b) { return upperbound[a] > upperbound[b]; });
    }
    return candidate_nodes;
}

// -------------------------------------------------------------------------
// make_candidate_nodes_v2 (T2 but not T1)
// -------------------------------------------------------------------------
std::vector<int> make_candidate_nodes_v2(Graph& G, int s, int budget_left,
                                          const Coreness& coreness,
                                          std::unordered_map<int,double>& upperbound,
                                          bool T2_upperbound,
                                          const std::string& delta_tactic) {
    std::vector<int> candidate_nodes;
    for (int u = 0; u < G.n; u++) {
        // coreness.get(u, (s,0))[0]
        int cu = G.label[u] ? s : coreness[u].first;
        if (s - cu > budget_left) continue;

        if (!G.label[u]) {
            upperbound[u] = Upperbound(G, u, coreness, s, delta_tactic);
        } else {
            if (T2_upperbound) {
                upperbound[u] = 0.0;
            }
        }
        candidate_nodes.push_back(u);
    }
    // Always sort by -upperbound
    std::sort(candidate_nodes.begin(), candidate_nodes.end(),
              [&](int a, int b) { return upperbound[a] > upperbound[b]; });
    return candidate_nodes;
}

// -------------------------------------------------------------------------
// make_candidate_edges (no T1, no T2)
// -------------------------------------------------------------------------
std::vector<std::pair<int,int>> make_candidate_edges(const Graph& G, int s, int budget_left,
                                                      const Coreness& coreness) {
    std::vector<int> non_s_core, s_core;
    for (int u = 0; u < G.n; u++) {
        int cu = G.label[u] ? s : coreness[u].first;
        if (s - cu > budget_left) continue;
        if (!G.label[u]) {
            non_s_core.push_back(u);
        } else {
            s_core.push_back(u);
        }
    }

    std::vector<std::pair<int,int>> edges;
    int nn = (int)non_s_core.size();
    for (int i = 0; i < nn; i++) {
        for (int j = i + 1; j < nn; j++) {
            edges.push_back({non_s_core[i], non_s_core[j]});
        }
    }
    for (int u : non_s_core) {
        for (int v : s_core) {
            edges.push_back({u, v});
        }
    }
    return edges;
}

// -------------------------------------------------------------------------
// iteration_nodes_upperbound
// -------------------------------------------------------------------------
BestEdge iteration_nodes_upperbound(Graph& G, const std::vector<int>& candidate_nodes,
                                     int s, int b, int spent,
                                     Coreness& coreness,
                                     std::unordered_map<int,double>& upperbound,
                                     int s_cand,
                                     bool T1_self_edge, bool T3_reuse,
                                     const std::string& delta_tactic,
                                     const std::unordered_map<int,int>& comp_of,
                                     const BestEdge& initial_best,
                                     long long* counter) {
    BestEdge best = initial_best;
    double most_FR = best.valid ? best.FR : 0.0;
    // We need to track best_edge as optional-like
    bool has_best = best.valid;
    std::pair<int,int> best_edge = best.edge;
    int best_delta = best.delta;
    int best_followers = best.num_followers;

    int c = (int)candidate_nodes.size();
    for (int i = 0; i < c; i++) {
        int u = candidate_nodes[i];
        double ub_u = U_single(u, upperbound);

        if (most_FR > ub_u * 2.0) break;

        int j_start = T1_self_edge ? i : i + 1;
        for (int j = j_start; j < c; j++) {
            int v = candidate_nodes[j];

            if (T3_reuse) {
                auto itu = comp_of.find(u);
                auto itv = comp_of.find(v);
                if (itu != comp_of.end() && itv != comp_of.end() &&
                    itu->second == itv->second) {
                    continue;
                }
            }

            double ub_v = U_single(v, upperbound);
            if (most_FR > ub_u + ub_v) break;

            double udb = U_double(u, v, upperbound, coreness, G, s);
            if (most_FR > udb) continue;

            // Compute actual edge
            int eu = u, ev = v;
            if (u == v) ev = s_cand;

            int delta_e = compute_delta_e(G, s, eu, ev, delta_tactic, coreness, u, v);

            if (delta_e > 0 && spent + delta_e <= b) {
                auto followers = FindFollowers(eu, ev, delta_e, G, s, coreness);
                if (counter) (*counter)++;
                double FR = (double)followers.size() / delta_e;

                // Canonical edge (eu <= ev)
                std::pair<int,int> e_canon = {std::min(eu, ev), std::max(eu, ev)};

                if (FR > most_FR || (FR == most_FR && (!has_best || e_canon < best_edge))) {
                    best_edge = e_canon;
                    best_delta = delta_e;
                    most_FR = FR;
                    best_followers = (int)followers.size();
                    has_best = true;
                }
            }
        }
    }

    if (!has_best) return BestEdge{};
    return BestEdge{best_edge, best_delta, most_FR, best_followers};
}

// -------------------------------------------------------------------------
// iteration_nodes_no_upperbound
// -------------------------------------------------------------------------
BestEdge iteration_nodes_no_upperbound(Graph& G, const std::vector<int>& candidate_nodes,
                                        int s, int b, int spent,
                                        Coreness& coreness,
                                        int s_cand,
                                        bool T3_reuse,
                                        const std::string& delta_tactic,
                                        const std::unordered_map<int,int>& comp_of,
                                        const BestEdge& initial_best,
                                        long long* counter) {
    BestEdge best = initial_best;
    double most_FR = best.valid ? best.FR : 0.0;
    bool has_best = best.valid;
    std::pair<int,int> best_edge = best.edge;
    int best_delta = best.delta;
    int best_followers = best.num_followers;

    int c = (int)candidate_nodes.size();
    for (int i = 0; i < c; i++) {
        int u = candidate_nodes[i];
        for (int j = i; j < c; j++) {
            int v = candidate_nodes[j];

            if (T3_reuse) {
                auto itu = comp_of.find(u);
                auto itv = comp_of.find(v);
                if (itu != comp_of.end() && itv != comp_of.end() &&
                    itu->second == itv->second) {
                    continue;
                }
            }

            int eu = u, ev = v;
            if (u == v) ev = s_cand;

            int delta_e = compute_delta_e(G, s, eu, ev, delta_tactic, coreness, u, v);

            if (delta_e > 0 && spent + delta_e <= b) {
                auto followers = FindFollowers(eu, ev, delta_e, G, s, coreness);
                if (counter) (*counter)++;
                double FR = (double)followers.size() / delta_e;

                std::pair<int,int> e_canon = {std::min(eu, ev), std::max(eu, ev)};

                if (FR > most_FR || (FR == most_FR && (!has_best || e_canon < best_edge))) {
                    best_edge = e_canon;
                    best_delta = delta_e;
                    most_FR = FR;
                    best_followers = (int)followers.size();
                    has_best = true;
                }
            }
        }
    }

    if (!has_best) return BestEdge{};
    return BestEdge{best_edge, best_delta, most_FR, best_followers};
}

// -------------------------------------------------------------------------
// iteration_edges_no_upperbound
// -------------------------------------------------------------------------
BestEdge iteration_edges_no_upperbound(Graph& G,
                                        const std::vector<std::pair<int,int>>& candidate_edges,
                                        int s, int b, int spent,
                                        Coreness& coreness,
                                        bool T3_reuse,
                                        const std::string& delta_tactic,
                                        const std::unordered_map<int,int>& comp_of,
                                        const BestEdge& initial_best,
                                        long long* counter) {
    BestEdge best = initial_best;
    double most_FR = best.valid ? best.FR : 0.0;
    bool has_best = best.valid;
    std::pair<int,int> best_edge = best.edge;
    int best_delta = best.delta;
    int best_followers = best.num_followers;

    for (auto& [u, v] : candidate_edges) {
        if (T3_reuse) {
            auto itu = comp_of.find(u);
            auto itv = comp_of.find(v);
            if (itu != comp_of.end() && itv != comp_of.end() &&
                itu->second == itv->second) {
                continue;
            }
        }

        int delta_e = compute_delta_e(G, s, u, v, delta_tactic, coreness);

        if (delta_e > 0 && spent + delta_e <= b) {
            auto followers = FindFollowers(u, v, delta_e, G, s, coreness);
            if (counter) (*counter)++;
            double FR = (double)followers.size() / delta_e;

            std::pair<int,int> e_canon = {std::min(u, v), std::max(u, v)};

            if (FR > most_FR || (FR == most_FR && (!has_best || e_canon < best_edge))) {
                best_edge = e_canon;
                best_delta = delta_e;
                most_FR = FR;
                best_followers = (int)followers.size();
                has_best = true;
            }
        }
    }

    if (!has_best) return BestEdge{};
    return BestEdge{best_edge, best_delta, most_FR, best_followers};
}
