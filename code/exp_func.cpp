#include "exp_func.h"
#include "functions.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

// -------------------------------------------------------------------------
// compute_delta_e
// -------------------------------------------------------------------------
int compute_delta_e(Graph& G, int s, int u, int v,
                    const std::string& delta_tactic,
                    const Coreness& coreness,
                    int orig_u, int orig_v) {
    int hu = (orig_u >= 0) ? orig_u : u;
    int hv = (orig_v >= 0) ? orig_v : v;

    if (delta_tactic == "compute") {
        return computeDelta(G, s, u, v, coreness);
    }
    if (delta_tactic == "mincost") {
        throw std::logic_error("compute_delta_e: mincost has no a-priori delta, use evaluate_candidate");
    }

    if (G.has_edge(hu, hv)) {
        return s - G.get_weight(hu, hv);
    } else {
        return s;
    }
}

// -------------------------------------------------------------------------
// is_better
// -------------------------------------------------------------------------
bool is_better(const BestEdge& candidate, const BestEdge& current) {
    if (!current.valid) return candidate.valid;
    if (!candidate.valid) return false;
    long long lhs = (long long)candidate.num_followers * current.delta;
    long long rhs = (long long)current.num_followers * candidate.delta;
    if (lhs != rhs) return lhs > rhs;
    return candidate.edge < current.edge;
}

// -------------------------------------------------------------------------
// exceeds
// -------------------------------------------------------------------------
static bool exceeds(double fr, double bound) {
    return fr > bound + 1e-9 * std::max(1.0, std::fabs(bound));
}

// -------------------------------------------------------------------------
// evaluate_candidate
// -------------------------------------------------------------------------
CandidateEval evaluate_candidate(Graph& G, int s, int budget_left, int eu, int ev,
                                 const std::string& delta_tactic, Coreness& coreness,
                                 int orig_u, int orig_v, long long* counter) {
    CandidateEval res;

    if (delta_tactic != "mincost") {
        int delta_e = compute_delta_e(G, s, eu, ev, delta_tactic, coreness, orig_u, orig_v);
        if (delta_e <= 0 || delta_e > budget_left) return res;
        auto followers = FindFollowers(eu, ev, delta_e, G, s, coreness);
        if (counter) (*counter)++;
        res.feasible = true;
        res.lower_bound = computeDelta(G, s, eu, ev, coreness);
        res.probe_delta = delta_e;
        res.delta = delta_e;
        res.gain = (int)followers.size();
        return res;
    }

    if (eu == ev || budget_left <= 0) return res;
    if (G.label[eu] && G.label[ev]) return res;

    int lower_bound = computeDelta(G, s, eu, ev, coreness);
    if (lower_bound > budget_left) return res;

    int w_before = G.get_weight(eu, ev);
    int probe_delta = std::min(s - w_before, budget_left);
    if (probe_delta < 1) {
        throw std::logic_error("evaluate_candidate: probe_delta < 1 for a non-core candidate");
    }

    auto F = FindFollowers(eu, ev, probe_delta, G, s, coreness);
    if (counter) (*counter)++;
    if (F.empty()) return res;

    auto support_in_H = [&](int x) {
        long long d = 0;
        for (auto& [nb, w] : G.adj[x]) {
            if (G.label[nb] || F.count(nb)) d += w;
        }
        return d;
    };
    long long deg_u = support_in_H(eu);
    long long deg_v = support_in_H(ev);
    long long need = std::max({0LL, (long long)s - deg_u, (long long)s - deg_v});
    int exact_delta = (int)need;

    bool u_in_H = G.label[eu] || F.count(eu);
    bool v_in_H = G.label[ev] || F.count(ev);
    if (!u_in_H || !v_in_H || exact_delta < 1 ||
        exact_delta < lower_bound || exact_delta > probe_delta) {
        throw std::logic_error("evaluate_candidate: inconsistent mincost result for (" +
                               std::to_string(eu) + "," + std::to_string(ev) + ")");
    }

    res.feasible = true;
    res.lower_bound = lower_bound;
    res.probe_delta = probe_delta;
    res.delta = exact_delta;
    res.gain = (int)F.size();
    return res;
}

static BestEdge make_best(int eu, int ev, const CandidateEval& ce, bool core_rep) {
    BestEdge be({std::min(eu, ev), std::max(eu, ev)}, ce.delta,
                (double)ce.gain / ce.delta, ce.gain);
    be.lower_bound = ce.lower_bound;
    be.core_rep = core_rep;
    return be;
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
    return {non_s_core, s_cand};
}

// -------------------------------------------------------------------------
// make_candidate_nodes
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
// make_candidate_nodes_v2
// -------------------------------------------------------------------------
std::vector<int> make_candidate_nodes_v2(Graph& G, int s, int budget_left,
                                          const Coreness& coreness,
                                          std::unordered_map<int,double>& upperbound,
                                          bool T2_upperbound,
                                          const std::string& delta_tactic) {
    std::vector<int> candidate_nodes;
    for (int u = 0; u < G.n; u++) {
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
    std::sort(candidate_nodes.begin(), candidate_nodes.end(),
              [&](int a, int b) { return upperbound[a] > upperbound[b]; });
    return candidate_nodes;
}

// -------------------------------------------------------------------------
// make_candidate_edges
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

    int c = (int)candidate_nodes.size();
    for (int i = 0; i < c; i++) {
        int u = candidate_nodes[i];
        double ub_u = U_single(u, upperbound);

        if (exceeds(most_FR, ub_u * 2.0)) break;

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
            if (exceeds(most_FR, ub_u + ub_v)) break;

            if (u == v && s_cand < 0) continue;

            double udb = U_double(u, v, upperbound, coreness, G, s);
            if (exceeds(most_FR, udb)) continue;

            int eu = u, ev = v;
            if (u == v) ev = s_cand;

            CandidateEval ce = evaluate_candidate(G, s, b - spent, eu, ev, delta_tactic,
                                                  coreness, u, v, counter);
            if (!ce.feasible) continue;

            BestEdge cand = make_best(eu, ev, ce, u == v);
            if (is_better(cand, best)) {
                best = cand;
                most_FR = best.FR;
            }
        }
    }

    return best;
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

            if (u == v && s_cand < 0) continue;

            int eu = u, ev = v;
            if (u == v) ev = s_cand;

            CandidateEval ce = evaluate_candidate(G, s, b - spent, eu, ev, delta_tactic,
                                                  coreness, u, v, counter);
            if (!ce.feasible) continue;

            BestEdge cand = make_best(eu, ev, ce, u == v);
            if (is_better(cand, best)) best = cand;
        }
    }

    return best;
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

    for (auto& [u, v] : candidate_edges) {
        if (T3_reuse) {
            auto itu = comp_of.find(u);
            auto itv = comp_of.find(v);
            if (itu != comp_of.end() && itv != comp_of.end() &&
                itu->second == itv->second) {
                continue;
            }
        }

        CandidateEval ce = evaluate_candidate(G, s, b - spent, u, v, delta_tactic,
                                              coreness, -1, -1, counter);
        if (!ce.feasible) continue;

        BestEdge cand = make_best(u, v, ce, false);
        if (is_better(cand, best)) best = cand;
    }

    return best;
}
