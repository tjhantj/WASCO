#include "exp_func.h"
#include "functions.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <unordered_set>

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
// MEI method
// -------------------------------------------------------------------------
MeiMethod g_mei_method = MeiMethod::Mei;
long long g_ff_calls = 0;

bool parse_mei_method(const std::string& name, MeiMethod& m) {
    if (name == "mei") m = MeiMethod::Mei;
    else if (name == "lbfirst") m = MeiMethod::LbFirst;
    else if (name == "binary") m = MeiMethod::Binary;
    else if (name == "linear") m = MeiMethod::Linear;
    else if (name == "galloping") m = MeiMethod::Galloping;
    else return false;
    return true;
}

const char* mei_method_name(MeiMethod m) {
    switch (m) {
        case MeiMethod::Mei: return "mei";
        case MeiMethod::LbFirst: return "lbfirst";
        case MeiMethod::Binary: return "binary";
        case MeiMethod::Linear: return "linear";
        case MeiMethod::Galloping: return "galloping";
    }
    return "?";
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
        g_ff_calls++;
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

    if (counter) (*counter)++;
    auto follow = [&](int d) {
        g_ff_calls++;
        return FindFollowers(eu, ev, d, G, s, coreness);
    };

    // MEI from the endpoint support in H = C u F
    auto mei_from_support = [&](const std::unordered_set<int>& F) {
        auto support_in_H = [&](int x) {
            long long d = 0;
            for (auto& [nb, w] : G.adj[x]) {
                if (G.label[nb] || F.count(nb)) d += w;
            }
            return d;
        };
        bool u_in_H = G.label[eu] || F.count(eu);
        bool v_in_H = G.label[ev] || F.count(ev);
        if (!u_in_H || !v_in_H) return -1;
        long long need = std::max({0LL, (long long)s - support_in_H(eu),
                                   (long long)s - support_in_H(ev)});
        return (int)need;
    };

    // Smallest d in [lo, hi] with followers, given F(lo-1) is empty and F = F(hi)
    auto binary_search = [&](int lo, int hi, std::unordered_set<int>& F) {
        while (lo < hi) {
            int mid = lo + (hi - lo) / 2;
            auto Fm = follow(mid);
            if (Fm.empty()) {
                lo = mid + 1;
            } else {
                hi = mid;
                F = std::move(Fm);
            }
        }
        return hi;
    };

    std::unordered_set<int> F;
    int exact_delta = 0;
    switch (g_mei_method) {
        case MeiMethod::Mei:
            F = follow(probe_delta);
            if (F.empty()) return res;
            exact_delta = mei_from_support(F);
            break;
        case MeiMethod::LbFirst:
            F = follow(lower_bound);
            if (!F.empty()) {
                exact_delta = lower_bound;
                break;
            }
            if (probe_delta == lower_bound) return res;
            F = follow(probe_delta);
            if (F.empty()) return res;
            exact_delta = mei_from_support(F);
            break;
        case MeiMethod::Binary:
            F = follow(probe_delta);
            if (F.empty()) return res;
            exact_delta = binary_search(lower_bound, probe_delta, F);
            break;
        case MeiMethod::Linear:
            for (int d = lower_bound; d <= probe_delta; d++) {
                F = follow(d);
                if (!F.empty()) { exact_delta = d; break; }
            }
            if (F.empty()) return res;
            break;
        case MeiMethod::Galloping: {
            int prev = lower_bound - 1, d = lower_bound, gap = 1;
            while (true) {
                F = follow(d);
                if (!F.empty()) break;
                if (d == probe_delta) return res;
                prev = d;
                d = std::min(probe_delta, d + gap);
                gap *= 2;
            }
            exact_delta = binary_search(prev + 1, d, F);
            break;
        }
    }

    if (exact_delta < 1 || exact_delta < lower_bound || exact_delta > probe_delta) {
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
