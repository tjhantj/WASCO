#include "experiment.h"
#include "functions.h"
#include "exp_func.h"
#include <algorithm>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <stdexcept>

// -------------------------------------------------------------------------
// Connected components via BFS (for the non-s-core subgraph or full graph)
// -------------------------------------------------------------------------
static std::vector<std::unordered_set<int>> connected_components_subset(
    const Graph& G, const std::unordered_set<int>& node_set) {
    std::unordered_set<int> visited;
    std::vector<std::unordered_set<int>> components;

    // Iterate in G.node_insertion_order order (matches Python's G.nodes() order).
    // This ensures CC CID assignment matches Python's enumerate(nx.connected_components()).
    for (int start : G.node_insertion_order) {
        if (!node_set.count(start)) continue;
        if (visited.count(start)) continue;
        std::unordered_set<int> component;
        std::vector<int> stack = {start};
        visited.insert(start);
        while (!stack.empty()) {
            int cur = stack.back(); stack.pop_back();
            component.insert(cur);
            for (auto& [nb, w] : G.adj[cur]) {
                if (node_set.count(nb) && !visited.count(nb)) {
                    visited.insert(nb);
                    stack.push_back(nb);
                }
            }
        }
        components.push_back(std::move(component));
    }
    return components;
}

static std::vector<std::unordered_set<int>> connected_components_all(const Graph& G) {
    std::vector<bool> visited(G.n, false);
    std::vector<std::unordered_set<int>> components;
    for (int start = 0; start < G.n; start++) {
        if (visited[start]) continue;
        std::unordered_set<int> component;
        std::vector<int> stack = {start};
        visited[start] = true;
        while (!stack.empty()) {
            int cur = stack.back(); stack.pop_back();
            component.insert(cur);
            for (auto& [nb, w] : G.adj[cur]) {
                if (!visited[nb]) {
                    visited[nb] = true;
                    stack.push_back(nb);
                }
            }
        }
        components.push_back(std::move(component));
    }
    return components;
}

// -------------------------------------------------------------------------
// build_initial_caches
// -------------------------------------------------------------------------
CCCache build_initial_caches(Graph& G, int s, int b, int spent,
                               Coreness& coreness,
                               std::unordered_map<int,double>& upperbound,
                               bool T1_self_edge, bool T2_upperbound,
                               int& s_cand,
                               const std::string& delta_tactic,
                               long long* counter) {
    CCCache cache;
    s_cand = -1;

    if (T1_self_edge) {
        std::unordered_set<int> non_s_core_nodes;
        for (int v = 0; v < G.n; v++) {
            if (!G.label[v]) {
                non_s_core_nodes.insert(v);
            } else {
                if (s_cand == -1 || v < s_cand) s_cand = v;
            }
        }
        if (s_cand == -1) {
            throw std::runtime_error("No node in s-core. Change s value");
        }

        auto ccs = connected_components_subset(G, non_s_core_nodes);
        int cid = 0;
        for (auto& nodes : ccs) {
            cache.nodes_in[cid] = nodes;
            for (int v : nodes) {
                cache.comp_of[v] = cid;
            }
            cid++;
        }
    } else {
        auto ccs = connected_components_all(G);
        int cid = 0;
        for (auto& nodes : ccs) {
            cache.nodes_in[cid] = nodes;
            for (int v : nodes) {
                cache.comp_of[v] = cid;
            }
            cid++;
        }
    }

    // Build intra-best for each CC
    for (auto& [cid, nodes] : cache.nodes_in) {
        BestEdge be = find_intra_best(G, nodes, coreness, s, b, spent,
                                       upperbound, T1_self_edge, T2_upperbound,
                                       s_cand, delta_tactic, counter);
        if (be.valid) {
            cache.intra_best[cid] = be;
        }
    }

    return cache;
}

// -------------------------------------------------------------------------
// find_intra_best
// -------------------------------------------------------------------------
BestEdge find_intra_best(Graph& G, const std::unordered_set<int>& nodes_set,
                          Coreness& coreness, int s, int b, int spent,
                          std::unordered_map<int,double>& upperbound,
                          bool T1_self_edge, bool T2_upperbound,
                          int s_cand, const std::string& delta_tactic,
                          long long* counter) {
    std::vector<int> nodes_vec(nodes_set.begin(), nodes_set.end());
    std::unordered_map<int,int> empty_comp;

    if (T1_self_edge) {
        auto cands = make_candidate_nodes(G, nodes_vec, s, b - spent, coreness,
                                          upperbound, T2_upperbound, delta_tactic);
        if (T2_upperbound) {
            return iteration_nodes_upperbound(G, cands, s, b, spent, coreness,
                                               upperbound, s_cand,
                                               T1_self_edge, false,
                                               delta_tactic, empty_comp, BestEdge{},
                                               counter);
        } else {
            return iteration_nodes_no_upperbound(G, cands, s, b, spent, coreness,
                                                  s_cand, false, delta_tactic,
                                                  empty_comp, BestEdge{},
                                                  counter);
        }
    } else {
        if (T2_upperbound) {
            // Need to pass node set as "all nodes" - build v2 candidates from nodes_set
            std::unordered_map<int,double> local_ub = upperbound;
            std::vector<int> cands;
            for (int u : nodes_vec) {
                int cu = G.label[u] ? s : coreness[u].first;
                if (s - cu > b - spent) continue;
                if (!G.label[u]) {
                    local_ub[u] = Upperbound(G, u, coreness, s, delta_tactic);
                } else {
                    local_ub[u] = 0.0;
                }
                cands.push_back(u);
            }
            std::sort(cands.begin(), cands.end(),
                      [&](int a, int bb) { return local_ub[a] > local_ub[bb]; });
            // Merge local_ub back
            for (auto& [k, v] : local_ub) upperbound[k] = v;
            return iteration_nodes_upperbound(G, cands, s, b, spent, coreness,
                                               upperbound, s_cand,
                                               T1_self_edge, false,
                                               delta_tactic, empty_comp, BestEdge{},
                                               counter);
        } else {
            // Build candidate edges from nodes_set
            std::vector<int> non_s_core, s_core_local;
            for (int u : nodes_vec) {
                int cu = G.label[u] ? s : coreness[u].first;
                if (s - cu > b - spent) continue;
                if (!G.label[u]) non_s_core.push_back(u);
                else s_core_local.push_back(u);
            }
            std::vector<std::pair<int,int>> cand_edges;
            int nn = (int)non_s_core.size();
            for (int i = 0; i < nn; i++)
                for (int j = i+1; j < nn; j++)
                    cand_edges.push_back({non_s_core[i], non_s_core[j]});
            for (int u : non_s_core)
                for (int v : s_core_local)
                    cand_edges.push_back({u, v});
            return iteration_edges_no_upperbound(G, cand_edges, s, b, spent,
                                                  coreness, false, delta_tactic,
                                                  empty_comp, BestEdge{},
                                                  counter);
        }
    }
}

// -------------------------------------------------------------------------
// Helper: get nodes as vector from graph
// -------------------------------------------------------------------------
static std::vector<int> all_nodes(const Graph& G) {
    std::vector<int> v;
    v.reserve(G.n);
    for (int i = 0; i < G.n; i++) v.push_back(i);
    return v;
}

// -------------------------------------------------------------------------
// run (standard, no iter counting)
// -------------------------------------------------------------------------
RunResult run(Graph& G, int s, int b,
              bool T1_self_edge, bool T2_upperbound, bool T3_reuse,
              const std::string& delta_tactic) {
    RunResult result;
    result.FT = 0.0;
    result.UT = 0.0;
    result.total_follower = 0;

    Graph& G_prime = G; // We work in-place on a copy (caller passes copy)
    std::vector<int> all_n = all_nodes(G_prime);

    Coreness coreness(G_prime.n, {-1, -1});
    int s_core_num = calculate_s_core(G_prime, all_n, s, coreness);

    int spent = 0;
    int s_cand = -1;
    std::vector<int> non_s_core;

    if (T1_self_edge) {
        auto [ns, sc] = self_edge_pruning(G_prime);
        non_s_core = ns;
        s_cand = sc;
    }

    std::unordered_map<int,double> upperbound;
    CCCache cache;

    if (T3_reuse) {
        cache = build_initial_caches(G_prime, s, b, spent, coreness, upperbound,
                                      T1_self_edge, T2_upperbound, s_cand, delta_tactic);
    }

    while (spent < b) {
        BestEdge best{};
        BestEdge best_intra{};

        if (T3_reuse && !cache.intra_best.empty()) {
            for (auto& [cid, be] : cache.intra_best) {
                // Tie-break must match the single-pass scan (e_canon < best_edge):
                // higher FR wins; on equal FR, the smaller canonical edge wins.
                // (Previously first-wins on FR only, which diverged from TTF/FTF/TFF.)
                if (!best_intra.valid || be.FR > best_intra.FR ||
                    (be.FR == best_intra.FR && be.edge < best_intra.edge)) {
                    best_intra = be;
                }
            }
            if (best_intra.valid) best = best_intra;
        }

        if (T1_self_edge) {
            auto cands = make_candidate_nodes(G_prime, non_s_core, s, b - spent,
                                              coreness, upperbound, T2_upperbound,
                                              delta_tactic);
            if (T2_upperbound) {
                best = iteration_nodes_upperbound(G_prime, cands, s, b, spent,
                                                   coreness, upperbound, s_cand,
                                                   T1_self_edge, T3_reuse,
                                                   delta_tactic, cache.comp_of, best);
            } else {
                best = iteration_nodes_no_upperbound(G_prime, cands, s, b, spent,
                                                      coreness, s_cand, T3_reuse,
                                                      delta_tactic, cache.comp_of, best);
            }
        } else {
            if (T2_upperbound) {
                auto cands = make_candidate_nodes_v2(G_prime, s, b - spent,
                                                      coreness, upperbound,
                                                      T2_upperbound, delta_tactic);
                best = iteration_nodes_upperbound(G_prime, cands, s, b, spent,
                                                   coreness, upperbound, s_cand,
                                                   T1_self_edge, T3_reuse,
                                                   delta_tactic, cache.comp_of, best);
            } else {
                auto cand_edges = make_candidate_edges(G_prime, s, b - spent, coreness);
                best = iteration_edges_no_upperbound(G_prime, cand_edges, s, b, spent,
                                                      coreness, T3_reuse, delta_tactic,
                                                      cache.comp_of, best);
            }
        }

        if (!best.valid) break;

        auto [u, v] = best.edge;
        int best_delta = best.delta;

        // Check if budget exceeded (from reuse cache)
        if (best_delta > b - spent) {
            if (T3_reuse && best.edge == best_intra.edge && best.delta == best_intra.delta) {
                // Find which CC this came from
                int c = -1;
                if (G_prime.label[u]) {
                    auto it = cache.comp_of.find(v);
                    if (it != cache.comp_of.end()) c = it->second;
                } else {
                    auto it = cache.comp_of.find(u);
                    if (it != cache.comp_of.end()) c = it->second;
                }
                if (c != -1) {
                    cache.intra_best.erase(c);
                    auto it = cache.nodes_in.find(c);
                    if (it != cache.nodes_in.end()) {
                        BestEdge be2 = find_intra_best(G_prime, it->second, coreness,
                                                        s, b, spent, upperbound,
                                                        T1_self_edge, T2_upperbound,
                                                        s_cand, delta_tactic);
                        if (be2.valid) cache.intra_best[c] = be2;
                    }
                }
            }
            continue;
        }

        // Anchor best edge - resolve self-edge
        int eu = u, ev = v;
        if (eu == ev) {
            ev = s_cand;
        }

        // Update graph
        if (G_prime.has_edge(eu, ev)) {
            G_prime.increase_weight(eu, ev, best_delta);
        } else {
            G_prime.add_edge(eu, ev, best_delta);
        }

        spent += best_delta;
        result.total_follower += best.num_followers;
        result.A.push_back({best.edge, best_delta, best.FR, best.num_followers});

        // Update coreness
        if (T3_reuse) {
            bool is_intra = (best.edge == best_intra.edge && best.delta == best_intra.delta && best_intra.valid);

            if (is_intra) {
                // Same CC update
                int c = -1;
                if (G_prime.label[eu]) {
                    auto it = cache.comp_of.find(ev);
                    if (it != cache.comp_of.end()) c = it->second;
                } else {
                    auto it = cache.comp_of.find(eu);
                    if (it != cache.comp_of.end()) c = it->second;
                }
                if (c == -1) {
                    // fallback: global recalculation
                    coreness.assign(G_prime.n, {-1,-1});
                    s_core_num = calculate_s_core(G_prime, all_nodes(G_prime), s, coreness);
                    continue;
                }

                auto it = cache.nodes_in.find(c);
                if (it == cache.nodes_in.end()) continue;

                std::vector<int> cc_nodes_vec(it->second.begin(), it->second.end());
                s_core_num = calculate_s_core(G_prime, cc_nodes_vec, s, coreness);

                // Find non-s-core nodes still in CC
                std::unordered_set<int> new_nodes;
                for (int nd : cc_nodes_vec) {
                    if (!G_prime.label[nd]) new_nodes.insert(nd);
                }

                cache.intra_best.erase(c);

                if (new_nodes.empty()) {
                    cache.nodes_in.erase(c);
                } else {
                    if (T1_self_edge) {
                        cache.nodes_in[c] = new_nodes;
                    }
                    BestEdge be2 = find_intra_best(G_prime, cache.nodes_in[c], coreness,
                                                    s, b, spent, upperbound,
                                                    T1_self_edge, T2_upperbound,
                                                    s_cand, delta_tactic);
                    if (be2.valid) cache.intra_best[c] = be2;
                }
            } else {
                // Inter-CC: union
                int c1 = -1, c2 = -1;
                {
                    auto it = cache.comp_of.find(eu);
                    if (it != cache.comp_of.end()) c1 = it->second;
                }
                {
                    auto it = cache.comp_of.find(ev);
                    if (it != cache.comp_of.end()) c2 = it->second;
                }

                std::unordered_set<int> union_nodes;
                if (c1 != -1) {
                    auto it = cache.nodes_in.find(c1);
                    if (it != cache.nodes_in.end())
                        union_nodes.insert(it->second.begin(), it->second.end());
                }
                if (c2 != -1 && c2 != c1) {
                    auto it = cache.nodes_in.find(c2);
                    if (it != cache.nodes_in.end())
                        union_nodes.insert(it->second.begin(), it->second.end());
                }

                std::vector<int> union_vec(union_nodes.begin(), union_nodes.end());
                s_core_num = calculate_s_core(G_prime, union_vec, s, coreness);

                std::unordered_set<int> new_nodes;
                for (int nd : union_vec) {
                    if (!G_prime.label[nd]) new_nodes.insert(nd);
                }

                // Delete old caches
                cache.intra_best.erase(c1);
                cache.intra_best.erase(c2);
                if (c1 != -1) cache.nodes_in.erase(c1);
                if (c2 != -1 && c2 != c1) cache.nodes_in.erase(c2);

                if (!new_nodes.empty()) {
                    int new_c = 0;
                    if (!cache.nodes_in.empty()) {
                        for (auto& [k, v2] : cache.nodes_in) {
                            if (k >= new_c) new_c = k + 1;
                        }
                    }
                    cache.nodes_in[new_c] = new_nodes;
                    for (int nd : new_nodes) {
                        cache.comp_of[nd] = new_c;
                    }
                    BestEdge be2 = find_intra_best(G_prime, new_nodes, coreness,
                                                    s, b, spent, upperbound,
                                                    T1_self_edge, T2_upperbound,
                                                    s_cand, delta_tactic);
                    if (be2.valid) cache.intra_best[new_c] = be2;
                }
            }
        } else {
            // Global recalculation
            coreness.assign(G_prime.n, {-1,-1});
            s_core_num = calculate_s_core(G_prime, all_nodes(G_prime), s, coreness);
        }
    }

    result.G_prime = G_prime;
    return result;
}

// -------------------------------------------------------------------------
// run_iter (with GLOBAL_CNT)
// -------------------------------------------------------------------------
IterRunResult run_iter(Graph& G, int s, int b,
                        bool T1_self_edge, bool T2_upperbound, bool T3_reuse,
                        const std::string& delta_tactic) {
    IterRunResult result;
    result.FT = 0.0;
    result.UT = 0.0;
    result.total_follower = 0;
    result.global_cnt = 0;

    Graph& G_prime = G;
    std::vector<int> all_n = all_nodes(G_prime);

    Coreness coreness(G_prime.n, {-1, -1});
    int s_core_num = calculate_s_core(G_prime, all_n, s, coreness);

    int spent = 0;
    int s_cand = -1;
    std::vector<int> non_s_core;

    if (T1_self_edge) {
        auto [ns, sc] = self_edge_pruning(G_prime);
        non_s_core = ns;
        s_cand = sc;
    }

    std::unordered_map<int,double> upperbound;
    CCCache cache;

    if (T3_reuse) {
        cache = build_initial_caches(G_prime, s, b, spent, coreness, upperbound,
                                      T1_self_edge, T2_upperbound, s_cand, delta_tactic,
                                      &result.global_cnt);
    }

    while (spent < b) {
        BestEdge best{};
        BestEdge best_intra{};

        if (T3_reuse && !cache.intra_best.empty()) {
            for (auto& [cid, be] : cache.intra_best) {
                // Tie-break must match the single-pass scan (e_canon < best_edge):
                // higher FR wins; on equal FR, the smaller canonical edge wins.
                // (Previously first-wins on FR only, which diverged from TTF/FTF/TFF.)
                if (!best_intra.valid || be.FR > best_intra.FR ||
                    (be.FR == best_intra.FR && be.edge < best_intra.edge)) {
                    best_intra = be;
                }
            }
            if (best_intra.valid) best = best_intra;
        }

        if (T1_self_edge) {
            auto cands = make_candidate_nodes(G_prime, non_s_core, s, b - spent,
                                              coreness, upperbound, T2_upperbound,
                                              delta_tactic);
            if (T2_upperbound) {
                best = iteration_nodes_upperbound(G_prime, cands, s, b, spent,
                                                   coreness, upperbound, s_cand,
                                                   T1_self_edge, T3_reuse,
                                                   delta_tactic, cache.comp_of, best,
                                                   &result.global_cnt);
            } else {
                best = iteration_nodes_no_upperbound(G_prime, cands, s, b, spent,
                                                      coreness, s_cand, T3_reuse,
                                                      delta_tactic, cache.comp_of, best,
                                                      &result.global_cnt);
            }
        } else {
            if (T2_upperbound) {
                auto cands = make_candidate_nodes_v2(G_prime, s, b - spent,
                                                      coreness, upperbound,
                                                      T2_upperbound, delta_tactic);
                best = iteration_nodes_upperbound(G_prime, cands, s, b, spent,
                                                   coreness, upperbound, s_cand,
                                                   T1_self_edge, T3_reuse,
                                                   delta_tactic, cache.comp_of, best,
                                                   &result.global_cnt);
            } else {
                auto cand_edges = make_candidate_edges(G_prime, s, b - spent, coreness);
                best = iteration_edges_no_upperbound(G_prime, cand_edges, s, b, spent,
                                                      coreness, T3_reuse, delta_tactic,
                                                      cache.comp_of, best,
                                                      &result.global_cnt);
            }
        }

        if (!best.valid) break;

        auto [u, v] = best.edge;
        int best_delta = best.delta;

        if (best_delta > b - spent) {
            if (T3_reuse && best_intra.valid &&
                best.edge == best_intra.edge && best.delta == best_intra.delta) {
                int c = -1;
                if (G_prime.label[u]) {
                    auto it = cache.comp_of.find(v);
                    if (it != cache.comp_of.end()) c = it->second;
                } else {
                    auto it = cache.comp_of.find(u);
                    if (it != cache.comp_of.end()) c = it->second;
                }
                if (c != -1) {
                    cache.intra_best.erase(c);
                    auto it = cache.nodes_in.find(c);
                    if (it != cache.nodes_in.end()) {
                        BestEdge be2 = find_intra_best(G_prime, it->second, coreness,
                                                        s, b, spent, upperbound,
                                                        T1_self_edge, T2_upperbound,
                                                        s_cand, delta_tactic,
                                                        &result.global_cnt);
                        if (be2.valid) cache.intra_best[c] = be2;
                    }
                }
            }
            continue;
        }

        // Anchor with self/reinforce tracking
        bool self_flag = false;
        bool reinforce_flag = false;

        int eu = u, ev = v;
        if (eu == ev) {
            ev = s_cand;
            self_flag = true;
        }

        if (G_prime.has_edge(eu, ev)) {
            G_prime.increase_weight(eu, ev, best_delta);
            reinforce_flag = true;
        } else {
            G_prime.add_edge(eu, ev, best_delta);
        }

        spent += best_delta;
        result.total_follower += best.num_followers;
        result.A.push_back({best.edge, best_delta, best.FR, best.num_followers, self_flag, reinforce_flag});

        if (T3_reuse) {
            bool is_intra = (best_intra.valid && best.edge == best_intra.edge && best.delta == best_intra.delta);

            if (is_intra) {
                int c = -1;
                if (G_prime.label[eu]) {
                    auto it = cache.comp_of.find(ev);
                    if (it != cache.comp_of.end()) c = it->second;
                } else {
                    auto it = cache.comp_of.find(eu);
                    if (it != cache.comp_of.end()) c = it->second;
                }
                if (c == -1) {
                    coreness.assign(G_prime.n, {-1,-1});
                    s_core_num = calculate_s_core(G_prime, all_nodes(G_prime), s, coreness);
                    continue;
                }

                auto it = cache.nodes_in.find(c);
                if (it == cache.nodes_in.end()) continue;
                std::vector<int> cc_nodes_vec(it->second.begin(), it->second.end());
                s_core_num = calculate_s_core(G_prime, cc_nodes_vec, s, coreness);

                std::unordered_set<int> new_nodes;
                for (int nd : cc_nodes_vec) {
                    if (!G_prime.label[nd]) new_nodes.insert(nd);
                }

                cache.intra_best.erase(c);

                if (new_nodes.empty()) {
                    cache.nodes_in.erase(c);
                } else {
                    if (T1_self_edge) cache.nodes_in[c] = new_nodes;
                    BestEdge be2 = find_intra_best(G_prime, cache.nodes_in[c], coreness,
                                                    s, b, spent, upperbound,
                                                    T1_self_edge, T2_upperbound,
                                                    s_cand, delta_tactic,
                                                    &result.global_cnt);
                    if (be2.valid) cache.intra_best[c] = be2;
                }
            } else {
                int c1 = -1, c2 = -1;
                {
                    auto it = cache.comp_of.find(eu);
                    if (it != cache.comp_of.end()) c1 = it->second;
                }
                {
                    auto it = cache.comp_of.find(ev);
                    if (it != cache.comp_of.end()) c2 = it->second;
                }

                std::unordered_set<int> union_nodes;
                if (c1 != -1) {
                    auto it = cache.nodes_in.find(c1);
                    if (it != cache.nodes_in.end())
                        union_nodes.insert(it->second.begin(), it->second.end());
                }
                if (c2 != -1 && c2 != c1) {
                    auto it = cache.nodes_in.find(c2);
                    if (it != cache.nodes_in.end())
                        union_nodes.insert(it->second.begin(), it->second.end());
                }

                std::vector<int> union_vec(union_nodes.begin(), union_nodes.end());
                s_core_num = calculate_s_core(G_prime, union_vec, s, coreness);

                std::unordered_set<int> new_nodes;
                for (int nd : union_vec) {
                    if (!G_prime.label[nd]) new_nodes.insert(nd);
                }

                cache.intra_best.erase(c1);
                cache.intra_best.erase(c2);
                if (c1 != -1) cache.nodes_in.erase(c1);
                if (c2 != -1 && c2 != c1) cache.nodes_in.erase(c2);

                if (!new_nodes.empty()) {
                    int new_c = 0;
                    if (!cache.nodes_in.empty()) {
                        for (auto& [k, vv] : cache.nodes_in) {
                            if (k >= new_c) new_c = k + 1;
                        }
                    }
                    cache.nodes_in[new_c] = new_nodes;
                    for (int nd : new_nodes) {
                        cache.comp_of[nd] = new_c;
                    }
                    BestEdge be2 = find_intra_best(G_prime, new_nodes, coreness,
                                                    s, b, spent, upperbound,
                                                    T1_self_edge, T2_upperbound,
                                                    s_cand, delta_tactic,
                                                    &result.global_cnt);
                    if (be2.valid) cache.intra_best[new_c] = be2;
                }
            }
        } else {
            coreness.assign(G_prime.n, {-1,-1});
            s_core_num = calculate_s_core(G_prime, all_nodes(G_prime), s, coreness);
        }
    }

    result.G_prime = G_prime;
    return result;
}
