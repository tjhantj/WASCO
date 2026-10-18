#include "experiment.h"
#include "functions.h"
#include "exp_func.h"
#include <algorithm>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <stdexcept>

// -------------------------------------------------------------------------
// Connected components
// -------------------------------------------------------------------------
static std::vector<std::unordered_set<int>> connected_components_subset(
    const Graph& G, const std::unordered_set<int>& node_set) {
    std::unordered_set<int> visited;
    std::vector<std::unordered_set<int>> components;

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
    for (int v = 0; v < G.n; v++) {
        if (G.label[v]) { s_cand = v; break; }
    }

    std::vector<std::unordered_set<int>> ccs;
    if (T1_self_edge) {
        std::unordered_set<int> non_s_core_nodes;
        for (int v = 0; v < G.n; v++) {
            if (!G.label[v]) non_s_core_nodes.insert(v);
        }
        ccs = connected_components_subset(G, non_s_core_nodes);
    } else {
        ccs = connected_components_all(G);
    }

    int cid = 0;
    for (auto& nodes : ccs) {
        for (int v : nodes) {
            cache.comp_of[v] = cid;
        }
        cache.nodes_in[cid] = std::move(nodes);
        cid++;
    }
    cache.next_cid = cid;

    for (auto& [cid2, nodes] : cache.nodes_in) {
        BestEdge be = find_intra_best(G, nodes, coreness, s, b, spent,
                                       upperbound, T1_self_edge, T2_upperbound,
                                       s_cand, delta_tactic, counter);
        if (be.valid) {
            be.region = cid2;
            cache.intra_best[cid2] = be;
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
            for (auto& [k, v] : local_ub) upperbound[k] = v;
            return iteration_nodes_upperbound(G, cands, s, b, spent, coreness,
                                               upperbound, s_cand,
                                               T1_self_edge, false,
                                               delta_tactic, empty_comp, BestEdge{},
                                               counter);
        } else {
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
// all_nodes
// -------------------------------------------------------------------------
static std::vector<int> all_nodes(const Graph& G) {
    std::vector<int> v;
    v.reserve(G.n);
    for (int i = 0; i < G.n; i++) v.push_back(i);
    return v;
}

// -------------------------------------------------------------------------
// min_core_id
// -------------------------------------------------------------------------
static int min_core_id(const Graph& G) {
    for (int i = 0; i < G.n; i++) {
        if (G.label[i]) return i;
    }
    return -1;
}

static void erase_region(CCCache& cache, int cid) {
    if (cid == -1) return;
    auto it = cache.nodes_in.find(cid);
    if (it != cache.nodes_in.end()) {
        for (int nd : it->second) {
            auto m = cache.comp_of.find(nd);
            if (m != cache.comp_of.end() && m->second == cid) cache.comp_of.erase(m);
        }
        cache.nodes_in.erase(it);
    }
    cache.intra_best.erase(cid);
}

// -------------------------------------------------------------------------
// run_impl
// -------------------------------------------------------------------------
static IterRunResult run_impl(Graph& G, int s, int b,
                              bool T1_self_edge, bool T2_upperbound, bool T3_reuse,
                              const std::string& delta_tactic,
                              long long* counter) {
    IterRunResult result;
    result.FT = 0.0;
    result.UT = 0.0;
    result.total_follower = 0;
    result.global_cnt = 0;

    Graph& G_prime = G;

    Coreness coreness(G_prime.n, {-1, -1});
    result.initial_s_core = calculate_s_core(G_prime, all_nodes(G_prime), s, coreness);

    int spent = 0;
    int s_cand = min_core_id(G_prime);
    std::vector<int> non_s_core;

    if (T1_self_edge) {
        auto [ns, sc] = self_edge_pruning(G_prime);
        non_s_core = ns;
        s_cand = sc;
    }

    std::unordered_map<int,double> upperbound;
    CCCache cache;

    auto recompute_intra = [&](int cid) {
        cache.intra_best.erase(cid);
        auto it = cache.nodes_in.find(cid);
        if (it == cache.nodes_in.end()) return;
        BestEdge be = find_intra_best(G_prime, it->second, coreness, s, b, spent,
                                      upperbound, T1_self_edge, T2_upperbound,
                                      s_cand, delta_tactic, counter);
        if (be.valid) {
            be.region = cid;
            cache.intra_best[cid] = be;
        }
    };

    auto rebuild_all = [&]() {
        coreness.assign(G_prime.n, {-1,-1});
        calculate_s_core(G_prime, all_nodes(G_prime), s, coreness);
        s_cand = min_core_id(G_prime);
        if (T3_reuse) {
            cache = build_initial_caches(G_prime, s, b, spent, coreness, upperbound,
                                         T1_self_edge, T2_upperbound, s_cand,
                                         delta_tactic, counter);
        }
    };

    if (T3_reuse) {
        cache = build_initial_caches(G_prime, s, b, spent, coreness, upperbound,
                                      T1_self_edge, T2_upperbound, s_cand, delta_tactic,
                                      counter);
    }

    while (spent < b) {
        int budget_left = b - spent;
        BestEdge best{};

        if (T3_reuse) {
            while (true) {
                BestEdge top{};
                for (auto& [cid, be] : cache.intra_best) {
                    if (is_better(be, top)) top = be;
                }
                if (top.valid && top.delta > budget_left) {
                    recompute_intra(top.region);
                    continue;
                }
                best = top;
                break;
            }
        }

        if (T1_self_edge) {
            auto cands = make_candidate_nodes(G_prime, non_s_core, s, budget_left,
                                              coreness, upperbound, T2_upperbound,
                                              delta_tactic);
            if (T2_upperbound) {
                best = iteration_nodes_upperbound(G_prime, cands, s, b, spent,
                                                   coreness, upperbound, s_cand,
                                                   T1_self_edge, T3_reuse,
                                                   delta_tactic, cache.comp_of, best,
                                                   counter);
            } else {
                best = iteration_nodes_no_upperbound(G_prime, cands, s, b, spent,
                                                      coreness, s_cand, T3_reuse,
                                                      delta_tactic, cache.comp_of, best,
                                                      counter);
            }
        } else {
            if (T2_upperbound) {
                auto cands = make_candidate_nodes_v2(G_prime, s, budget_left,
                                                      coreness, upperbound,
                                                      T2_upperbound, delta_tactic);
                best = iteration_nodes_upperbound(G_prime, cands, s, b, spent,
                                                   coreness, upperbound, s_cand,
                                                   T1_self_edge, T3_reuse,
                                                   delta_tactic, cache.comp_of, best,
                                                   counter);
            } else {
                auto cand_edges = make_candidate_edges(G_prime, s, budget_left, coreness);
                best = iteration_edges_no_upperbound(G_prime, cand_edges, s, b, spent,
                                                      coreness, T3_reuse, delta_tactic,
                                                      cache.comp_of, best,
                                                      counter);
            }
        }

        if (!best.valid) break;
        if (best.delta <= 0 || best.delta > budget_left) {
            throw std::logic_error("run: selected edge does not fit the remaining budget");
        }

        auto [eu, ev] = best.edge;
        bool reinforce_flag = G_prime.has_edge(eu, ev);
        if (reinforce_flag) {
            G_prime.increase_weight(eu, ev, best.delta);
        } else {
            G_prime.add_edge(eu, ev, best.delta);
        }

        spent += best.delta;
        result.total_follower += best.num_followers;
        result.A.push_back({best.edge, best.delta, best.FR, best.num_followers,
                            best.lower_bound, best.core_rep, reinforce_flag});

        if (!T3_reuse) {
            coreness.assign(G_prime.n, {-1,-1});
            calculate_s_core(G_prime, all_nodes(G_prime), s, coreness);
            s_cand = min_core_id(G_prime);
            continue;
        }

        std::vector<int> dirty;
        bool resync = false;

        if (best.region != -1) {
            int c = best.region;
            auto it = cache.nodes_in.find(c);
            if (it == cache.nodes_in.end()) {
                resync = true;
            } else {
                std::vector<int> cc_nodes_vec(it->second.begin(), it->second.end());
                calculate_s_core(G_prime, cc_nodes_vec, s, coreness);

                std::unordered_set<int> new_nodes;
                for (int nd : cc_nodes_vec) {
                    if (!G_prime.label[nd]) new_nodes.insert(nd);
                }

                cache.intra_best.erase(c);

                if (new_nodes.empty()) {
                    erase_region(cache, c);
                } else {
                    if (T1_self_edge) {
                        for (int nd : cc_nodes_vec) {
                            if (G_prime.label[nd]) cache.comp_of.erase(nd);
                        }
                        cache.nodes_in[c] = std::move(new_nodes);
                    }
                    dirty.push_back(c);
                }
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

            if (union_nodes.empty()) {
                resync = true;
            } else {
                std::vector<int> union_vec(union_nodes.begin(), union_nodes.end());
                calculate_s_core(G_prime, union_vec, s, coreness);

                std::unordered_set<int> new_nodes;
                for (int nd : union_vec) {
                    if (!G_prime.label[nd]) new_nodes.insert(nd);
                }

                erase_region(cache, c1);
                if (c2 != c1) erase_region(cache, c2);

                if (!new_nodes.empty()) {
                    int new_c = cache.next_cid++;
                    for (int nd : new_nodes) {
                        cache.comp_of[nd] = new_c;
                    }
                    cache.nodes_in[new_c] = std::move(new_nodes);
                    dirty.push_back(new_c);
                }
            }
        }

        if (resync) {
            rebuild_all();
            continue;
        }

        int new_s_cand = min_core_id(G_prime);
        if (new_s_cand != s_cand) {
            s_cand = new_s_cand;
            if (T1_self_edge) {
                dirty.clear();
                for (auto& [cid, nodes] : cache.nodes_in) dirty.push_back(cid);
            }
        }

        for (int cid : dirty) recompute_intra(cid);
    }

    result.spent = spent;
    result.G_prime = G_prime;
    return result;
}

// -------------------------------------------------------------------------
// run
// -------------------------------------------------------------------------
RunResult run(Graph& G, int s, int b,
              bool T1_self_edge, bool T2_upperbound, bool T3_reuse,
              const std::string& delta_tactic) {
    IterRunResult r = run_impl(G, s, b, T1_self_edge, T2_upperbound, T3_reuse,
                               delta_tactic, nullptr);
    RunResult result;
    result.A = std::move(r.A);
    result.FT = r.FT;
    result.UT = r.UT;
    result.G_prime = std::move(r.G_prime);
    result.total_follower = r.total_follower;
    result.initial_s_core = r.initial_s_core;
    result.spent = r.spent;
    return result;
}

// -------------------------------------------------------------------------
// run_iter
// -------------------------------------------------------------------------
IterRunResult run_iter(Graph& G, int s, int b,
                        bool T1_self_edge, bool T2_upperbound, bool T3_reuse,
                        const std::string& delta_tactic) {
    long long cnt = 0;
    IterRunResult result = run_impl(G, s, b, T1_self_edge, T2_upperbound, T3_reuse,
                                    delta_tactic, &cnt);
    result.global_cnt = cnt;
    return result;
}
