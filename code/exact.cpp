#include "exact.h"
#include "functions.h"
#include <algorithm>
#include <vector>
#include <map>

// -------------------------------------------------------------------------
// Generate all combinations with replacement of size b from candidate_edges
// This is exponential; only feasible for small b and small candidate sets.
// -------------------------------------------------------------------------
static void combinations_with_replacement_helper(
    const std::vector<std::pair<int,int>>& pool, int r, int start,
    std::vector<int>& current,
    std::vector<std::vector<int>>& result) {
    if ((int)current.size() == r) {
        result.push_back(current);
        return;
    }
    for (int i = start; i < (int)pool.size(); i++) {
        current.push_back(i);
        combinations_with_replacement_helper(pool, r, i, current, result);
        current.pop_back();
    }
}

ExactResult exact_run(Graph& G, int s, int b) {
    ExactResult result;

    Graph G_prime = G.copy();
    std::vector<int> all_n;
    for (int i = 0; i < G_prime.n; i++) all_n.push_back(i);

    Coreness coreness(G_prime.n, {-1,-1});
    int s_core_num = calculate_s_core(G_prime, all_n, s, coreness);
    int initial_size = s_core_num;

    // Find non-s-core and s_cand
    std::vector<int> non_s_core;
    int s_cand = -1;
    for (int nd = 0; nd < G_prime.n; nd++) {
        if (!G_prime.label[nd]) {
            non_s_core.push_back(nd);
        } else {
            if (s_cand == -1) s_cand = nd;
        }
    }
    if (s_cand == -1) {
        result.final_size = s_core_num;
        result.followers = 0;
        return result;
    }

    // Generate candidate edges
    std::vector<std::pair<int,int>> candidate_edges;
    int nn = (int)non_s_core.size();
    for (int i = 0; i < nn; i++) {
        for (int j = i+1; j < nn; j++) {
            candidate_edges.push_back({non_s_core[i], non_s_core[j]});
        }
    }
    for (int u : non_s_core) {
        candidate_edges.push_back({u, s_cand});
    }

    // Generate all combinations_with_replacement of size b
    std::vector<std::vector<int>> combos;
    std::vector<int> cur;
    combinations_with_replacement_helper(candidate_edges, b, 0, cur, combos);

    int best_gain = -1;
    std::map<std::pair<int,int>, int> best_combo;

    for (auto& combo_indices : combos) {
        // Count occurrences
        std::map<std::pair<int,int>, int> combo_map;
        for (int idx : combo_indices) {
            combo_map[candidate_edges[idx]]++;
        }

        // Build temp graph
        Graph G_temp = G.copy();
        for (auto& [e, delta] : combo_map) {
            auto [u, v] = e;
            if (G_temp.has_edge(u, v)) {
                G_temp.increase_weight(u, v, delta);
            } else {
                G_temp.add_edge(u, v, delta);
            }
        }

        Coreness temp_cor(G_temp.n, {-1,-1});
        std::vector<int> an;
        for (int i = 0; i < G_temp.n; i++) an.push_back(i);
        int new_size = calculate_s_core(G_temp, an, s, temp_cor);
        int gain = new_size - initial_size;

        if (gain > best_gain) {
            best_gain = gain;
            best_combo = combo_map;
        }
    }

    // Apply best combo
    if (!best_combo.empty()) {
        for (auto& [e, delta] : best_combo) {
            auto [u, v] = e;
            if (G_prime.has_edge(u, v)) {
                G_prime.increase_weight(u, v, delta);
            } else {
                G_prime.add_edge(u, v, delta);
            }
            result.A.push_back({e, delta});
        }
    }

    Coreness final_cor(G_prime.n, {-1,-1});
    std::vector<int> an2;
    for (int i = 0; i < G_prime.n; i++) an2.push_back(i);
    result.final_size = calculate_s_core(G_prime, an2, s, final_cor);
    result.followers = result.final_size - initial_size;

    return result;
}
