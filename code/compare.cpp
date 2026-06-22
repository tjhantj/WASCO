#include "compare.h"
#include "functions.h"
#include <algorithm>
#include <random>
#include <stdexcept>

CompareResult compare_run(Graph& G, int s, int b,
                           const std::string& compare_tactic,
                           const std::string& delta_tactic) {
    CompareResult result;
    result.budget_used = 0;
    result.total_followers = 0;

    Graph G_prime = G.copy();
    std::vector<int> all_n;
    for (int i = 0; i < G_prime.n; i++) all_n.push_back(i);

    Coreness coreness(G_prime.n, {-1,-1});
    int s_core_num = calculate_s_core(G_prime, all_n, s, coreness);
    int initial_size = s_core_num;

    // Partition
    std::vector<int> non_s_core;
    int s_cand = -1;
    for (int nd = 0; nd < G_prime.n; nd++) {
        if (!G_prime.label[nd]) {
            non_s_core.push_back(nd);
        } else {
            if (s_cand == -1) s_cand = nd;
        }
    }
    if (s_cand == -1 || b <= 0) {
        result.final_size = 0;
        return result;
    }

    // Candidates: only (s-1)-shell nodes
    std::vector<int> candidates;
    for (int u : non_s_core) {
        if (coreness[u].first == s - 1) {
            candidates.push_back(u);
        }
    }

    if (candidates.empty()) {
        Coreness fc(G_prime.n, {-1,-1});
        result.final_size = calculate_s_core(G_prime, all_n, s, fc);
        result.total_followers = result.final_size - initial_size;
        return result;
    }

    // Score functions
    auto score_degree = [&](int u) -> double {
        return (double)G_prime.adj[u].size();
    };
    auto score_weight_sum = [&](int u) -> double {
        double ws = 0;
        for (auto& [nb, w] : G_prime.adj[u]) ws += w;
        return ws;
    };
    auto score_high_layer_degree = [&](int u) -> double {
        double cnt = 0;
        for (auto& [nb, w] : G_prime.adj[u]) {
            auto cu = coreness[u];
            auto cnb = G_prime.label[nb] ? std::make_pair(s, 0) : coreness[nb];
            if (cnb > cu) cnt++;
        }
        return cnt;
    };
    auto score_high_layer_weight = [&](int u) -> double {
        double ws = 0;
        for (auto& [nb, w] : G_prime.adj[u]) {
            auto cu = coreness[u];
            auto cnb = G_prime.label[nb] ? std::make_pair(s, 0) : coreness[nb];
            if (cnb > cu) ws += w;
        }
        return ws;
    };

    // Sort candidates
    if (compare_tactic == "random") {
        std::random_device rd;
        std::mt19937 rng(rd());
        std::shuffle(candidates.begin(), candidates.end(), rng);
    } else if (compare_tactic == "degree") {
        std::stable_sort(candidates.begin(), candidates.end(),
                  [&](int a, int bb) { return score_degree(a) > score_degree(bb); });
    } else if (compare_tactic == "high_degree") {
        std::stable_sort(candidates.begin(), candidates.end(),
                  [&](int a, int bb) { return score_high_layer_degree(a) > score_high_layer_degree(bb); });
    } else if (compare_tactic == "weight_sum") {
        std::stable_sort(candidates.begin(), candidates.end(),
                  [&](int a, int bb) { return score_weight_sum(a) > score_weight_sum(bb); });
    } else if (compare_tactic == "high_weight_sum") {
        std::stable_sort(candidates.begin(), candidates.end(),
                  [&](int a, int bb) { return score_high_layer_weight(a) > score_high_layer_weight(bb); });
    }

    int top_count = std::min((int)candidates.size(), b);
    std::vector<int> top_nodes(candidates.begin(), candidates.begin() + top_count);

    // Connect each to s_cand with delta=1
    int budget_used = 0;
    for (int u : top_nodes) {
        if (budget_used + 1 > b) break;
        int delta = 1;

        if (G_prime.has_edge(u, s_cand)) {
            G_prime.increase_weight(u, s_cand, delta);
        } else {
            G_prime.add_edge(u, s_cand, delta);
        }
        budget_used += delta;
    }

    // Re-calculate final s-core size
    Coreness fc(G_prime.n, {-1,-1});
    result.final_size = calculate_s_core(G_prime, all_n, s, fc);
    result.budget_used = budget_used;
    result.total_followers = result.final_size - initial_size;

    return result;
}
