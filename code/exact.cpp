#include "exact.h"
#include "functions.h"
#include <algorithm>
#include <atomic>
#include <iostream>
#include <map>
#include <thread>
#include <vector>

// -------------------------------------------------------------------------
// Exhaustive search over all combinations with replacement of size b from
// the candidate edges.
// -------------------------------------------------------------------------

namespace {

struct Candidate {
    int a;
    int c;  // -1: edge to the s-core representative
};

struct PartBest {
    int gain = -1;
    std::vector<int> combo;
};

// -------------------------------------------------------------------------
// Evaluator
// -------------------------------------------------------------------------
class Evaluator {
public:
    Evaluator(int nn, int s, const std::vector<int>& deg0, const std::vector<int>& w0,
              const std::vector<Candidate>& cand, int b)
        : nn_(nn), s_(s), b_(b), deg_(deg0), w_(w0), cand_(cand),
          d_(nn), alive_(nn), stack_(nn), cur_(b) {}

    PartBest run_part(int first) {
        best_ = PartBest();
        cur_[0] = first;
        apply(first, +1);
        rec(1, first);
        apply(first, -1);
        return best_;
    }

private:
    int nn_, s_, b_;
    std::vector<int> deg_;
    std::vector<int> w_;
    const std::vector<Candidate>& cand_;
    std::vector<int> d_, alive_, stack_, cur_;
    PartBest best_;

    void apply(int idx, int dw) {
        const Candidate& e = cand_[idx];
        deg_[e.a] += dw;
        if (e.c >= 0) {
            deg_[e.c] += dw;
            w_[e.a * nn_ + e.c] += dw;
            w_[e.c * nn_ + e.a] += dw;
        }
    }

    void rec(int depth, int start) {
        if (depth == b_) {
            int gain = survivors();
            if (gain > best_.gain) {
                best_.gain = gain;
                best_.combo = cur_;
            }
            return;
        }
        int m = (int)cand_.size();
        for (int i = start; i < m; i++) {
            cur_[depth] = i;
            apply(i, +1);
            rec(depth + 1, i);
            apply(i, -1);
        }
    }

    int survivors() {
        int top = 0;
        for (int i = 0; i < nn_; i++) {
            d_[i] = deg_[i];
            alive_[i] = d_[i] >= s_;
            if (!alive_[i]) stack_[top++] = i;
        }
        int left = nn_ - top;
        while (top > 0 && left > 0) {
            int x = stack_[--top];
            const int* row = &w_[x * nn_];
            for (int y = 0; y < nn_; y++) {
                if (!alive_[y] || row[y] == 0) continue;
                d_[y] -= row[y];
                if (d_[y] < s_) {
                    alive_[y] = 0;
                    stack_[top++] = y;
                    left--;
                }
            }
        }
        return left;
    }
};

}  // namespace

// -------------------------------------------------------------------------
// exact_run
// -------------------------------------------------------------------------
ExactResult exact_run(Graph& G, int s, int b, int threads) {
    ExactResult result;

    Graph G_prime = G.copy();
    std::vector<int> all_n;
    for (int i = 0; i < G_prime.n; i++) all_n.push_back(i);

    Coreness coreness(G_prime.n, {-1,-1});
    int s_core_num = calculate_s_core(G_prime, all_n, s, coreness);
    int initial_size = s_core_num;

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

    int nn = (int)non_s_core.size();
    std::vector<std::pair<int,int>> candidate_edges;
    std::vector<Candidate> cand;
    for (int i = 0; i < nn; i++) {
        for (int j = i+1; j < nn; j++) {
            candidate_edges.push_back({non_s_core[i], non_s_core[j]});
            cand.push_back({i, j});
        }
    }
    for (int i = 0; i < nn; i++) {
        candidate_edges.push_back({non_s_core[i], s_cand});
        cand.push_back({i, -1});
    }

    std::vector<int> pos(G.n, -1);
    for (int i = 0; i < nn; i++) pos[non_s_core[i]] = i;
    std::vector<int> deg0(nn, 0), w0((size_t)nn * nn, 0);
    for (int i = 0; i < nn; i++) {
        for (auto& [nb, w] : G.adj[non_s_core[i]]) {
            deg0[i] += w;
            if (pos[nb] >= 0) w0[(size_t)i * nn + pos[nb]] = w;
        }
    }

    int best_gain = -1;
    std::vector<int> best_indices;

    int m = (int)cand.size();
    if (b == 0) {
        best_gain = 0;
    } else if (m > 0) {
        std::vector<PartBest> parts(m);
        std::atomic<int> next(0);
        auto worker = [&]() {
            Evaluator ev(nn, s, deg0, w0, cand, b);
            for (int f = next++; f < m; f = next++) parts[f] = ev.run_part(f);
        };
        int nt = std::max(1, std::min(threads, m));
        std::vector<std::thread> pool;
        for (int t = 1; t < nt; t++) pool.emplace_back(worker);
        worker();
        for (auto& th : pool) th.join();

        for (int f = 0; f < m; f++) {
            if (parts[f].gain > best_gain) {
                best_gain = parts[f].gain;
                best_indices = parts[f].combo;
            }
        }
    }

    std::map<std::pair<int,int>, int> best_combo;
    for (int idx : best_indices) best_combo[candidate_edges[idx]]++;

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

    if (!best_combo.empty() && result.followers != best_gain) {
        std::cerr << "Warning: exact gain " << best_gain << " != followers "
                  << result.followers << std::endl;
    }

    return result;
}
