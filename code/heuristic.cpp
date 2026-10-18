#include "heuristic.h"
#include "exp_func.h"
#include <algorithm>
#include <climits>
#include <numeric>
#include <random>
#include <stdexcept>

// -------------------------------------------------------------------------
// Names
// -------------------------------------------------------------------------
bool parse_heuristic_name(const std::string& name, ScoreFamily& fam, IncrementPolicy& pol) {
    std::string base;
    if (name.size() > 3 && name.compare(name.size() - 3, 3, "-LB") == 0) {
        pol = IncrementPolicy::InitialLB;
        base = name.substr(0, name.size() - 3);
    } else if (name.size() > 4 && name.compare(name.size() - 4, 4, "-MEI") == 0) {
        pol = IncrementPolicy::CurrentMEI;
        base = name.substr(0, name.size() - 4);
    } else {
        return false;
    }
    if (base == "Deg") fam = ScoreFamily::Deg;
    else if (base == "Weight") fam = ScoreFamily::Weight;
    else if (base == "Deg-L") fam = ScoreFamily::DegL;
    else if (base == "Weight-L") fam = ScoreFamily::WeightL;
    else if (base == "Random") fam = ScoreFamily::Random;
    else return false;
    return true;
}

std::string heuristic_name(ScoreFamily fam, IncrementPolicy pol) {
    std::string base;
    switch (fam) {
        case ScoreFamily::Deg: base = "Deg"; break;
        case ScoreFamily::Weight: base = "Weight"; break;
        case ScoreFamily::DegL: base = "Deg-L"; break;
        case ScoreFamily::WeightL: base = "Weight-L"; break;
        case ScoreFamily::Random: base = "Random"; break;
    }
    return base + (pol == IncrementPolicy::InitialLB ? "-LB" : "-MEI");
}

// -------------------------------------------------------------------------
// Initial node scores
// -------------------------------------------------------------------------
std::vector<long long> compute_initial_node_scores(const Graph& G0, int, const Coreness& c0,
                                                   ScoreFamily fam, unsigned long long seed) {
    if (fam == ScoreFamily::Random) {
        std::vector<int> nc;
        for (int u = 0; u < G0.n; u++)
            if (!G0.label[u]) nc.push_back(u);
        std::mt19937_64 rng(seed);
        std::shuffle(nc.begin(), nc.end(), rng);
        std::vector<long long> q(G0.n, 0);
        for (size_t k = 0; k < nc.size(); k++) q[nc[k]] = (long long)k + 1;
        return q;
    }
    const std::pair<int,int> CORE = {INT_MAX, 0};
    auto cl = [&](int x) { return G0.label[x] ? CORE : c0[x]; };
    const bool layered = (fam == ScoreFamily::DegL || fam == ScoreFamily::WeightL);
    const bool weighted = (fam == ScoreFamily::Weight || fam == ScoreFamily::WeightL);

    std::vector<long long> q(G0.n, 0);
    for (int u = 0; u < G0.n; u++) {
        if (G0.label[u]) continue;
        auto cu = c0[u];
        long long sc = 0;
        for (auto& [nb, w] : G0.adj[u]) {
            auto cv = cl(nb);
            bool counted = layered ? (cu < cv) : true;
            if (counted) sc += weighted ? (long long)w : 1LL;
        }
        q[u] = sc;
    }
    return q;
}

// -------------------------------------------------------------------------
// FrozenPairOrder
// -------------------------------------------------------------------------
static long long choose2(long long n) { return n * (n - 1) / 2; }

FrozenPairOrder::FrozenPairOrder(const std::vector<long long>& qprime, const std::vector<char>& c0)
    : core0(c0) {
    int n = (int)qprime.size();
    vals = qprime;
    std::sort(vals.begin(), vals.end(), std::greater<long long>());
    vals.erase(std::unique(vals.begin(), vals.end()), vals.end());
    int D = (int)vals.size();
    all_ids.assign(D, {});
    core_ids.assign(D, {});
    nc_ids.assign(D, {});
    bucket_of.assign(n, -1);
    for (int v = 0; v < n; v++) {
        if (core0[v] && qprime[v] != 0) throw std::logic_error("FrozenPairOrder: core node with q' != 0");
        int i = (int)(std::lower_bound(vals.begin(), vals.end(), qprime[v], std::greater<long long>()) - vals.begin());
        bucket_of[v] = i;
        all_ids[i].push_back(v);
        (core0[v] ? core_ids[i] : nc_ids[i]).push_back(v);
    }
}

bool FrozenPairOrder::bucket_pair_nonempty(int i, int j) const {
    return bucket_pair_count(i, j) > 0;
}

long long FrozenPairOrder::bucket_pair_count(int i, int j) const {
    if (i != j) return (long long)all_ids[i].size() * (long long)all_ids[j].size();
    return choose2((long long)all_ids[i].size()) - choose2((long long)core_ids[i].size());
}

// -------------------------------------------------------------------------
// FrozenPairStream
// -------------------------------------------------------------------------
FrozenPairStream::FrozenPairStream(const FrozenPairOrder& ord) : ord_(ord) {}

void FrozenPairStream::rebuild_gen(bool capped, long long cap) {
    live_.clear();
    for (int i = 0; i < ord_.num_buckets(); i++) {
        if (!ord_.core_ids[i].empty() || !nc_f_[i].empty()) live_.push_back(i);
    }
    gen_heap_.clear();
    const int A = (int)live_.size();
    for (int ii = 0; ii < A; ii++) {
        int jj = ii;
        if (capped) {
            long long t = cap - ord_.vals[live_[ii]];
            jj = (int)(std::partition_point(live_.begin() + ii, live_.end(),
                                            [&](int bk) { return ord_.vals[bk] > t; }) - live_.begin());
            if (jj >= A) continue;
        }
        gen_heap_.push_back({ord_.vals[live_[ii]] + ord_.vals[live_[jj]], ii, jj});
    }
    std::make_heap(gen_heap_.begin(), gen_heap_.end());
}

bool FrozenPairStream::next_level() {
    level_pairs_.clear();
    const int A = (int)live_.size();
    while (!gen_heap_.empty()) {
        Row g = gen_heap_.front();
        if (!level_pairs_.empty() && g.sum != level_score_) break;
        std::pop_heap(gen_heap_.begin(), gen_heap_.end());
        gen_heap_.pop_back();
        if (g.jj + 1 < A) {
            gen_heap_.push_back({ord_.vals[live_[g.ii]] + ord_.vals[live_[g.jj + 1]], g.ii, g.jj + 1});
            std::push_heap(gen_heap_.begin(), gen_heap_.end());
        }
        int i = live_[g.ii], j = live_[g.jj];
        if (!ord_.bucket_pair_nonempty(i, j)) continue;
        if (level_pairs_.empty()) level_score_ = g.sum;
        level_pairs_.push_back({i, j});
    }
    return !level_pairs_.empty();
}

void FrozenPairStream::reset_filter(const OrderFilter& f) {
    f_ = f;
    const int D = ord_.num_buckets();
    nc_f_.assign(D, {});
    for (int i = 0; i < D; i++) {
        for (int p : ord_.nc_ids[i]) {
            if (f_.nc_partner_ok(p)) nc_f_[i].push_back(p);
        }
    }
    if (!started_) {
        started_ = true;
        rebuild_gen(false, 0);
        return;
    }
    if (!level_active_) return;
    const long long S = level_score_;
    rebuild_gen(true, S);
    if (!next_level()) {
        level_active_ = false;
        subs_.clear();
        heap_.clear();
        return;
    }
    if (level_score_ != S) has_last_ = false;
    build_level();
}

void FrozenPairStream::build_level() {
    subs_.clear();
    heap_.clear();
    auto add_sub = [&](Sub sb) {
        if (has_last_) {
            const PairKey last = last_key_;
            const int a = sb.anchor;
            if (sb.core_list) {
                auto it = std::partition_point(sb.core_list->begin(), sb.core_list->end(),
                                               [&](int p) { return canon(a, p) <= last; });
                sb.ci = (size_t)(it - sb.core_list->begin());
            }
            auto it = std::partition_point(sb.nc_list->begin() + sb.ni, sb.nc_list->end(),
                                           [&](int p) { return canon(a, p) <= last; });
            sb.ni = (size_t)(it - sb.nc_list->begin());
        }
        subs_.push_back(sb);
        push_head((int)subs_.size() - 1);
    };

    for (auto [i, j] : level_pairs_) {
        if (i == j) {
            const auto* cl = ord_.core_ids[i].empty() ? nullptr : &ord_.core_ids[i];
            const auto& ncl = nc_f_[i];
            for (int a : ord_.nc_ids[i]) {
                if (!f_.anchor_alive(a)) continue;
                size_t start = (size_t)(std::upper_bound(ncl.begin(), ncl.end(), a) - ncl.begin());
                add_sub(Sub{a, cl, 0, &ncl, start});
            }
        } else {
            int A = (ord_.all_ids[i].size() <= ord_.all_ids[j].size()) ? i : j;
            int P = (A == i) ? j : i;
            const auto* cl = ord_.core_ids[P].empty() ? nullptr : &ord_.core_ids[P];
            const auto& ncl = nc_f_[P];
            for (int a : ord_.all_ids[A]) {
                if (!f_.anchor_alive(a)) continue;
                add_sub(Sub{a, cl, 0, &ncl, 0});
            }
        }
    }
}

bool FrozenPairStream::head(int sub_idx, Item& out) {
    Sub& sb = subs_[sub_idx];
    while (true) {
        bool use_core = sb.core_list && sb.ci < sb.core_list->size() &&
                        f_.core_partners_enabled(sb.anchor);
        bool use_nc = sb.ni < sb.nc_list->size();
        if (!use_core && !use_nc) return false;
        int pc = use_core ? (*sb.core_list)[sb.ci] : INT_MAX;
        int pn = use_nc ? (*sb.nc_list)[sb.ni] : INT_MAX;
        bool from_core = pc < pn;
        int p = from_core ? pc : pn;
        if (f_.pair_alive(sb.anchor, p)) {
            out = Item{canon(sb.anchor, p), sub_idx, from_core};
            return true;
        }
        if (from_core) sb.ci++; else sb.ni++;
    }
}

void FrozenPairStream::push_head(int sub_idx) {
    Item it;
    if (!head(sub_idx, it)) return;
    heap_.push_back(it);
    std::push_heap(heap_.begin(), heap_.end(), std::greater<Item>());
}

bool FrozenPairStream::next(long long& score, PairKey& key) {
    while (true) {
        if (!level_active_) {
            if (!next_level()) return false;
            level_active_ = true;
            has_last_ = false;
            build_level();
        }
        while (!heap_.empty()) {
            std::pop_heap(heap_.begin(), heap_.end(), std::greater<Item>());
            Item it = heap_.back();
            heap_.pop_back();
            Sub& sb = subs_[it.sub];
            int p = (it.key.a == sb.anchor) ? it.key.b : it.key.a;

            bool core_ok = !it.from_core || f_.core_partners_enabled(sb.anchor);
            if (!core_ok) { push_head(it.sub); continue; }
            if (!f_.pair_alive(sb.anchor, p)) {
                if (it.from_core) sb.ci++; else sb.ni++;
                push_head(it.sub);
                continue;
            }

            if (it.from_core) sb.ci++; else sb.ni++;
            push_head(it.sub);
            has_last_ = true;
            last_key_ = it.key;
            score = level_score_;
            key = it.key;
            return true;
        }
        level_active_ = false;
    }
}

// -------------------------------------------------------------------------
// run_heuristic
// -------------------------------------------------------------------------
HeuristicResult run_heuristic(const Graph& G0, int s, int b, ScoreFamily fam, IncrementPolicy pol,
                              unsigned long long seed) {
    HeuristicResult res;
    res.name = heuristic_name(fam, pol);

    Graph G = G0.copy();
    const int n = G.n;
    std::vector<int> all(n);
    std::iota(all.begin(), all.end(), 0);

    Coreness c0(n, {-1, -1});
    res.initial_core = calculate_s_core(G, all, s, c0);
    std::vector<char> core0(n);
    for (int v = 0; v < n; v++) core0[v] = G.label[v] ? 1 : 0;

    std::vector<long long> q = compute_initial_node_scores(G, s, c0, fam, seed);
    FrozenPairOrder ord(q, core0);

    auto c0_of = [&](int x) { return core0[x] ? s : c0[x].first; };
    auto lb0 = [&](int x, int y) { return s - std::min(c0_of(x), c0_of(y)); };

    FrozenPairStream st(ord);
    int spent = 0;
    long long score = 0;
    PairKey key;

    auto apply = [&](int u, int v, int d) {
        if (G.has_edge(u, v)) G.increase_weight(u, v, d);
        else G.add_edge(u, v, d);
    };

    if (pol == IncrementPolicy::InitialLB) {
        auto set_filter = [&]() {
            int thr = s - (b - spent);
            OrderFilter f;
            f.anchor_alive = [&, thr](int a) { return core0[a] || c0[a].first >= thr; };
            f.nc_partner_ok = [&, thr](int p) { return c0[p].first >= thr; };
            f.core_partners_enabled = [](int) { return true; };
            f.pair_alive = [](int, int) { return true; };
            st.reset_filter(f);
        };
        if (b > 0) set_filter();
        while (spent < b && st.next(score, key)) {
            int d = lb0(key.a, key.b);
            if (d < 1 || d > b - spent) throw std::logic_error("LB: filtered pair is not affordable");
            apply(key.a, key.b, d);
            spent += d;
            res.num_actions++;
            if (spent < b) set_filter();
        }
        Coreness cf(n, {-1, -1});
        res.final_core = calculate_s_core(G, all, s, cf);
    } else {
        Coreness cur = c0;
        int cur_core = res.initial_core;
        std::vector<char> memo(n, 0);
        auto set_filter = [&]() {
            int thr = s - (b - spent);
            std::fill(memo.begin(), memo.end(), 0);
            OrderFilter f;
            f.anchor_alive = [&, thr](int a) { return G.label[a] || cur[a].first >= thr; };
            f.nc_partner_ok = [&, thr](int p) { return G.label[p] || cur[p].first >= thr; };
            f.core_partners_enabled = [&](int a) { return !G.label[a] && !memo[a]; };
            f.pair_alive = [&](int a, int p) {
                bool la = G.label[a], lp = G.label[p];
                if (la && lp) return false;
                if (la && memo[p]) return false;
                if (lp && memo[a]) return false;
                return true;
            };
            st.reset_filter(f);
        };
        if (b > 0) set_filter();
        while (spent < b && st.next(score, key)) {
            int r = b - spent;
            CandidateEval ce = evaluate_candidate(G, s, r, key.a, key.b, "mincost", cur);
            if (!ce.feasible) {
                if (G.label[key.a]) memo[key.b] = 1;
                else if (G.label[key.b]) memo[key.a] = 1;
                continue;
            }
            if (ce.delta < ce.lower_bound || ce.delta > r || ce.gain <= 0) {
                throw std::logic_error("MEI: inconsistent increment");
            }
            apply(key.a, key.b, ce.delta);
            spent += ce.delta;
            cur.assign(n, {-1, -1});
            int now = calculate_s_core(G, all, s, cur);
            if (now - cur_core != ce.gain) {
                throw std::logic_error("MEI: follower count differs from the s-core update");
            }
            cur_core = now;
            res.num_actions++;
            if (spent < b) set_filter();
        }
        res.final_core = cur_core;
    }

    res.spent = spent;
    res.gain = res.final_core - res.initial_core;
    return res;
}
