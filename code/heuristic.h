#pragma once
#include "graph.h"
#include "functions.h"
#include <functional>
#include <unordered_set>
#include <string>
#include <utility>
#include <vector>

enum class ScoreFamily { Deg, Weight, DegL, WeightL, Random };
enum class IncrementPolicy { FullInitialLB, CostCurrentMEI };

bool parse_heuristic_name(const std::string& name, ScoreFamily& fam, IncrementPolicy& pol);
std::string heuristic_name(ScoreFamily fam, IncrementPolicy pol);

std::vector<long long> compute_initial_node_scores(const Graph& G0, int s, const Coreness& c0,
                                                   ScoreFamily fam, unsigned long long seed = 0);

struct PairKey {
    int a = -1, b = -1;
    bool operator<(const PairKey& o) const { return a != o.a ? a < o.a : b < o.b; }
    bool operator==(const PairKey& o) const { return a == o.a && b == o.b; }
    bool operator<=(const PairKey& o) const { return !(o < *this); }
};
inline PairKey canon(int x, int y) { return x < y ? PairKey{x, y} : PairKey{y, x}; }

struct FrozenPairOrder {
    std::vector<long long> vals;
    std::vector<std::vector<int>> all_ids;
    std::vector<std::vector<int>> core_ids;
    std::vector<std::vector<int>> nc_ids;
    std::vector<int> bucket_of;
    std::vector<char> core0;

    FrozenPairOrder(const std::vector<long long>& qprime, const std::vector<char>& core0);

    int num_buckets() const { return (int)vals.size(); }
    bool bucket_pair_nonempty(int i, int j) const;
    long long bucket_pair_count(int i, int j) const;
};

struct OrderFilter {
    std::function<bool(int)> anchor_alive;
    std::function<bool(int)> nc_partner_ok;
    std::function<bool(int)> core_partners_enabled;
    std::function<bool(int, int)> pair_alive;
};

class FrozenPairStream {
public:
    explicit FrozenPairStream(const FrozenPairOrder& ord);
    void reset_filter(const OrderFilter& f);
    bool next(long long& score, PairKey& key);

private:
    struct Sub {
        int anchor;
        const std::vector<int>* core_list;
        size_t ci;
        const std::vector<int>* nc_list;
        size_t ni;
    };
    struct Item {
        PairKey key;
        int sub;
        bool from_core;
        bool operator>(const Item& o) const { return o.key < key; }
    };

    const FrozenPairOrder& ord_;
    OrderFilter f_;
    std::vector<std::vector<int>> nc_f_;

    struct Row { long long sum; int ii, jj; bool operator<(const Row& o) const { return sum < o.sum; } };
    std::vector<int> live_;
    std::vector<Row> gen_heap_;
    bool started_ = false;
    bool level_active_ = false;
    long long level_score_ = 0;
    std::vector<std::pair<int,int>> level_pairs_;
    bool has_last_ = false;
    PairKey last_key_;
    std::vector<Sub> subs_;
    std::vector<Item> heap_;

    bool next_level();
    void rebuild_gen(bool capped, long long cap);
    void build_level();
    bool head(int sub_idx, Item& out);
    void push_head(int sub_idx);
};

struct HeuristicResult {
    std::string name;
    int initial_core = 0;
    int final_core = 0;
    int gain = 0;
    int spent = 0;
    int num_actions = 0;
};

HeuristicResult run_heuristic(const Graph& G0, int s, int b, ScoreFamily fam, IncrementPolicy pol,
                              unsigned long long seed = 0);
