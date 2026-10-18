#pragma once
#include "graph.h"
#include "functions.h"
#include "exp_func.h"
#include <vector>
#include <string>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <utility>

struct AnchorRecord {
    std::pair<int,int> edge;
    int delta;
    double FR;
    int num_followers;
    int lower_bound;
    bool self_flag;
    bool reinforce_flag;
};
using AnchorRecordIter = AnchorRecord;

struct RunResult {
    std::vector<AnchorRecord> A;
    double FT;
    double UT;
    Graph G_prime;
    int total_follower;
    int initial_s_core = 0;
    int spent = 0;
};

struct IterRunResult {
    std::vector<AnchorRecordIter> A;
    double FT;
    double UT;
    Graph G_prime;
    int total_follower;
    long long global_cnt;
    int initial_s_core = 0;
    int spent = 0;
};

struct CCCache {
    std::unordered_map<int,int> comp_of;
    std::unordered_map<int, std::unordered_set<int>> nodes_in;
    std::map<int, BestEdge> intra_best;
    int next_cid = 0;
};

CCCache build_initial_caches(Graph& G, int s, int b, int spent,
                               Coreness& coreness,
                               std::unordered_map<int,double>& upperbound,
                               bool T1_self_edge, bool T2_upperbound,
                               int& s_cand,
                               const std::string& delta_tactic,
                               long long* counter = nullptr);

BestEdge find_intra_best(Graph& G, const std::unordered_set<int>& nodes,
                          Coreness& coreness, int s, int b, int spent,
                          std::unordered_map<int,double>& upperbound,
                          bool T1_self_edge, bool T2_upperbound,
                          int s_cand, const std::string& delta_tactic,
                          long long* counter = nullptr);

RunResult run(Graph& G, int s, int b,
              bool T1_self_edge, bool T2_upperbound, bool T3_reuse,
              const std::string& delta_tactic);

IterRunResult run_iter(Graph& G, int s, int b,
                        bool T1_self_edge, bool T2_upperbound, bool T3_reuse,
                        const std::string& delta_tactic);
