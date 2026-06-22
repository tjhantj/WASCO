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
};

struct AnchorRecordIter {
    std::pair<int,int> edge;
    int delta;
    double FR;
    int num_followers;
    bool self_flag;
    bool reinforce_flag;
};

struct RunResult {
    std::vector<AnchorRecord> A;
    double FT;
    double UT;
    Graph G_prime;
    int total_follower;
};

struct IterRunResult {
    std::vector<AnchorRecordIter> A;
    double FT;
    double UT;
    Graph G_prime;
    int total_follower;
    long long global_cnt;
};

// CC cache structures
struct CCCache {
    std::unordered_map<int,int> comp_of;           // node -> CC id
    std::unordered_map<int, std::unordered_set<int>> nodes_in;  // CC id -> node set
    std::map<int, BestEdge> intra_best;             // CC id -> best edge (ordered for Python-matching iteration)
};

// Build initial CC caches
CCCache build_initial_caches(Graph& G, int s, int b, int spent,
                               Coreness& coreness,
                               std::unordered_map<int,double>& upperbound,
                               bool T1_self_edge, bool T2_upperbound,
                               int& s_cand,
                               const std::string& delta_tactic,
                               long long* counter = nullptr);

// Find best intra-CC edge
BestEdge find_intra_best(Graph& G, const std::unordered_set<int>& nodes,
                          Coreness& coreness, int s, int b, int spent,
                          std::unordered_map<int,double>& upperbound,
                          bool T1_self_edge, bool T2_upperbound,
                          int s_cand, const std::string& delta_tactic,
                          long long* counter = nullptr);

// Main GPA run function
RunResult run(Graph& G, int s, int b,
              bool T1_self_edge, bool T2_upperbound, bool T3_reuse,
              const std::string& delta_tactic);

// Iter version (tracks GLOBAL_CNT)
IterRunResult run_iter(Graph& G, int s, int b,
                        bool T1_self_edge, bool T2_upperbound, bool T3_reuse,
                        const std::string& delta_tactic);
