#pragma once
#include "graph.h"
#include "functions.h"
#include <vector>
#include <unordered_map>
#include <string>
#include <utility>
#include <tuple>

struct BestEdge {
    std::pair<int,int> edge;
    int delta;
    double FR;
    int num_followers;
    bool valid;
    int lower_bound = 0;
    bool core_rep = false;
    int region = -1;

    BestEdge() : edge({-1,-1}), delta(0), FR(0.0), num_followers(0), valid(false) {}
    BestEdge(std::pair<int,int> e, int d, double fr, int nf)
        : edge(e), delta(d), FR(fr), num_followers(nf), valid(true) {}
};

struct CandidateEval {
    bool feasible = false;
    int lower_bound = 0;
    int probe_delta = 0;
    int delta = 0;
    int gain = 0;
};

CandidateEval evaluate_candidate(Graph& G, int s, int budget_left, int eu, int ev,
                                 const std::string& delta_tactic, Coreness& coreness,
                                 int orig_u = -1, int orig_v = -1,
                                 long long* counter = nullptr);

std::pair<std::vector<int>, int> self_edge_pruning(const Graph& G);

std::vector<int> make_candidate_nodes(Graph& G, const std::vector<int>& nodes, int s,
                                       int budget_left, const Coreness& coreness,
                                       std::unordered_map<int,double>& upperbound,
                                       bool T2_upperbound,
                                       const std::string& delta_tactic);

std::vector<int> make_candidate_nodes_v2(Graph& G, int s, int budget_left,
                                          const Coreness& coreness,
                                          std::unordered_map<int,double>& upperbound,
                                          bool T2_upperbound,
                                          const std::string& delta_tactic);

std::vector<std::pair<int,int>> make_candidate_edges(const Graph& G, int s, int budget_left,
                                                      const Coreness& coreness);

BestEdge iteration_nodes_upperbound(Graph& G, const std::vector<int>& candidate_nodes,
                                     int s, int b, int spent,
                                     Coreness& coreness,
                                     std::unordered_map<int,double>& upperbound,
                                     int s_cand,
                                     bool T1_self_edge, bool T3_reuse,
                                     const std::string& delta_tactic,
                                     const std::unordered_map<int,int>& comp_of,
                                     const BestEdge& initial_best,
                                     long long* counter = nullptr);

BestEdge iteration_nodes_no_upperbound(Graph& G, const std::vector<int>& candidate_nodes,
                                        int s, int b, int spent,
                                        Coreness& coreness,
                                        int s_cand,
                                        bool T3_reuse,
                                        const std::string& delta_tactic,
                                        const std::unordered_map<int,int>& comp_of,
                                        const BestEdge& initial_best,
                                        long long* counter = nullptr);

BestEdge iteration_edges_no_upperbound(Graph& G,
                                        const std::vector<std::pair<int,int>>& candidate_edges,
                                        int s, int b, int spent,
                                        Coreness& coreness,
                                        bool T3_reuse,
                                        const std::string& delta_tactic,
                                        const std::unordered_map<int,int>& comp_of,
                                        const BestEdge& initial_best,
                                        long long* counter = nullptr);

int compute_delta_e(Graph& G, int s, int u, int v, const std::string& delta_tactic, const Coreness& coreness,
                    int orig_u = -1, int orig_v = -1);

bool is_better(const BestEdge& candidate, const BestEdge& current);
