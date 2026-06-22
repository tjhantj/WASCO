#pragma once
#include "graph.h"
#include "functions.h"
#include <vector>
#include <unordered_map>
#include <string>
#include <utility>
#include <tuple>
#include <optional>

// Best edge result tuple
struct BestEdge {
    std::pair<int,int> edge;  // (u, v) canonical (u <= v)
    int delta;
    double FR;
    int num_followers;
    bool valid;  // false = no edge found

    BestEdge() : edge({-1,-1}), delta(0), FR(0.0), num_followers(0), valid(false) {}
    BestEdge(std::pair<int,int> e, int d, double fr, int nf)
        : edge(e), delta(d), FR(fr), num_followers(nf), valid(true) {}
};

// Self-edge pruning: returns non_s_core list and s_cand (min-id s-core node)
std::pair<std::vector<int>, int> self_edge_pruning(const Graph& G);

// Candidate node lists (T1 tactic: only non-s-core)
std::vector<int> make_candidate_nodes(Graph& G, const std::vector<int>& nodes, int s,
                                       int budget_left, const Coreness& coreness,
                                       std::unordered_map<int,double>& upperbound,
                                       bool T2_upperbound,
                                       const std::string& delta_tactic);

// Candidate node list (for non-T1 with T2: includes s-core nodes with ub=0)
std::vector<int> make_candidate_nodes_v2(Graph& G, int s, int budget_left,
                                          const Coreness& coreness,
                                          std::unordered_map<int,double>& upperbound,
                                          bool T2_upperbound,
                                          const std::string& delta_tactic);

// Candidate edges (non-T1, non-T2)
std::vector<std::pair<int,int>> make_candidate_edges(const Graph& G, int s, int budget_left,
                                                      const Coreness& coreness);

// Iteration with upperbound pruning
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

// Iteration without upperbound (T1 only, no T2)
BestEdge iteration_nodes_no_upperbound(Graph& G, const std::vector<int>& candidate_nodes,
                                        int s, int b, int spent,
                                        Coreness& coreness,
                                        int s_cand,
                                        bool T3_reuse,
                                        const std::string& delta_tactic,
                                        const std::unordered_map<int,int>& comp_of,
                                        const BestEdge& initial_best,
                                        long long* counter = nullptr);

// Iteration over explicit edge list (no T1, no T2)
BestEdge iteration_edges_no_upperbound(Graph& G,
                                        const std::vector<std::pair<int,int>>& candidate_edges,
                                        int s, int b, int spent,
                                        Coreness& coreness,
                                        bool T3_reuse,
                                        const std::string& delta_tactic,
                                        const std::unordered_map<int,int>& comp_of,
                                        const BestEdge& initial_best,
                                        long long* counter = nullptr);

// Compute delta_e based on delta_tactic
// orig_u, orig_v: original (u,v) before self-edge transform, used for has_edge check in naive/percentage
int compute_delta_e(Graph& G, int s, int u, int v, const std::string& delta_tactic, const Coreness& coreness,
                    int orig_u = -1, int orig_v = -1);

// Update best if FR is better (or equal with lexicographically smaller edge)
bool is_better(const BestEdge& candidate, const BestEdge& current);
