#pragma once
#include "graph.h"
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <utility>
#include <tuple>

// coreness[v] = {c, l} where c=core level, l=layer; {-1,-1} means not computed (node is s-core)
using Coreness = std::vector<std::pair<int,int>>;

// Calculate s-core decomposition.
// Sets G.label[v]=true for nodes remaining in s-core, false for peeled ones.
// Updates coreness for peeled nodes.
// Returns count of s-core nodes.
int calculate_s_core(Graph& G, const std::vector<int>& nodes, int s, Coreness& coreness);

// Compute minimum edge delta needed to reach s-core
int computeDelta(Graph& G, int s, int u, int v, const Coreness& coreness);

// Find followers: nodes that would join s-core if edge (u,v) gets weight += delta_e
// Temporarily modifies G, then rolls back. Returns set of follower nodes.
std::unordered_set<int> FindFollowers(int u, int v, int delta_e, Graph& G, int s, Coreness& coreness);

// Upper bound on follower ratio for node u
double Upperbound(Graph& G, int u, const Coreness& coreness, int s, const std::string& delta_tactic);

// U_single: just return upperbound[u]
double U_single(int u, const std::unordered_map<int,double>& upperbound);

// U_double: combined upper bound for pair (u,v)
double U_double(int u, int v, const std::unordered_map<int,double>& upperbound,
                const Coreness& coreness, const Graph& G, int s);
