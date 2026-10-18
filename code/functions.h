#pragma once
#include "graph.h"
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <utility>
#include <tuple>

// coreness[v] = {core level, layer}
using Coreness = std::vector<std::pair<int,int>>;

int calculate_s_core(Graph& G, const std::vector<int>& nodes, int s, Coreness& coreness);

int computeDelta(Graph& G, int s, int u, int v, const Coreness& coreness);

std::unordered_set<int> FindFollowers(int u, int v, int delta_e, Graph& G, int s, Coreness& coreness);

double Upperbound(Graph& G, int u, const Coreness& coreness, int s, const std::string& delta_tactic);

double U_single(int u, const std::unordered_map<int,double>& upperbound);

double U_double(int u, int v, const std::unordered_map<int,double>& upperbound,
                const Coreness& coreness, const Graph& G, int s);
