#pragma once
#include "graph.h"
#include "functions.h"
#include <vector>
#include <utility>
#include <string>

struct CompareRecord {
    std::pair<int,int> edge;
    int delta;
};

struct CompareResult {
    int budget_used;
    int final_size;
    int total_followers;
};

CompareResult compare_run(Graph& G, int s, int b,
                           const std::string& compare_tactic,
                           const std::string& delta_tactic);
