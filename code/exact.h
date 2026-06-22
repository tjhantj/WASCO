#pragma once
#include "graph.h"
#include "functions.h"
#include <vector>
#include <utility>
#include <unordered_map>
#include <unordered_set>

struct ExactRecord {
    std::pair<int,int> edge;
    int delta;
};

struct ExactResult {
    std::vector<ExactRecord> A;
    int final_size;
    int followers;
};

ExactResult exact_run(Graph& G, int s, int b);
