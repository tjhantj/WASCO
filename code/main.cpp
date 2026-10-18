#include "graph.h"
#include "functions.h"
#include "exp_func.h"
#include "experiment.h"
#include "exact.h"
#include "heuristic.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <chrono>
#include <cstring>
#include <algorithm>
#include <stdexcept>
#include <sys/stat.h>

// -------------------------------------------------------------------------
// Memory usage (MB)
// -------------------------------------------------------------------------
static double get_memory_mb() {
    std::ifstream f("/proc/self/status");
    if (!f.is_open()) return 0.0;
    std::string line;
    while (std::getline(f, line)) {
        if (line.rfind("VmRSS:", 0) == 0) {
            std::istringstream ss(line.substr(6));
            double kb;
            ss >> kb;
            return kb / 1024.0;
        }
    }
    return 0.0;
}

// -------------------------------------------------------------------------
// CSV helper
// -------------------------------------------------------------------------
static bool file_exists(const std::string& path) {
    std::ifstream f(path);
    return f.good();
}

static bool file_nonempty(const std::string& path) {
    std::ifstream f(path);
    return f.good() && f.peek() != std::ifstream::traits_type::eof();
}

static std::string join_cols(const std::vector<std::string>& cols) {
    std::string s;
    for (int i = 0; i < (int)cols.size(); i++) {
        if (i) s += ",";
        s += cols[i];
    }
    return s;
}

static bool csv_header_ok(const std::string& path, const std::vector<std::string>& cols) {
    std::ifstream f(path);
    if (!f.good()) return true;
    std::string first;
    if (!std::getline(f, first) || first.empty()) return true;
    if (!first.empty() && first.back() == '\r') first.pop_back();
    if (first == join_cols(cols)) return true;
    std::cerr << "Error: " << path << " has a different CSV header.\n"
              << "  found:    " << first << "\n"
              << "  expected: " << join_cols(cols) << "\n"
              << "Use a new --output_path." << std::endl;
    return false;
}

static void write_csv_header(const std::string& path, const std::vector<std::string>& cols) {
    std::ofstream f(path, std::ios::app);
    f << join_cols(cols) << "\n";
}

static void write_csv_row(const std::string& path, const std::vector<std::string>& vals) {
    std::ofstream f(path, std::ios::app);
    f << join_cols(vals) << "\n";
}

// -------------------------------------------------------------------------
// Argument parsing
// -------------------------------------------------------------------------
struct Args {
    std::string algorithm = "exp";
    std::string network = "../dataset/real/exact/karate/network.dat";
    int s = 15;
    int b = 60;
    std::string tactics = "TTT";
    std::string calculating_iter = "F";
    std::string delta_tactic = "mincost";
    std::string heuristic = "Deg-MEI";
    int seed = 1;
    std::string output_path = "../output/results.csv";
};

static Args parse_args(int argc, char* argv[]) {
    Args args;
    for (int i = 1; i < argc; i++) {
        std::string key = argv[i];
        if (i + 1 >= argc) continue;
        std::string val = argv[i + 1];
        if (key == "--algorithm") { args.algorithm = val; i++; }
        else if (key == "--network") { args.network = val; i++; }
        else if (key == "--s") { args.s = std::stoi(val); i++; }
        else if (key == "--b") { args.b = std::stoi(val); i++; }
        else if (key == "--tactics") { args.tactics = val; i++; }
        else if (key == "--calculating_iter") { args.calculating_iter = val; i++; }
        else if (key == "--delta_tactic") { args.delta_tactic = val; i++; }
        else if (key == "--heuristic") { args.heuristic = val; i++; }
        else if (key == "--seed") { args.seed = std::stoi(val); i++; }
        else if (key == "--output_path") { args.output_path = val; i++; }
    }
    return args;
}

static std::string extract_data_name(const std::string& path) {
    std::vector<std::string> parts;
    std::stringstream ss(path);
    std::string token;
    while (std::getline(ss, token, '/')) {
        parts.push_back(token);
    }
    if (parts.size() > 4) return parts[4];
    if (!parts.empty()) return parts.back();
    return path;
}

// -------------------------------------------------------------------------
// main
// -------------------------------------------------------------------------
int main(int argc, char* argv[]) {
    Args args = parse_args(argc, argv);

    if (args.delta_tactic != "mincost" && args.delta_tactic != "compute" &&
        args.delta_tactic != "naive") {
        std::cerr << "Unknown --delta_tactic: " << args.delta_tactic
                  << " (use mincost, compute or naive)" << std::endl;
        return 1;
    }

    double memory_before = get_memory_mb();

    Graph G_orig;
    try {
        G_orig = load_graph(args.network);
    } catch (const std::exception& e) {
        std::cerr << "Error loading graph: " << e.what() << std::endl;
        return 1;
    }

    std::cout << "data: " << extract_data_name(args.network)
              << ", nodes: " << G_orig.n << std::endl;

    std::string data = extract_data_name(args.network);

    bool T1 = (args.tactics.size() > 0 && args.tactics[0] == 'T');
    bool T2 = (args.tactics.size() > 1 && args.tactics[1] == 'T');
    bool T3 = (args.tactics.size() > 2 && args.tactics[2] == 'T');

    {
        std::string dir = args.output_path;
        auto pos = dir.rfind('/');
        if (pos != std::string::npos) {
            dir = dir.substr(0, pos);
            std::string cmd = "mkdir -p \"" + dir + "\"";
            system(cmd.c_str());
        }
    }

    auto t_start = std::chrono::high_resolution_clock::now();

    if (args.algorithm == "exp") {
        const std::vector<std::string> iter_cols =
            {"data","s","b","tactics","num_s_core","follower","num_iter",
             "num_anchored_edges","iter","delta","num_follower","self","reinforce",
             "total_time","memory","delta_tactic","edge_u","edge_v","lower_bound"};
        const std::vector<std::string> std_cols =
            {"data","s","b","tactics","num_s_core","follower","num_anchored_edges",
             "iter","delta","num_follower","total_time","memory",
             "delta_tactic","edge_u","edge_v","lower_bound"};
        const bool iter_mode = (args.calculating_iter == "T");
        const auto& cols = iter_mode ? iter_cols : std_cols;
        if (!csv_header_ok(args.output_path, cols)) return 1;

        Graph G_copy = G_orig.copy();
        IterRunResult res;
        if (iter_mode) {
            res = run_iter(G_copy, args.s, args.b, T1, T2, T3, args.delta_tactic);
        } else {
            RunResult r = run(G_copy, args.s, args.b, T1, T2, T3, args.delta_tactic);
            res.A = std::move(r.A);
            res.FT = r.FT;
            res.UT = r.UT;
            res.G_prime = std::move(r.G_prime);
            res.total_follower = r.total_follower;
            res.global_cnt = 0;
            res.initial_s_core = r.initial_s_core;
            res.spent = r.spent;
        }
        auto t_end = std::chrono::high_resolution_clock::now();
        double total_time = std::chrono::duration<double>(t_end - t_start).count();

        double memory_after = get_memory_mb();
        double memory_usage = memory_after - memory_before;

        int s_core_num = 0;
        for (int i = 0; i < res.G_prime.n; i++) {
            if (res.G_prime.label[i]) s_core_num++;
        }

        std::cout << "delta_tactic=" << args.delta_tactic
                  << " num_anchored=" << res.A.size()
                  << " initial_s_core=" << res.initial_s_core
                  << " s_core=" << s_core_num
                  << " total_follower=" << res.total_follower
                  << " spent=" << res.spent << "/" << args.b;
        if (iter_mode) std::cout << " global_cnt=" << res.global_cnt;
        std::cout << std::endl;
        if (s_core_num - res.initial_s_core != res.total_follower) {
            std::cerr << "Warning: total_follower (" << res.total_follower
                      << ") != final s-core - initial s-core ("
                      << (s_core_num - res.initial_s_core) << ")" << std::endl;
        }

        bool write_header = !file_nonempty(args.output_path);
        if (write_header && !res.A.empty()) {
            write_csv_header(args.output_path, cols);
        }
        for (int idx = 0; idx < (int)res.A.size(); idx++) {
            auto& rec = res.A[idx];
            std::vector<std::string> row;
            if (iter_mode) {
                row = {
                    data,
                    std::to_string(args.s),
                    std::to_string(args.b),
                    args.tactics,
                    std::to_string(s_core_num),
                    std::to_string(res.total_follower),
                    std::to_string(res.global_cnt),
                    std::to_string((int)res.A.size()),
                    std::to_string(idx+1),
                    std::to_string(rec.delta),
                    std::to_string(rec.num_followers),
                    rec.self_flag ? "True" : "False",
                    rec.reinforce_flag ? "True" : "False",
                    std::to_string(total_time),
                    std::to_string(memory_usage)
                };
            } else {
                row = {
                    data,
                    std::to_string(args.s),
                    std::to_string(args.b),
                    args.tactics,
                    std::to_string(s_core_num),
                    std::to_string(res.total_follower),
                    std::to_string((int)res.A.size()),
                    std::to_string(idx+1),
                    std::to_string(rec.delta),
                    std::to_string(rec.num_followers),
                    std::to_string(total_time),
                    std::to_string(memory_usage)
                };
            }
            row.push_back(args.delta_tactic);
            row.push_back(std::to_string(rec.edge.first));
            row.push_back(std::to_string(rec.edge.second));
            row.push_back(std::to_string(rec.lower_bound));
            write_csv_row(args.output_path, row);
        }
        std::cout << "Saved results to (" << data << ", " << args.s << ", " << args.b << ", " << args.tactics << ")" << std::endl;

    } else if (args.algorithm == "exact") {
        Graph G_copy = G_orig.copy();
        ExactResult res = exact_run(G_copy, args.s, args.b);
        auto t_end = std::chrono::high_resolution_clock::now();
        double total_time = std::chrono::duration<double>(t_end - t_start).count();

        bool write_header = !file_exists(args.output_path);
        if (write_header) {
            write_csv_header(args.output_path,
                {"data","s","b","num_s_core","follower","num_anchored_edges","total_time"});
        }
        write_csv_row(args.output_path, {
            data,
            std::to_string(args.s),
            std::to_string(args.b),
            std::to_string(res.final_size),
            std::to_string(res.followers),
            std::to_string((int)res.A.size()),
            std::to_string(total_time)
        });
        std::cout << "Saved results to (" << data << ", " << args.s << ", " << args.b << ")" << std::endl;

    } else if (args.algorithm == "heuristic") {
        ScoreFamily fam;
        IncrementPolicy pol;
        if (!parse_heuristic_name(args.heuristic, fam, pol)) {
            std::cerr << "Unknown --heuristic: " << args.heuristic << std::endl;
            return 1;
        }
        const std::vector<std::string> cols =
            {"data","s","b","heuristic","seed","num_s_core","follower",
             "num_anchored_edges","budget_used","total_time"};
        if (!csv_header_ok(args.output_path, cols)) return 1;

        const bool is_random = (fam == ScoreFamily::Random);
        unsigned long long seed = is_random ? (unsigned long long)args.seed : 0;
        HeuristicResult res = run_heuristic(G_orig, args.s, args.b, fam, pol, seed);
        auto t_end = std::chrono::high_resolution_clock::now();
        double total_time = std::chrono::duration<double>(t_end - t_start).count();

        std::cout << res.name;
        if (is_random) std::cout << " seed=" << args.seed;
        std::cout << " initial_s_core=" << res.initial_core
                  << " s_core=" << res.final_core
                  << " follower=" << res.gain
                  << " spent=" << res.spent << "/" << args.b << std::endl;

        if (!file_nonempty(args.output_path)) {
            write_csv_header(args.output_path, cols);
        }
        write_csv_row(args.output_path, {
            data,
            std::to_string(args.s),
            std::to_string(args.b),
            res.name,
            is_random ? std::to_string(args.seed) : "",
            std::to_string(res.final_core),
            std::to_string(res.gain),
            std::to_string(res.num_actions),
            std::to_string(res.spent),
            std::to_string(total_time)
        });
        std::cout << "Saved results to (" << data << ", " << args.s << ", " << args.b << ", "
                  << res.name << ")" << std::endl;

    } else {
        std::cerr << "Unknown algorithm: " << args.algorithm << std::endl;
        return 1;
    }

    std::string log_path = "../output/cmdlog.txt";
    {
        std::ofstream logf(log_path, std::ios::app);
        for (int i = 0; i < argc; i++) {
            if (i) logf << " ";
            logf << argv[i];
        }
        logf << "\n";
    }

    return 0;
}
