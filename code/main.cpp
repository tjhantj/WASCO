#include "graph.h"
#include "functions.h"
#include "exp_func.h"
#include "experiment.h"
#include "exact.h"
#include "compare.h"

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
// Memory usage via /proc/self/status (VmRSS in MB)
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

static void write_csv_header(const std::string& path, const std::vector<std::string>& cols) {
    std::ofstream f(path, std::ios::app);
    for (int i = 0; i < (int)cols.size(); i++) {
        if (i) f << ",";
        f << cols[i];
    }
    f << "\n";
}

static void write_csv_row(const std::string& path, const std::vector<std::string>& vals) {
    std::ofstream f(path, std::ios::app);
    for (int i = 0; i < (int)vals.size(); i++) {
        if (i) f << ",";
        f << vals[i];
    }
    f << "\n";
}

// -------------------------------------------------------------------------
// Argument parsing
// -------------------------------------------------------------------------
struct Args {
    std::string algorithm = "exp";
    std::string network = "../dataset/test/new_network.dat";
    int s = 15;
    int b = 60;
    std::string tactics = "TTT";
    std::string calculating_iter = "F";
    std::string compare_tactic = "random";
    std::string delta_tactic = "compute";
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
        else if (key == "--compare_tactic") { args.compare_tactic = val; i++; }
        else if (key == "--delta_tactic") { args.delta_tactic = val; i++; }
        else if (key == "--output_path") { args.output_path = val; i++; }
    }
    return args;
}

// Extract data name from network path (split by '/', take index [4])
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

    double memory_before = get_memory_mb();

    // Load graph
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

    // Parse tactics
    bool T1 = (args.tactics.size() > 0 && args.tactics[0] == 'T');
    bool T2 = (args.tactics.size() > 1 && args.tactics[1] == 'T');
    bool T3 = (args.tactics.size() > 2 && args.tactics[2] == 'T');

    // Ensure output directory exists (best-effort)
    {
        std::string dir = args.output_path;
        auto pos = dir.rfind('/');
        if (pos != std::string::npos) {
            dir = dir.substr(0, pos);
            // mkdir is system-dependent, but on Linux:
            std::string cmd = "mkdir -p \"" + dir + "\"";
            system(cmd.c_str());
        }
    }

    auto t_start = std::chrono::high_resolution_clock::now();

    if (args.algorithm == "exp") {
        if (args.calculating_iter == "T") {
            Graph G_copy = G_orig.copy();
            IterRunResult res = run_iter(G_copy, args.s, args.b, T1, T2, T3, args.delta_tactic);
            auto t_end = std::chrono::high_resolution_clock::now();
            double total_time = std::chrono::duration<double>(t_end - t_start).count();

            double memory_after = get_memory_mb();
            double memory_usage = memory_after - memory_before;

            int s_core_num = 0;
            for (int i = 0; i < res.G_prime.n; i++) {
                if (res.G_prime.label[i]) s_core_num++;
            }

            std::cout << "num_anchored=" << res.A.size()
                      << " s_core=" << s_core_num
                      << " total_follower=" << res.total_follower
                      << " global_cnt=" << res.global_cnt << std::endl;

            // CSV output (iter version)
            bool write_header = !file_exists(args.output_path);
            if (write_header) {
                write_csv_header(args.output_path,
                    {"data","s","b","tactics","num_s_core","follower","num_iter",
                     "num_anchored_edges","iter","delta","num_follower","self","reinforce",
                     "total_time","memory"});
            }
            for (int idx = 0; idx < (int)res.A.size(); idx++) {
                auto& rec = res.A[idx];
                write_csv_row(args.output_path, {
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
                });
            }
            std::cout << "Saved results to (" << data << ", " << args.s << ", " << args.b << ", " << args.tactics << ")" << std::endl;

        } else {
            Graph G_copy = G_orig.copy();
            RunResult res = run(G_copy, args.s, args.b, T1, T2, T3, args.delta_tactic);
            auto t_end = std::chrono::high_resolution_clock::now();
            double total_time = std::chrono::duration<double>(t_end - t_start).count();

            double memory_after = get_memory_mb();
            double memory_usage = memory_after - memory_before;

            int s_core_num = 0;
            for (int i = 0; i < res.G_prime.n; i++) {
                if (res.G_prime.label[i]) s_core_num++;
            }

            std::cout << "num_anchored=" << res.A.size()
                      << " s_core=" << s_core_num
                      << " total_follower=" << res.total_follower << std::endl;

            // CSV output (standard)
            bool write_header = !file_exists(args.output_path);
            if (write_header) {
                write_csv_header(args.output_path,
                    {"data","s","b","tactics","num_s_core","follower","num_anchored_edges",
                     "iter","delta","num_follower","total_time","memory"});
            }
            for (int idx = 0; idx < (int)res.A.size(); idx++) {
                auto& rec = res.A[idx];
                write_csv_row(args.output_path, {
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
                });
            }
            std::cout << "Saved results to (" << data << ", " << args.s << ", " << args.b << ", " << args.tactics << ")" << std::endl;
        }

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

    } else if (args.algorithm == "compare") {
        Graph G_copy = G_orig.copy();
        CompareResult res = compare_run(G_copy, args.s, args.b, args.compare_tactic, args.delta_tactic);
        auto t_end = std::chrono::high_resolution_clock::now();
        double total_time = std::chrono::duration<double>(t_end - t_start).count();

        bool write_header = !file_exists(args.output_path);
        if (write_header) {
            write_csv_header(args.output_path,
                {"data","s","b","compare","num_s_core","follower","num_anchored_edges","total_time"});
        }
        write_csv_row(args.output_path, {
            data,
            std::to_string(args.s),
            std::to_string(args.b),
            args.compare_tactic,
            std::to_string(res.final_size),
            std::to_string(res.total_followers),
            std::to_string(res.budget_used),
            std::to_string(total_time)
        });
        std::cout << "Saved results to (" << data << ", " << args.s << ", " << args.b << ")" << std::endl;

    } else {
        std::cerr << "Unknown algorithm: " << args.algorithm << std::endl;
        return 1;
    }

    // Write cmdlog
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
