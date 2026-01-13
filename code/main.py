import argparse
import networkx as nx
import time
import psutil
import os
import math
import csv
import sys

import experiment, experiment_iter, exp_func, exp_func_iter 
import exact
import compare

parser = argparse.ArgumentParser()
parser.add_argument('--s', type=int, default=15,
                    help='user parameter s')
parser.add_argument('--b', type=int, default=30,
                    help='budget b')
parser.add_argument('--t', type=str, default='',
                    help='a folder name for bound function') # ignore this argument
parser.add_argument('--algorithm', default="exp",
                    help='specify algorithm name')
parser.add_argument('--network', default="../dataset/test/new_network.dat",
                    help='a folder name containing network.dat')
parser.add_argument('--tactics', default="TTT",
                    help='ON/OFF for Tactics')
parser.add_argument('--calculating_iter', default="F",
                    help='calculating iteration of calculating followers or not')
parser.add_argument('--compare_tactic', default="random",
                    help='How to select anchor edge')
parser.add_argument('--delta_tactic', default="compute",
                    help='How to calculate delta')
parser.add_argument('--output_path', default="../output/results.csv",
                    help='csv file for saving results')
args = parser.parse_args()

process = psutil.Process(os.getpid())
memory_before = process.memory_info().rss / (1024 * 1024)  # Convert to MB

G = nx.read_weighted_edgelist(f'{args.network}', nodetype=int)

original_nodes = sorted(G.nodes())
relabel_map = {old_label: new_label for new_label, old_label in enumerate(original_nodes)}

G = nx.relabel_nodes(G, relabel_map)
for u, v, data in G.edges(data=True):
    data["weight"] = int(math.ceil(data["weight"]))

print(f"data: {args.network.split('/')[2]}, nodes: {len(G.nodes())}, edges: {len(G.edges())}")

if args.algorithm == "exp":
    if args.tactics[0] == 'T':
        T1_self_edge = True
    else:
        T1_self_edge = False
    if args.tactics[1] == 'T':
        T2_upperbound = True
    else:
        T2_upperbound = False
    if args.tactics[2] == 'T':
        T3_reuse = True
    else:
        T3_reuse = False
    
    if args.calculating_iter == "T":
        start_time = time.time()
        A, FT, UT, G_prime, total_follower, num_iteration = experiment_iter.run(G, args.s, args.b, args.t, T1_self_edge, T2_upperbound, T3_reuse, args.delta_tactic)
        end_time = time.time()
        
    else:
        start_time = time.time()
        A, FT, UT, G_prime, total_follower = experiment.run(G, args.s, args.b, args.t, T1_self_edge, T2_upperbound, T3_reuse, args.delta_tactic)
        end_time = time.time()
    total_time = end_time - start_time
    s_core_num = 0
    for i in G_prime.nodes:
        if G_prime.nodes[i]['label']:
            s_core_num += 1

    memory_after = process.memory_info().rss / (1024 * 1024)  # Convert to MB
    memory_usage = memory_after - memory_before  # Calculate memory used

    
    if args.calculating_iter == "T":
        print(A)
        print(len(A), s_core_num, total_follower, num_iteration)
        exp_func_iter.save_result_to_csv(A, s_core_num, total_time, memory_usage, args, total_follower, num_iteration)
    else:
        print(A)
        print(len(A), s_core_num, total_follower)
        exp_func.save_result_to_csv(A, s_core_num, total_time, memory_usage, args, total_follower)
    

elif args.algorithm == "exact":
    start_time = time.time()
    A, new_score_size, follower = exact.run(G, args.s, args.b, args.t)
    end_time = time.time()
    total_time = end_time - start_time
    
    write_header = not os.path.exists(args.output_path)
    data = args.network.split('/')[4]

    with open(args.output_path, 'a', newline='') as f:
        writer = csv.writer(f)
        if write_header:
            writer.writerow(["data", "s", "b", "num_s_core", "follower", "num_anchored_edges", "total_time"])
        writer.writerow(
                [data, args.s, args.b, new_score_size, follower, len(A), total_time])
    print(f"Saved results to {data, args.s, args.b}")

elif args.algorithm == "compare":
    start_time = time.time()
    budget, new_score_size, total_follower = compare.run(G, args.s, args.b, args.t, args.compare_tactic, args.delta_tactic)
    end_time = time.time()
    total_time = end_time - start_time
    
    write_header = not os.path.exists(args.output_path)
    data = args.network.split('/')[4]

    with open(args.output_path, 'a', newline='') as f:
        writer = csv.writer(f)
        if write_header:
            writer.writerow(["data", "s", "b", "compare", "num_s_core", "follower", "num_anchored_edges", "total_time"])
        writer.writerow(
                [data, args.s, args.b, args.compare_tactic, new_score_size, total_follower, budget, total_time])
    print(f"Saved results to {data, args.s, args.b}")


cmd = ' '.join(sys.argv)  
log_path = "../output/cmdlog.txt"


with open(log_path, 'a') as f:
    f.write(cmd + '\n')