import functions
import time
import csv
import os
import sys
import math

import experiment_iter

def self_edge_pruning(G):
    non_s_core = []
    s_cand = None
    for n, d in G.nodes(data=True):
        if not d['label']:
            non_s_core.append(n)
        else:
            if s_cand is None or (s_cand is not None and n < s_cand):
                s_cand = n
    
    if s_cand is None:
        print("No node in s-core. Change s value")
        sys.exit(1)
    
    return non_s_core, s_cand


def make_candidate_nodes(G_prime, nodes, s, budget_left, coreness, upperbound, UT, T2_upperbound, delta_Tactic):
    '''
    This function is only for self_edge tactic. S-core can be perfectly ignored.
    So just store non-s-core as candidate_nodes.
    If upperbound tactic is used, sort candidate_nodes with upperbound for pruning.
    '''
    candidate_nodes = []
    for u in nodes:
        # Pruning considering remain budget
        # print(budget_left)
        if not G_prime.nodes[u]['label'] and s - coreness[u][0] <= budget_left:
            candidate_nodes.append(u)
            # Compute upperbound
            if T2_upperbound:
                temp_start = time.time()
                upperbound[u] = functions.Upperbound(G_prime, u, coreness, s, delta_Tactic)
                temp_end = time.time()
                UT += temp_end - temp_start
    
    if T2_upperbound:
        candidate_nodes.sort(key = lambda x : -upperbound[x])
    
    return candidate_nodes

def iteration_nodes_upperbound(G_prime, candidate_nodes, s, b, t, spent, coreness, upperbound, s_cand, FT, T1_self_edge, T3_reuse, delta_Tactic, comp_of={}, best=(None, 0, 0.0, 0)):
    '''
    When using the upper-bound strategy, construct 'candidate_nodes' to apply the bound, 
    even if the self-edge strategy is not used. 
    Initialize the s-core upper bound to 0 and sort the candidates.
    '''
    
    # initial setting
    best_edge, best_delta, most_FR, most_follower = best

    c = len(candidate_nodes)
    for i in range(c):
        u = candidate_nodes[i]
        if most_FR > functions.U_single(u, upperbound) * 2:
            break
        for j in range(i if T1_self_edge else i+1, c):
            v = candidate_nodes[j]
            if T3_reuse:
                if comp_of[u] == comp_of[v]:    # We considered self edge in cache
                    continue
            
            e = (u, v)

            if most_FR > functions.U_single(u, upperbound) + functions.U_single(v, upperbound):
                break
            if most_FR > functions.U_double(u, v, upperbound, coreness, G_prime, s):
                continue
            else:
                e = (u, v)
                
                if u == v:
                    e = (u, s_cand)

                if delta_Tactic == "compute":
                    delta_e = functions.computeDelta(G_prime, s, e, t, coreness)

                elif delta_Tactic.isdigit(): 
                    ratio = int(delta_Tactic) / 100.0 
                    
                    mee = functions.computeDelta(G_prime, s, e, t, coreness)
                    
                    if G_prime.has_edge(u, v):
                        naive = s - G_prime[u][v]['weight']
                    else:
                        naive = s
                        
                    diff = naive - mee
                    delta_e = mee + math.ceil(diff * ratio)

                else:
                    if G_prime.has_edge(u, v):
                        delta_e = s - G_prime[u][v]['weight']
                        print("check", delta_e)
                    else:
                        delta_e = s

                if delta_e <= 0:
                    
                    pass
                
                if delta_e > 0 and spent + delta_e <= b:
                    # calculate the follower
                    temp_start = time.time()
                    followers = functions.FindFollowers(e, delta_e, G_prime, s, coreness)
                    experiment_iter.GLOBAL_CNT += 1
                    temp_end = time.time()
                    FT += temp_end - temp_start

                    FR = len(followers) / delta_e  # follower rate

                    # renew the maximum value
                    if e[0] > e[1]:
                        e = (e[1], e[0])
                    if FR > most_FR or (FR == most_FR and (best_edge is None or e < best_edge)):
                        best_edge = e
                        best_delta = delta_e
                        most_FR = FR
                        most_follower = len(followers)
    return best_edge, best_delta, most_FR, most_follower

def iteration_nodes_no_upperbound(G_prime, candidate_nodes, s, b, t, spent, coreness, s_cand, FT, T3_reuse, delta_Tactic, comp_of={}, best=(None, 0, 0.0, 0)):
    # initial setting
    best_edge, best_delta, most_FR, most_follower = best

    c = len(candidate_nodes)
    for i in range(c):
        u = candidate_nodes[i]
        for j in range(i, c):
            v = candidate_nodes[j]
            if T3_reuse:
                if comp_of[u] == comp_of[v]:
                    continue
            
            e = (u, v)
            if u == v:
                e = (u, s_cand)

            if delta_Tactic == "compute":
                delta_e = functions.computeDelta(G_prime, s, e, t, coreness)

            elif delta_Tactic.isdigit(): 
                ratio = int(delta_Tactic) / 100.0  
                
                mee = functions.computeDelta(G_prime, s, e, t, coreness)
                
                if G_prime.has_edge(u, v):
                    naive = s - G_prime[u][v]['weight']
                else:
                    naive = s
                    
                diff = naive - mee
                delta_e = mee + math.ceil(diff * ratio)

            else:
                if G_prime.has_edge(u, v):
                    delta_e = s - G_prime[u][v]['weight']
                    print("check", delta_e)
                else:
                    delta_e = s

            if delta_e > 0 and spent + delta_e <= b:
                # calculate the follower
                temp_start = time.time()
                followers = functions.FindFollowers(e, delta_e, G_prime, s, coreness)
                experiment_iter.GLOBAL_CNT += 1
                temp_end = time.time()
                FT += temp_end - temp_start

                FR = len(followers) / delta_e  # follower rate

                # renew the maximum value
                if e[0] > e[1]:
                    e = (e[1], e[0])
                if FR > most_FR or (FR == most_FR and (best_edge is None or e < best_edge)):
                    best_edge = e
                    best_delta = delta_e
                    most_FR = FR
                    most_follower = len(followers)
    return best_edge, best_delta, most_FR, most_follower


def make_candidate_edges(G_prime, nodes, s, budget_left, coreness):
    '''
    When self-edge and upper-bound strategies are not used, connections to the s-core must be considered.
    Thus, 'candidate_edges' include both edges between non-s-core nodes and edges connecting non-s-core to s-core nodes.
    '''
    candidate_edges = []

    non_s_core = []
    s_core = []
    for u in nodes:
        if s - coreness.get(u, (s, 0))[0] > budget_left:
            continue
        if not G_prime.nodes[u]['label']:
            non_s_core.append(u)
        else:
            s_core.append(u)

    # Intra non-core
    non_len = len(non_s_core)
    for i in range(non_len):
        u = non_s_core[i]
        for j in range(i+1, non_len):
            v = non_s_core[j]
            candidate_edges.append((u, v))

    # core <-> non-core
    for u in non_s_core:
        for v in s_core:
            candidate_edges.append((u, v))

    return candidate_edges

def make_candidate_nodes_v2(G_prime, nodes, s, budget_left, coreness, upperbound, UT, T2_upperbound, delta_Tactic):
    '''
    When using the upper-bound strategy, construct 'candidate_nodes' to apply the bound, 
    even if the self-edge strategy is not used. 
    Initialize the s-core upper bound to 0 and sort the candidates.
    '''
    candidate_nodes = []

    for u in nodes:
        if s - coreness.get(u, (s, 0))[0] > budget_left:
            continue

        if not G_prime.nodes[u]['label']:
            temp_start = time.time()
            upperbound[u] = functions.Upperbound(G_prime, u, coreness, s, delta_Tactic)
            temp_end = time.time()
            UT += temp_end - temp_start
        else:
            if T2_upperbound:
                upperbound[u] = 0
        candidate_nodes.append(u)
    
    candidate_nodes.sort(key = lambda x : -upperbound[x])

    return candidate_nodes

def iteration_edges_no_upperbound(G_prime, candidate_edges, s, b, t, spent, coreness, FT, T3_reuse, delta_Tactic, comp_of={}, best=(None, 0, 0.0, 0)):
    '''
    naive algorithm.
    '''

    # initial setting
    best_edge, best_delta, most_FR, most_follower = best

    for (u,v) in candidate_edges:
        if T3_reuse:
            if comp_of[u] == comp_of[v]:
                continue
            
        e = (u, v)
        if delta_Tactic == "compute":
            delta_e = functions.computeDelta(G_prime, s, e, t, coreness)
        else:
            if G_prime.has_edge(u, v):
                delta_e = s - G_prime[u][v]['weight']
            else:
                delta_e = s

        if delta_e > 0 and spent + delta_e <= b:
            # calculate the follower
            temp_start = time.time()
            followers = functions.FindFollowers(e, delta_e, G_prime, s, coreness)
            experiment_iter.GLOBAL_CNT += 1
            temp_end = time.time()
            FT += temp_end - temp_start

            FR = len(followers) / delta_e  # follower rate

            # renew the maximum value
            if e[0] > e[1]:
                e = (e[1], e[0])
            if FR > most_FR or (FR == most_FR and (best_edge is None or e < best_edge)):
                best_edge = e
                best_delta = delta_e
                most_FR = FR
                most_follower = len(followers)
    return best_edge, best_delta, most_FR, most_follower

def save_result_to_csv(A, s_core_num, total_time, memory, args, total_follower, num_iter):
    write_header = not os.path.exists(args.output_path)
    data = args.network.split('/')[4]
    with open(args.output_path, 'a', newline='') as f:
        writer = csv.writer(f)
        if write_header:
            writer.writerow(["data", "s", "b", "tactics", "num_s_core", "follower", "num_iter", "num_anchored_edges", "iter", "delta", "num_follower", "self", "reinforce", "total_time", "memory"])
        for idx, (_, delta, _, num_follower, self, reinforce) in enumerate(A, start=1):
            writer.writerow(
                [data, args.s, args.b, args.tactics, s_core_num, total_follower, num_iter, len(A), idx, delta, num_follower, self, reinforce, total_time, memory])
    print(f"Saved results to {data, args.s, args.b, args.tactics}")