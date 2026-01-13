import functions
import networkx as nx
import matplotlib.pyplot as plt
import sys
from itertools import combinations_with_replacement
from collections import Counter

def run(G, s, b, t):

    G_prime = G.copy()
    A = set()  # the set of (increased edge, delta) pair
    
    coreness = {}
    s_core_num = functions.calculate_s_core(G_prime, G_prime.nodes, s, coreness)  # Calculate s-core and coreness

    # Filter candidate edges
    non_s_core = []
    s_cand = None
    for n, d in G_prime.nodes(data=True):
        if not d['label']:
            non_s_core.append(n)
        else:
            if not s_cand:
                s_cand = n
    if s_cand is None:
        print("No node in s-core. Change s value")
        sys.exit(1)

    # generate candidate edge
    candidate_edges = []
    # 1) Between non-core nodes
    for i in range(len(non_s_core)):
        for j in range(i+1, len(non_s_core)):
            candidate_edges.append((non_s_core[i], non_s_core[j]))
    # 2) Between non-core nodes and s-core nodes
    for u in non_s_core:
        candidate_edges.append((u, s_cand))

    C = list(combinations_with_replacement(candidate_edges, b))

    best_gain = -1
    best_combo = None

    # Iterate over all combinations
    for combo in C:
        combo_map = Counter(combo)  # e.g. Counter({(u,v):2, (x,y):1, ...})
        # Create a copy for s-core calculation
        G_temp = G.copy()
        for (u,v), delta in combo_map.items():
            if G_temp.has_edge(u, v):
                G_temp[u][v]['weight'] += delta
            else:
                G_temp.add_edge(u, v, weight=delta)

        # Re-calculate s-core on the modified G_temp
        temp_coreness = {}
        new_size = functions.calculate_s_core(G_temp, G_temp.nodes, s, temp_coreness)

        gain = new_size - s_core_num
        if gain > best_gain:
            best_gain  = gain
            best_combo = combo_map

    # Apply the best combination
    for (u, v), delta in best_combo.items():
        if G_prime.has_edge(u, v):
            G_prime[u][v]['weight'] += delta
        else:
            G_prime.add_edge(u, v, weight=delta)
        A.add(((u, v), delta))

    # Final s-core size
    final_coreness = {}
    final_size = functions.calculate_s_core(G_prime, G_prime.nodes, s, final_coreness)

    return A, final_size, final_size-s_core_num