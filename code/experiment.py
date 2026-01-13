import networkx as nx
from itertools import combinations
import time
import sys

import functions
import exp_func

def run(G, s, b, t, T1_self_edge = True, T2_upperbound = True, T3_reuse = True, delta_Tactic = "compute"):
    # Tactic 1 self_edge and Tactic 2 upperbound as a parameter

    total_follower = 0

    FT = 0.0    # Follower computing time
    UT = 0.0    # Upperbound computing time

    G_prime = G.copy()
    A = []  # The anchored edges
    
    # Calculate s-core, coreness, layer
    # The information whether the node is in s-core is stored as a attribute in graph.
    coreness = {}
    s_core_num = functions.calculate_s_core(G_prime, G_prime.nodes, s, coreness)

    spent = 0  # the budget used
    
    # Using self_edge theorem, s-core are permanently removed in candidates.
    s_cand = None   # Only store one node in s-core for anchoring self edge.
    if T1_self_edge:
        non_s_core, s_cand = exp_func.self_edge_pruning(G_prime) # pruned set

    upperbound = {}

    # Using reuse theorem, store the best edge of each connected component.
    if T3_reuse:
        comp_of, nodes_in, intra_best, s_cand = build_initial_caches(G_prime, s, t, b, spent, coreness, upperbound, T1_self_edge, T2_upperbound, UT, FT, delta_Tactic)
    else:
        comp_of = {}
    

    ### ITERATION START
    it = 0
    while spent < b:
        '''
        Considering self_edge tactic and upperbound tactic, two procedures below are differently progressed.
        1. Making Candidate
        2. Iterating Candidate to find best edge
        '''
        
        # init value
        best = (None, 0, 0.0, 0)
        # Reuse Logic : First, get best edge using Cache. Second, edges in same component can be pruned.
        if T3_reuse:
            best_intra = max(intra_best.values(), key=lambda x:x[2], default=None)
            if best_intra:
                best = best_intra
        
        if T1_self_edge:
            # 1. Making Candidate
            candidate_nodes = exp_func.make_candidate_nodes(G_prime, non_s_core, s, b-spent, coreness, upperbound, UT, T2_upperbound, delta_Tactic)

            # 2. Iterating Candidate to find best edge
            if T2_upperbound:
                # tactics TT
                best = exp_func.iteration_nodes_upperbound(G_prime, candidate_nodes, s, b, t, spent, coreness, upperbound, s_cand, FT, T1_self_edge, T3_reuse, delta_Tactic, comp_of, best)
            else:
                # tactics TF
                best = exp_func.iteration_nodes_no_upperbound(G_prime, candidate_nodes, s, b, t, spent, coreness, s_cand, FT, T3_reuse, delta_Tactic, comp_of, best)
        else:
            if T2_upperbound:
                # tactics FT
                candidate_nodes = exp_func.make_candidate_nodes_v2(G_prime, G_prime.nodes, s, b-spent, coreness, upperbound, UT, T2_upperbound, delta_Tactic)
                best = exp_func.iteration_nodes_upperbound(G_prime, candidate_nodes, s, b, t, spent, coreness, upperbound, s_cand, FT, T1_self_edge, T3_reuse, delta_Tactic, comp_of, best)
            else:
                # tactics FF
                candidate_edges = exp_func.make_candidate_edges(G_prime, G_prime.nodes, s, b-spent, coreness)
                best = exp_func.iteration_edges_no_upperbound(G_prime, candidate_edges, s, b, t, spent, coreness, FT, T3_reuse, delta_Tactic, comp_of, best)

        # No more possible anchor edge
        if not best[0]: break

        best_edge, best_delta, most_FR, most_follower = best
        (u,v) = best_edge

        # Left budget cannot handle delta (only edge from reuse logic)
        if best_delta > b - spent:
            if best == best_intra:
                if G_prime.nodes[u]['label']:
                    c = comp_of[v]
                else:
                    c = comp_of[u]
                intra_best.pop(c, None)

                # Renew the connected component where the best edge came from using budget_left
                edge2, d2, fr2, most_follower = find_intra_best(G_prime, nodes_in[c], coreness, s, t, b, spent, upperbound, T1_self_edge, T2_upperbound, UT, FT, s_cand, delta_Tactic)
                if edge2 is not None:
                    intra_best[c] = (edge2, d2, fr2, most_follower)

            else:
                # In the iteration, we check budget in selecting candidates.
                print("[ERROR] Cannot happen. Some error detected")

            continue 
    
        # debugging 3
        # print()
        # print(best)

        ### Anchor best edge
        # Consider self edge
        if u == v:
            v = s_cand

        # Update Graph
        if G_prime.has_edge(u,v):
            G_prime[u][v]['weight'] += best_delta
        else:
            G_prime.add_edge(u,v, weight=best_delta)
        
        # add budget
        spent += best_delta

        # add answer
        total_follower += most_follower
        A.append(((u,v), best_delta, most_FR, most_follower))

        # Update coreness (When using reuse tactic, update it locally)
        if T3_reuse:
            # ----- Renew Cache --------------------------------
            # intra vs inter
            if best == best_intra:
                # same CC
                if G_prime.nodes[u]['label']:
                    c = comp_of[v]
                else:
                    c = comp_of[u]

                # calculate locally
                s_core_num = functions.calculate_s_core(G_prime, nodes_in[c], s, coreness)

                new_nodes = {v for v in nodes_in[c] if not G_prime.nodes[v]['label']}

                # Check whteher non-s-core is left in CC
                if not new_nodes:
                    # Delete Cache
                    invalidate({c}, intra_best)
                    nodes_in.pop(c)
                else:
                    # Using self_edge tactic, only store non-s-core in CC.
                    if T1_self_edge:
                        nodes_in[c] = new_nodes
                    
                    # Change Cache with new value
                    invalidate({c}, intra_best)
                    edge, delta2, FR2, most_follower = find_intra_best(G_prime, nodes_in[c], coreness, s, t, b, spent, upperbound, T1_self_edge, T2_upperbound, UT, FT, s_cand, delta_Tactic)
                    if edge:
                        intra_best[c] = (edge, delta2, FR2, most_follower)

            else:
                # inter : Union CC
                c1, c2 = comp_of[u], comp_of[v]
                union_nodes = nodes_in[c1] | nodes_in[c2]

                # calculate locally
                s_core_num = functions.calculate_s_core(G_prime, union_nodes, s, coreness)

                new_nodes = {x for x in union_nodes if not G_prime.nodes[x]['label']}

                # Delete old Cache
                invalidate({c1,c2}, intra_best)
                nodes_in.pop(c1); nodes_in.pop(c2)

                # # Check whteher non-s-core is left in CC
                if new_nodes:
                    new_c = max(nodes_in)+1  # New id
                    nodes_in[new_c] = new_nodes
                    for x in new_nodes:
                        comp_of[x] = new_c

                    # New Cache
                    # G, nodes, coreness, s, t, b, spent, upperbound, T1_self_edge, T2_upperbound, UT, FT, s_cand
                    edge, delta2, FR2, most_follower = find_intra_best(G_prime, new_nodes, coreness, s, t, b, spent, upperbound, T1_self_edge, T2_upperbound, UT, FT, s_cand, delta_Tactic)
                    if edge:
                        intra_best[new_c] = (edge, delta2, FR2, most_follower)
        else:   # Not using reuse tactic. Calculate globally.
            # calculate s-core again
            coreness = {}
            s_core_num = functions.calculate_s_core(G_prime, G_prime.nodes, s, coreness)
    return A, FT, UT, G_prime, total_follower


def build_initial_caches(G, s, t, b, spent, coreness, upperbound, T1_self_edge, T2_upperbound, UT, FT, delta_Tactic):
    '''
    If self_edge tactic is used, only non-s-cores are saved in Connected Component.
    If not, s-core has to be saved in Connected Component.
    '''
    comp_of, nodes_in = {}, {}
    intra_best = {}
    s_cand = None

    if T1_self_edge:    # Only non-s-core
        # Extract non-s-core nodes
        non_s_core_nodes = set()
        for v in G.nodes:
            if not G.nodes[v]['label']:
                non_s_core_nodes.add(v)
            else:
                if s_cand is None or (s_cand is not None and v < s_cand):
                    s_cand = v
        
        # Create subgraph with only non-s-core nodes
        non_s_core_subgraph = G.subgraph(non_s_core_nodes)

        # Decompose CC
        for cid, nodes in enumerate(nx.connected_components(non_s_core_subgraph)):

            nodes_in[cid] = set(nodes)
            for v in nodes:
                comp_of[v] = cid

        # We need at least one s_core           
        if s_cand is None:
            print("No node in s-core. Change s value")
            sys.exit(1)

    else:   # Store every nodes
        # Decompose CC
        for cid, nodes in enumerate(nx.connected_components(G)):
            nodes_in[cid] = nodes
            for v in nodes:
                comp_of[v] = cid

    # intra
    for cid, nodes in nodes_in.items():
        edge, delta, FR, most_follower = find_intra_best(G, nodes, coreness, s, t, b, spent, upperbound, T1_self_edge, T2_upperbound, UT, FT, s_cand, delta_Tactic)
        if edge:
            intra_best[cid] = (edge, delta, FR, most_follower)

    return comp_of, nodes_in, intra_best, s_cand


def find_intra_best(G, nodes, coreness, s, t, b, spent, upperbound, T1_self_edge, T2_upperbound, UT, FT, s_cand, delta_Tactic):

    if T1_self_edge:
        candidate_nodes = exp_func.make_candidate_nodes(G, nodes, s, b-spent, coreness, upperbound, UT, T2_upperbound, delta_Tactic)

        if T2_upperbound:
            # tactics TT
            best = exp_func.iteration_nodes_upperbound(G, candidate_nodes, s, b, t, spent, coreness, upperbound, s_cand, FT, T1_self_edge, False, delta_Tactic)
        else:
            # tactics TF
            best = exp_func.iteration_nodes_no_upperbound(G, candidate_nodes, s, b, t, spent, coreness, s_cand, FT, False, delta_Tactic)
    
    else:
        if T2_upperbound:
            # tactics FT
            candidate_nodes = exp_func.make_candidate_nodes_v2(G, nodes, s, b-spent, coreness, upperbound, UT, T2_upperbound, delta_Tactic)
            best = exp_func.iteration_nodes_upperbound(G, candidate_nodes, s, b, t, spent, coreness, upperbound, s_cand, FT, T1_self_edge, False, delta_Tactic)
        else:
            # tactics FF
            candidate_edges = exp_func.make_candidate_edges(G, nodes, s, b-spent, coreness)
            best = exp_func.iteration_edges_no_upperbound(G, candidate_edges, s, b, t, spent, coreness, FT, False, delta_Tactic)
    
    return best   # best_edge is None → No appropriate edge in this CC


def invalidate(ccs, intra_best):
    for c in ccs:
        intra_best.pop(c, None)