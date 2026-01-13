import functions
import networkx as nx
import random
import sys
def run(G, s, b, t, compare_Tactic, delta_tatic):
    G_prime = G.copy()
    A = []

    # Calculate initial s-core and coreness
    coreness = {}
    s_core_num = functions.calculate_s_core(G_prime, G_prime.nodes, s, coreness)
    first = s_core_num
    
    # Partition nodes into s-core and non-s-core sets
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
        
    # Select candidates: strictly nodes in the (s-1)-shell
    candidates = [u for u in non_s_core if coreness[u][0] == s - 1]
    if not candidates or b <= 0:
        return A, 0
    
    # Define scoring functions for ranking
    def score_degree(G_, u, coreness_, s_):
        return G_.degree(u)
    def score_weight_sum(G_, u, coreness_, s_):
        return sum(data.get('weight', 1) for _, _, data in G_.edges(u, data=True))
    def score_high_layer_degree(G_, u, coreness_, s_):
        return sum(1 for w in G_.neighbors(u) if coreness_.get(w, (s_, 0)) > coreness_.get(u, (s_, 0)))
    def score_high_layer_weight(G_, u, coreness_, s_):
        return sum(data.get('weight', 1)
                for _, w, data in G_.edges(u, data=True)
                if coreness_.get(w, (s_, 0)) > coreness_.get(u, (s_, 0)))
    if compare_Tactic == 'degree':
        score_function = score_degree
    elif compare_Tactic == 'high_degree':
        score_function = score_high_layer_degree
    elif compare_Tactic == 'weight_sum':
        score_function = score_weight_sum
    elif compare_Tactic == 'high_weight_sum':
        score_function = score_high_layer_weight
    else:
        score_function = None  # random
        
    # Sort and select top candidates based on the scoring metric
    if score_function is None:  # random
        random.shuffle(candidates)
        top_nodes = candidates[:min(b, len(candidates))]
    else:
        candidates.sort(
            key=lambda u: score_function(G_prime, u, coreness, s),
            reverse=True
        )
        top_nodes = candidates[:min(b, len(candidates))]
        
    # Connect selected nodes to 's_cand' (increment delta=1)
    budget_used = 0
    for u in top_nodes:
        if budget_used + 1 > b:
            break
        delta = 1
        
        # Update graph structure
        if G_prime.has_edge(u, s_cand):
            G_prime[u][s_cand]['weight'] += delta
        else:
            G_prime.add_edge(u, s_cand, weight=delta)
        budget_used += delta
        
        # Record the modification results
        A.append(((u, s_cand), delta, None, None))
        
    # Re-calculate final s-core size
    coreness = {}
    final_size = functions.calculate_s_core(G_prime, G_prime.nodes, s, coreness)
    return len(A), final_size, final_size - first