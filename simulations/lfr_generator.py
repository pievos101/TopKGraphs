# Save as lfr_generator.py or just use inline via reticulate

import networkx as nx

def generate_lfr(
    n_nodes=50,
    tau1=2,
    tau2=1,
    mu=0.05,
    avg_degree=3,
    max_degree=5,
    min_community=5,
    max_community=20,
    seed=None
):
    """
    Generate an LFR benchmark graph using NetworkX.
    
    Returns:
        - G: networkx graph
        - membership: dict of node -> community
    """
    G = nx.LFR_benchmark_graph(
        n=n_nodes,
        tau1=tau1,
        tau2=tau2,
        mu=mu,
        average_degree=avg_degree,
        max_degree=max_degree,
        min_community=min_community,
        max_community=max_community,
        seed=seed
    )
    
    # Extract membership
    membership = {n: list(G.nodes[n]['community'])[0] 
                  if isinstance(G.nodes[n]['community'], set)
                  else G.nodes[n]['community'] 
                  for n in G.nodes}
    
    return G, membership