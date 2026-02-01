# ------------------------------
# soft_topkgraphs.py
# ------------------------------

import torch
import torch.nn as nn
import torch.nn.functional as F
import networkx as nx

class SoftTopKGraphs(nn.Module):
    """
    Jaccard-anchored soft TopKGraphs with learnable edge scaling.
    """
    def __init__(self, G, walk_length=10, eps=1e-3):
        """
        G : networkx.Graph (undirected)
        walk_length : number of soft propagation steps
        eps : additive smoothing for Jaccard
        """
        super().__init__()
        self.G = G
        self.n = G.number_of_nodes()
        self.walk_length = walk_length
        self.eps = eps

        # Precompute adjacency list
        self.adj_list = [list(G.neighbors(i)) for i in range(self.n)]

        # Learnable edge scaling per node
        # Each node has a parameter for scaling outgoing edges
        self.edge_params = nn.Parameter(torch.zeros(self.n))
    
    def fixed_jaccard(self, anchor):
        """
        Compute Jaccard similarity vector of all nodes to the anchor node
        """
        anchor_neighbors = set(self.adj_list[anchor])
        jaccard = torch.zeros(self.n)
        for v in range(self.n):
            neighbors_v = set(self.adj_list[v])
            inter = len(anchor_neighbors & neighbors_v)
            union = len(anchor_neighbors | neighbors_v)
            jaccard[v] = inter / union if union > 0 else 0.0
        return jaccard + self.eps

    def soft_walk(self, anchor):
        """
        Anchor-conditioned soft walk with learnable scaling
        """
        P = torch.zeros(self.n, self.n)  # transition matrix
        jaccard = self.fixed_jaccard(anchor)

        for u in range(self.n):
            neighbors_u = self.adj_list[u]
            if len(neighbors_u) == 0:
                continue
            # Anchor-anchored + learnable scaling
            scores = jaccard[neighbors_u] * torch.exp(self.edge_params[u])
            P[u, neighbors_u] = scores / scores.sum()

        # Initialize soft visitation
        h = torch.zeros(self.n)
        h[anchor] = 1.0
        V_cum = h.clone()

        # Propagate
        for _ in range(self.walk_length):
            h = P @ h
            V_cum += h

        return V_cum / V_cum.sum()

    def forward(self):
        """
        Generate full embeddings: row = anchor node soft visitation
        Returns:
            embeddings : n x n tensor
        """
        embeddings = torch.zeros(self.n, self.n)
        for anchor in range(self.n):
            embeddings[anchor, :] = self.soft_walk(anchor)
        return embeddings

# ------------------------------
# Contrastive learning loss (InfoNCE)
# ------------------------------
def contrastive_loss(embeddings, temperature=0.5):
    """
    embeddings : n x d tensor
    """
    norm_emb = F.normalize(embeddings, dim=1)
    sim_matrix = (norm_emb @ norm_emb.t()) / temperature
    pos = torch.diag(sim_matrix)
    loss = -torch.log(torch.exp(pos) / torch.exp(sim_matrix).sum(dim=1))
    return loss.mean()

# ------------------------------
# Example usage
# ------------------------------
if __name__ == "__main__":
    G = nx.cycle_graph(5)
    model = SoftTopKGraphs(G, walk_length=5)
    optimizer = torch.optim.Adam(model.parameters(), lr=0.01)

    for epoch in range(50):
        optimizer.zero_grad()
        emb = model()
        loss = contrastive_loss(emb)
        loss.backward()
        optimizer.step()
        if epoch % 10 == 0:
            print(f"Epoch {epoch}, Loss: {loss.item():.4f}")

    print("Final embeddings:\n", emb.detach().numpy())