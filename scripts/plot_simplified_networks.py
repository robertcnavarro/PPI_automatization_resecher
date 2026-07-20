import os
import networkx as nx
import matplotlib.pyplot as plt
import glob

RESULTS_DIR = "results"
CYTOSCAPE_DIR = os.path.join(RESULTS_DIR, "interactomes", "cytoscape")
FIGURES_DIR = os.path.join(RESULTS_DIR, "figures")
os.makedirs(FIGURES_DIR, exist_ok=True)

def plot_network(sif_path, timepoint):
    if not os.path.exists(sif_path):
        return
        
    G = nx.Graph()
    with open(sif_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                node_a, rel, node_b = parts[0], parts[1], parts[2]
                G.add_edge(node_a, node_b)
                
    if len(G) == 0:
        return

    # Find Top 30 Host Hubs + All RSV nodes
    degrees = dict(G.degree())
    rsv_nodes = [n for n in G.nodes() if n.startswith("RSV")]
    host_nodes = [n for n in G.nodes() if not n.startswith("RSV")]
    
    print(f"[{timepoint}] Found viral nodes: {rsv_nodes}")
    
    # Keep top 30 host nodes by degree to prevent hairball
    top_host = sorted(host_nodes, key=lambda x: degrees[x], reverse=True)[:30]
    
    # Also explicitly include the neighbors of the RSV nodes so their edges are drawn!
    for v in rsv_nodes:
        for neighbor in G.neighbors(v):
            if neighbor not in top_host:
                top_host.append(neighbor)
    
    sub_nodes = set(rsv_nodes + top_host)
    H = G.subgraph(sub_nodes)
    
    plt.figure(figsize=(10, 8))
    pos = nx.spring_layout(H, k=0.8, iterations=70, seed=42)
    
    # Draw edges
    nx.draw_networkx_edges(H, pos, alpha=0.3, edge_color='gray')
    
    # Draw host nodes, sizing by their degree in the FULL network
    h_nodes = [n for n in H.nodes() if not n.startswith("RSV")]
    h_sizes = [min(degrees[n]*20, 2000) + 100 for n in h_nodes]
    nx.draw_networkx_nodes(H, pos, nodelist=h_nodes, node_color='#3b82f6', node_size=h_sizes, alpha=0.8)
    
    # Draw viral nodes (Fixed large size)
    v_nodes = [n for n in H.nodes() if n.startswith("RSV")]
    if v_nodes:
        nx.draw_networkx_nodes(H, pos, nodelist=v_nodes, node_color='#ef4444', node_size=800, alpha=0.9, edgecolors='black', linewidths=2)
    
    # Draw labels
    nx.draw_networkx_labels(H, pos, font_size=10, font_weight='bold', font_family="sans-serif", font_color='black')
    
    plt.title(f"Red Simplificada de Interacción VSR-Huésped ({timepoint})\nTop 30 Hubs Celulares", fontsize=16, fontweight='bold', pad=20)
    plt.axis('off')
    
    out_path = os.path.join(FIGURES_DIR, f"network_preview_{timepoint}.png")
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"✅ Generada visualización para {timepoint}: {out_path}")

def main():
    print("Generando previsualizaciones de redes para el Dashboard...")
    for tp in ["12hrs", "24hrs", "48hrs", "72hrs", "120hrs"]:
        sif_path = os.path.join(CYTOSCAPE_DIR, f"network_{tp}.sif")
        plot_network(sif_path, tp)

if __name__ == "__main__":
    main()
