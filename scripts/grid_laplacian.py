from pathlib import Path
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_PROCESSED = REPO_ROOT / "data"
OUTPUTS_FIGURES = REPO_ROOT / "figures"
LATEX_PAPER2_FIGURES = REPO_ROOT / "figures"

DATA_PROCESSED.mkdir(parents=True, exist_ok=True)
OUTPUTS_FIGURES.mkdir(parents=True, exist_ok=True)
LATEX_PAPER2_FIGURES.mkdir(parents=True, exist_ok=True)

# Haversine formula to compute distance in km between two GPS coordinates
def haversine_distance(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    """Compute great-circle distance in kilometers using the Haversine formula."""
    R = 6371.0  # Earth radius in km
    
    phi1 = np.radians(lat1)
    phi2 = np.radians(lat2)
    delta_phi = np.radians(lat2 - lat1)
    delta_lambda = np.radians(lon2 - lon1)
    
    a = np.sin(delta_phi/2.0)**2 + np.cos(phi1) * np.cos(phi2) * np.sin(delta_lambda/2.0)**2
    c = 2.0 * np.arctan2(np.sqrt(a), np.sqrt(1.0 - a))
    
    return R * c

def build_bihar_grid() -> None:
    """Construct the BSPTCL Bihar transmission grid, compute Laplacian spectrum, and export topology."""
    print("================================================================================")
    print("      BSPTCL Bihar Transmission Grid Topology & Laplacian Connectivity         ")
    print("================================================================================")
    
    # 1. Define Nodes (substations & generation hubs in Bihar)
    # Mapping GPS coordinates: (Latitude, Longitude)
    nodes = {
        "Patna_Digha":     {"coords": (25.6418, 85.0976), "type": "Load_GSS", "color": "skyblue"},
        "Patna_Gaurichak": {"coords": (25.5235, 85.1985), "type": "Load_GSS", "color": "skyblue"},
        "Patna_Fatuha":    {"coords": (25.5085, 85.3125), "type": "Load_GSS", "color": "skyblue"},
        "Patna_Khagaul":   {"coords": (25.5785, 85.0452), "type": "Load_GSS", "color": "skyblue"},
        "Bihta_GSS":       {"coords": (25.5601, 84.8710), "type": "Load_GSS", "color": "skyblue"},
        
        "NTPC_Barh":       {"coords": (25.4775, 85.7092), "type": "Thermal_Gen", "color": "crimson"},
        "NTPC_Kahalgaon":  {"coords": (25.2677, 87.2187), "type": "Thermal_Gen", "color": "crimson"},
        "NTPC_Nabinagar":  {"coords": (24.6155, 84.1488), "type": "Thermal_Gen", "color": "crimson"},
        "Barauni_Thermal": {"coords": (25.4248, 86.0204), "type": "Thermal_Gen", "color": "crimson"},
        "Kanti_Thermal":   {"coords": (26.2009, 85.2647), "type": "Thermal_Gen", "color": "crimson"},
        
        "Kajra_Solar":     {"coords": (25.1764, 86.1664), "type": "Solar_Gen", "color": "gold"},
        "Banka_Solar":     {"coords": (24.8785, 86.9242), "type": "Solar_Gen", "color": "gold"},
        "Jamui_Solar":     {"coords": (24.9221, 86.2235), "type": "Solar_Gen", "color": "gold"},
        
        "Bodhgaya_GSS":    {"coords": (24.6958, 84.9913), "type": "Load_GSS", "color": "skyblue"},
        "Gaya_GSS":        {"coords": (24.7964, 85.0076), "type": "Load_GSS", "color": "skyblue"},
        "Dehri_GSS":       {"coords": (24.9121, 84.1865), "type": "Load_GSS", "color": "skyblue"},
        "Biharsharif_GSS": {"coords": (25.1982, 85.5149), "type": "Load_GSS", "color": "skyblue"},
        "Begusarai_GSS":   {"coords": (25.4182, 86.1272), "type": "Load_GSS", "color": "skyblue"},
        "Bhagalpur_GSS":   {"coords": (25.2425, 86.9718), "type": "Load_GSS", "color": "skyblue"},
        "Darbhanga_GSS":   {"coords": (26.1542, 85.8918), "type": "Load_GSS", "color": "skyblue"},
        "Samastipur_GSS":  {"coords": (25.8622, 85.7797), "type": "Load_GSS", "color": "skyblue"},
        "Kishanganj_GSS":  {"coords": (26.2652, 87.9409), "type": "Load_GSS", "color": "skyblue"},
        "Gopalganj_GSS":   {"coords": (26.4682, 84.4418), "type": "Load_GSS", "color": "skyblue"},
        "Nawada_GSS":      {"coords": (24.8894, 85.5414), "type": "Load_GSS", "color": "skyblue"}
    }
    
    # 2. Define Edges (corridors connecting substations circles)
    transmission_lines = [
        ("NTPC_Barh", "Patna_Fatuha"),
        ("NTPC_Barh", "Biharsharif_GSS"),
        ("NTPC_Barh", "Begusarai_GSS"),
        
        ("NTPC_Kahalgaon", "Bhagalpur_GSS"),
        ("NTPC_Kahalgaon", "Banka_Solar"),
        ("NTPC_Kahalgaon", "Kishanganj_GSS"),
        
        ("NTPC_Nabinagar", "Dehri_GSS"),
        ("NTPC_Nabinagar", "Gaya_GSS"),
        
        ("Barauni_Thermal", "Begusarai_GSS"),
        ("Barauni_Thermal", "Samastipur_GSS"),
        ("Barauni_Thermal", "Biharsharif_GSS"),
        
        ("Kanti_Thermal", "Darbhanga_GSS"),
        ("Kanti_Thermal", "Gopalganj_GSS"),
        ("Kanti_Thermal", "Patna_Digha"),
        
        # Patna Ring Topology
        ("Patna_Digha", "Patna_Khagaul"),
        ("Patna_Digha", "Patna_Gaurichak"),
        ("Patna_Gaurichak", "Patna_Fatuha"),
        ("Patna_Fatuha", "Biharsharif_GSS"),
        ("Patna_Khagaul", "Bihta_GSS"),
        ("Bihta_GSS", "Dehri_GSS"),
        
        # Southern & Eastern Corridors
        ("Bodhgaya_GSS", "Gaya_GSS"),
        ("Bodhgaya_GSS", "Nawada_GSS"),
        ("Bodhgaya_GSS", "Jamui_Solar"),
        
        ("Kajra_Solar", "Biharsharif_GSS"),
        ("Kajra_Solar", "Begusarai_GSS"),
        
        ("Banka_Solar", "Bhagalpur_GSS"),
        ("Banka_Solar", "Jamui_Solar"),
        
        ("Jamui_Solar", "Nawada_GSS"),
        ("Darbhanga_GSS", "Samastipur_GSS"),
        ("Samastipur_GSS", "Begusarai_GSS"),
        ("Nawada_GSS", "Biharsharif_GSS")
    ]
    
    # 3. Construct Graph
    G = nx.Graph()
    
    # Add nodes with metadata
    for node_name, info in nodes.items():
        G.add_node(node_name, coords=info["coords"], type=info["type"], color=info["color"])
        
    # Calculate physical coupling strengths as edge weights
    # Reactance X_ij = 0.32 * distance_km
    # Coupling K_ij = 50 / X_ij (normalized factor)
    edges_data = []
    for u, v in transmission_lines:
        lat1, lon1 = nodes[u]["coords"]
        lat2, lon2 = nodes[v]["coords"]
        
        dist = haversine_distance(lat1, lon1, lat2, lon2)
        reactance = 0.32 * dist
        weight = 50.0 / reactance  # Coupling strength
        
        G.add_edge(u, v, weight=weight, distance=dist)
        edges_data.append({"From": u, "To": v, "Distance_km": dist, "Coupling_K": weight})
        
    print(f"Nodes in Network: {G.number_of_nodes()}")
    print(f"Lines in Network: {G.number_of_edges()}")
    
    # 4. Construct Weighted Laplacian Matrix
    print("\n--- Computing Graph Laplacian Matrix ---")
    nodes_list = list(G.nodes())
    A = nx.adjacency_matrix(G, nodelist=nodes_list, weight='weight').toarray()
    
    # Degree Matrix D
    D = np.diag(np.sum(A, axis=1))
    
    # Laplacian Matrix L = D - A
    L = D - A
    
    # 5. Calculate Eigenvalues of Graph Laplacian
    eigenvalues = np.sort(np.linalg.eigvals(L))
    
    # Under standard theorem: second smallest eigenvalue (Fiedler eigenvalue) indicates synchronizability
    fiedler_val = eigenvalues[1]
    print("\nGraph Laplacian Eigenvalues:")
    for j, val in enumerate(eigenvalues[:6], 1):
        print(f"  lambda_{j} = {val:.6f}")
    print("  ...")
    print(f"Fiedler Eigenvalue (lambda_2, Algebraic Connectivity) = {fiedler_val:.6f}")
    
    if fiedler_val > 0:
        print("\nCRITICAL RESULT: lambda_2 > 0 => The BSPTCL network is mathematically GUARANTEED to achieve complete synchronization under sufficient coupling!")
    else:
        print("\nWARNING: lambda_2 <= 0 => The network graph is disconnected! Grid synchronization is physically impossible.")
        
    # Save network data to CSV
    df_edges = pd.DataFrame(edges_data)
    edges_csv = DATA_PROCESSED / "bihar_grid_edges.csv"
    df_edges.to_csv(edges_csv, index=False)
    print(f"\nSaved transmission lines parameters to: {edges_csv}")
    
    # 6. Visualize Geographic Network Graph
    print("\nPlotting geographic BSPTCL power grid topology...")
    plt.figure(figsize=(12, 10))
    
    # Plot background outline or grid
    pos = {node: (info["coords"][1], info["coords"][0]) for node, info in nodes.items()}  # (Lon, Lat) to match geographic representation
    
    # Get node colors
    node_colors = [info["color"] for info in nodes.values()]
    
    # Draw edges with line thickness proportional to coupling strength
    weights = [G[u][v]['weight'] for u, v in G.edges()]
    max_w = max(weights)
    edge_widths = [1.5 + 4.0 * (w / max_w) for w in weights]
    
    # Draw nodes
    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=350, edgecolors='black', linewidths=1.0)
    nx.draw_networkx_edges(G, pos, width=edge_widths, edge_color='gray', alpha=0.8)
    
    # Draw labels with small vertical offset for legibility
    label_pos = {node: (coords[0], coords[1] + 0.04) for node, coords in pos.items()}
    nx.draw_networkx_labels(G, label_pos, font_size=8, font_family='sans-serif', font_weight='bold')
    
    # Add custom legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', label='Grid Substations (Load GSS)', markerfacecolor='skyblue', markersize=10, markeredgecolor='black'),
        Line2D([0], [0], marker='o', color='w', label='NTPC Thermal Power plants', markerfacecolor='crimson', markersize=10, markeredgecolor='black'),
        Line2D([0], [0], marker='o', color='w', label='Renewable Solar Plants', markerfacecolor='gold', markersize=10, markeredgecolor='black')
    ]
    plt.legend(handles=legend_elements, loc='lower right', frameon=True, facecolor='white', edgecolor='black')
    
    plt.title("Bihar BSPTCL Intra-State Power Transmission Network Map (132kV / 220kV GSS)", fontsize=14, fontweight='bold')
    plt.xlabel("Longitude (°E)", fontsize=12)
    plt.ylabel("Latitude (°N)", fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.5)
    
    fig_file = OUTPUTS_FIGURES / "bihar_grid_topology.png"
    fig_paper2 = LATEX_PAPER2_FIGURES / "bihar_grid_topology.png"
    plt.savefig(fig_file, dpi=300)
    plt.savefig(fig_paper2, dpi=300)
    plt.close()
    
    print(f"Successfully saved BSPTCL grid visualization at: {fig_file} and {fig_paper2}")

if __name__ == "__main__":
    build_bihar_grid()
