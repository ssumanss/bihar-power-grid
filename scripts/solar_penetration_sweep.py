from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_PROCESSED = REPO_ROOT / "data"
OUTPUTS_FIGURES = REPO_ROOT / "figures"
LATEX_PAPER2_FIGURES = REPO_ROOT / "figures"

DATA_PROCESSED.mkdir(parents=True, exist_ok=True)
OUTPUTS_FIGURES.mkdir(parents=True, exist_ok=True)
LATEX_PAPER2_FIGURES.mkdir(parents=True, exist_ok=True)

def swing_ode(t: float, y: np.ndarray, num_nodes: int, M: np.ndarray, D: np.ndarray, P: np.ndarray, K_matrix: np.ndarray) -> np.ndarray:
    """Coupled Kuramoto-swing ODE system computing phase angle and frequency derivatives."""
    theta = y[:num_nodes]
    omega = y[num_nodes:]
    dtheta_dt = omega
    domega_dt = np.zeros(num_nodes)
    for i in range(num_nodes):
        coupling_sum = 0.0
        for j in range(num_nodes):
            if K_matrix[i, j] > 0:
                coupling_sum += K_matrix[i, j] * np.sin(theta[i] - theta[j])
        domega_dt[i] = (1.0 / M[i]) * (P[i] - D[i] * omega[i] - coupling_sum)
    return np.concatenate([dtheta_dt, domega_dt])

def run_penetration_sweep() -> None:
    """Evaluate grid synchronization stability across solar penetration steps from 0% to 100%."""
    print("================================================================================")
    print("      BSPTCL Bihar Grid Sweep: Solar Penetration vs. Synchronization Cohesion   ")
    print("================================================================================")
    
    # 1. Base Node Configurations
    nodes_config = {
        "Patna_Digha":     {"type": "Load_GSS", "inertia": 1.0, "damping": 0.1, "power": -0.8},
        "Patna_Gaurichak": {"type": "Load_GSS", "inertia": 1.0, "damping": 0.1, "power": -0.7},
        "Patna_Fatuha":    {"type": "Load_GSS", "inertia": 1.0, "damping": 0.1, "power": -0.6},
        "Patna_Khagaul":   {"type": "Load_GSS", "inertia": 1.0, "damping": 0.1, "power": -0.5},
        "Bihta_GSS":       {"type": "Load_GSS", "inertia": 1.0, "damping": 0.1, "power": -0.6},
        "NTPC_Barh":       {"type": "Thermal_Gen", "inertia": 6.0, "damping": 0.25, "power": 3.0},
        "NTPC_Kahalgaon":  {"type": "Thermal_Gen", "inertia": 5.5, "damping": 0.25, "power": 2.5},
        "NTPC_Nabinagar":  {"type": "Thermal_Gen", "inertia": 5.5, "damping": 0.25, "power": 2.5},
        "Barauni_Thermal": {"type": "Thermal_Gen", "inertia": 4.5, "damping": 0.2, "power": 1.5},
        "Kanti_Thermal":   {"type": "Thermal_Gen", "inertia": 4.0, "damping": 0.2, "power": 1.0},
        "Kajra_Solar":     {"type": "Solar_Gen", "inertia": 0.5, "damping": 0.05, "power": 1.0},
        "Banka_Solar":     {"type": "Solar_Gen", "inertia": 0.5, "damping": 0.05, "power": 0.8},
        "Jamui_Solar":     {"type": "Solar_Gen", "inertia": 0.5, "damping": 0.05, "power": 0.7},
        "Bodhgaya_GSS":    {"type": "Load_GSS", "inertia": 1.0, "damping": 0.1, "power": -0.8},
        "Gaya_GSS":        {"type": "Load_GSS", "inertia": 1.0, "damping": 0.1, "power": -0.7},
        "Dehri_GSS":       {"type": "Load_GSS", "inertia": 1.0, "damping": 0.1, "power": -0.8},
        "Biharsharif_GSS": {"type": "Load_GSS", "inertia": 1.2, "damping": 0.12, "power": -1.2},
        "Begusarai_GSS":   {"type": "Load_GSS", "inertia": 1.0, "damping": 0.1, "power": -0.7},
        "Bhagalpur_GSS":   {"type": "Load_GSS", "inertia": 1.0, "damping": 0.1, "power": -0.8},
        "Darbhanga_GSS":   {"type": "Load_GSS", "inertia": 1.0, "damping": 0.1, "power": -0.9},
        "Samastipur_GSS":  {"type": "Load_GSS", "inertia": 1.0, "damping": 0.1, "power": -0.7},
        "Kishanganj_GSS":  {"type": "Load_GSS", "inertia": 1.0, "damping": 0.1, "power": -0.6},
        "Gopalganj_GSS":   {"type": "Load_GSS", "inertia": 1.0, "damping": 0.1, "power": -0.5},
        "Nawada_GSS":      {"type": "Load_GSS", "inertia": 1.0, "damping": 0.1, "power": -0.5}
    }
    
    node_names = list(nodes_config.keys())
    num_nodes = len(node_names)
    
    # Load K matrix from CSV
    edges_csv = DATA_PROCESSED / "bihar_grid_edges.csv"
    if not edges_csv.exists():
        raise FileNotFoundError(f"Missing topology parameters CSV: {edges_csv}. Please run grid_laplacian.py first!")
    df_edges = pd.read_csv(edges_csv)
    
    K_matrix = np.zeros((num_nodes, num_nodes))
    for _, row in df_edges.iterrows():
        u = row["From"]
        v = row["To"]
        k_val = row["Coupling_K"]
        idx_u = node_names.index(u)
        idx_v = node_names.index(v)
        
        # Using a baseline coupling multiplier of 4.0
        K_matrix[idx_u, idx_v] = 4.0 * k_val
        K_matrix[idx_v, idx_u] = 4.0 * k_val
        
    # Order of generators we replace to increase solar penetration
    generators_to_replace = ["Kanti_Thermal", "Barauni_Thermal", "NTPC_Nabinagar", "NTPC_Kahalgaon", "NTPC_Barh"]
    
    sweep_results = []
    
    # We define 6 steps:
    # Step 0: 0% solar (all generators are high-inertia thermal)
    # Step 1: 20% solar (Kanti thermal replaced by solar)
    # Step 2: 40% solar (Kanti + Barauni replaced by solar)
    # Step 3: 60% solar (Kanti + Barauni + Nabinagar replaced by solar)
    # Step 4: 80% solar (Kanti + Barauni + Nabinagar + Kahalgaon replaced by solar)
    # Step 5: 100% solar (all 5 major generators replaced by low-inertia solar)
    
    for step in range(6):
        penetration_pct = step * 20
        print(f"\nEvaluating Solar Penetration Step {step}: {penetration_pct}% Solar...")
        
        # Copy base parameters
        M = np.array([nodes_config[name]["inertia"] for name in node_names])
        D = np.array([nodes_config[name]["damping"] for name in node_names])
        P = np.array([nodes_config[name]["power"] for name in node_names])
        
        # Replace generators based on current penetration step
        for idx_rep in range(step):
            name_rep = generators_to_replace[idx_rep]
            node_idx = node_names.index(name_rep)
            
            # Reduce inertia representing inverter-coupled solar (low virtual inertia)
            M[node_idx] = 0.4    # High-inertia turbine (4.0-6.0) -> Low-inertia solar (0.4)
            D[node_idx] = 0.05   # Lower damping
            
        # Balance Power Sum
        load_indices = [i for i, name in enumerate(node_names) if nodes_config[name]["type"] == "Load_GSS"]
        gen_power = np.sum([P[i] for i in range(num_nodes) if P[i] > 0])
        total_baseline_load = np.sum([P[i] for i in load_indices])
        scale_factor = -gen_power / total_baseline_load
        for idx in load_indices:
            P[idx] *= scale_factor
            
        # Define Initial Conditions with identical random seed for comparison consistency
        np.random.seed(24)
        theta_0 = np.random.uniform(-0.25, 0.25, num_nodes)
        omega_0 = np.zeros(num_nodes)
        y0 = np.concatenate([theta_0, omega_0])
        
        # Simulate swing equations for 30 seconds
        t_span = (0, 30)
        t_eval = np.linspace(0, 30, 600)
        sol = solve_ivp(
            swing_ode,
            t_span,
            y0,
            args=(num_nodes, M, D, P, K_matrix),
            t_eval=t_eval,
            method='RK45'
        )
        
        phases = sol.y[:num_nodes, :]
        frequencies = sol.y[num_nodes:, :]
        
        # Compute steady-state Kuramoto order parameter R at t_end
        phases_end = phases[:, -1]
        complex_sum = np.sum(np.exp(1j * phases_end))
        r_steady = np.abs(complex_sum) / num_nodes
        
        # Max frequency deviation at steady state
        max_freq_dev_end = np.max(np.abs(frequencies[:, -1]))
        
        print(f"  Steady-State Order Parameter R = {r_steady:.4f}")
        print(f"  Max Frequency Deviation        = {max_freq_dev_end:.6f} rad/s")
        
        sweep_results.append({
            "Penetration_pct": penetration_pct,
            "OrderParameter_R": r_steady,
            "MaxFreqDeviation": max_freq_dev_end
        })
        
    df_results = pd.DataFrame(sweep_results)
    
    # Save sweep results to CSV
    results_csv = DATA_PROCESSED / "solar_penetration_sweep_results.csv"
    df_results.to_csv(results_csv, index=False)
    print(f"\nSaved sweep results parameters to: {results_csv}")
    
    # Plotting Sweep Curves
    print("\nPlotting Solar Penetration Sweep analysis...")
    fig, ax1 = plt.subplots(figsize=(10, 6))
    
    color = 'tab:blue'
    ax1.set_xlabel('Solar Penetration Ratio in Bihar Grid (%)', fontweight='bold')
    ax1.set_ylabel('Steady-State Grid Cohesion R(t_end)', color=color, fontweight='bold')
    ax1.plot(df_results['Penetration_pct'], df_results['OrderParameter_R'], color=color, marker='o', linewidth=2.5, markersize=8)
    ax1.tick_params(axis='y', labelcolor=color)
    ax1.set_ylim([0.4, 1.05])
    ax1.grid(True, linestyle='--', alpha=0.5)
    
    ax2 = ax1.twinx()  
    color = 'tab:red'
    ax2.set_ylabel('Max Steady-State Frequency Deviation (rad/s)', color=color, fontweight='bold')
    ax2.plot(df_results['Penetration_pct'], df_results['MaxFreqDeviation'], color=color, marker='s', linestyle='--', linewidth=2.0, markersize=8)
    ax2.tick_params(axis='y', labelcolor=color)
    
    plt.title("Impact of Solar Penetration on Bihar Transmission Grid Synchronization Stability", fontsize=12, fontweight='bold')
    fig.tight_layout()  
    
    fig_file_png = OUTPUTS_FIGURES / "solar_penetration_sweep.png"
    fig_file_pdf = OUTPUTS_FIGURES / "solar_penetration_sweep.pdf"
    plt.savefig(fig_file_png, dpi=300)
    plt.savefig(fig_file_pdf, bbox_inches='tight')

    fig_paper2_png = LATEX_PAPER2_FIGURES / "solar_penetration_sweep.png"
    fig_paper2_pdf = LATEX_PAPER2_FIGURES / "solar_penetration_sweep.pdf"
    plt.savefig(fig_paper2_png, dpi=300)
    plt.savefig(fig_paper2_pdf, bbox_inches='tight')
    plt.close()
    
    print(f"Successfully saved sweep plot at: {fig_file_png} and {fig_paper2_png}")

if __name__ == "__main__":
    run_penetration_sweep()
