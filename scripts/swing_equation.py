from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_PROCESSED = REPO_ROOT / "data"
OUTPUTS_FIGURES = REPO_ROOT / "figures"

DATA_PROCESSED.mkdir(parents=True, exist_ok=True)
OUTPUTS_FIGURES.mkdir(parents=True, exist_ok=True)

def swing_ode(t: float, y: np.ndarray, num_nodes: int, M: np.ndarray, D: np.ndarray, P: np.ndarray, K_matrix: np.ndarray) -> np.ndarray:
    """Coupled Kuramoto-swing ODE system computing phase angle and frequency derivatives."""
    # State variables: y is a flat vector of length 2*num_nodes
    # First half are phase angles (theta), second half are frequency deviations (omega)
    theta = y[:num_nodes]
    omega = y[num_nodes:]
    
    dtheta_dt = omega
    domega_dt = np.zeros(num_nodes)
    
    for i in range(num_nodes):
        # Coupling term sum_j K_ij * sin(theta_i - theta_j)
        coupling_sum = 0.0
        for j in range(num_nodes):
            if K_matrix[i, j] > 0:
                coupling_sum += K_matrix[i, j] * np.sin(theta[i] - theta[j])
                
        # Swing Equation: d_omega_i / dt = (1/M_i) * (P_i - D_i * omega_i - coupling_sum)
        domega_dt[i] = (1.0 / M[i]) * (P[i] - D[i] * omega[i] - coupling_sum)
        
    return np.concatenate([dtheta_dt, domega_dt])

def simulate_grid_synchronization() -> None:
    """Simulate coupled swing dynamics on the BSPTCL grid and export frequency/phase plots."""
    print("================================================================================")
    print("        BSPTCL Bihar Transmission Grid: Coupled Swing Dynamics Simulation      ")
    print("================================================================================")
    
    # 1. Load Node Configuration
    # We define the 24 nodes to match the Laplacian script
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
    
    # 2. Extract vectors of parameters
    M = np.array([nodes_config[name]["inertia"] for name in node_names])
    D = np.array([nodes_config[name]["damping"] for name in node_names])
    P = np.array([nodes_config[name]["power"] for name in node_names])
    
    # Crucial physical check: Power Balance sum P_i = 0
    # Let's adjust the total generation or load exactly to achieve absolute balance
    p_sum = np.sum(P)
    print(f"Initial Power Sum: {p_sum:.3f}")
    if not np.isclose(p_sum, 0.0):
        # Adjust load nodes proportionally to balance out the active generation
        load_indices = [i for i, name in enumerate(node_names) if nodes_config[name]["type"] == "Load_GSS"]
        gen_power = np.sum([P[i] for i in range(num_nodes) if P[i] > 0])
        total_baseline_load = np.sum([P[i] for i in load_indices])
        
        # Scale loads
        scale_factor = -gen_power / total_baseline_load
        for idx in load_indices:
            P[idx] *= scale_factor
            
    print(f"Balanced Power Sum: {np.sum(P):.6f} (Generation balanced exactly with consumption)")
    
    # 3. Read topological coupling parameters from generated CSV
    edges_csv = DATA_PROCESSED / "bihar_grid_edges.csv"
    if not edges_csv.exists():
        raise FileNotFoundError(f"Missing topology parameters CSV: {edges_csv}. Please run grid_laplacian.py first!")
        
    df_edges = pd.read_csv(edges_csv)
    
    # Construct K matrix
    K_matrix = np.zeros((num_nodes, num_nodes))
    for _, row in df_edges.iterrows():
        u = row["From"]
        v = row["To"]
        k_val = row["Coupling_K"]
        
        idx_u = node_names.index(u)
        idx_v = node_names.index(v)
        
        # Scale coupling parameter to safe operating margins
        K_matrix[idx_u, idx_v] = 4.0 * k_val
        K_matrix[idx_v, idx_u] = 4.0 * k_val
        
    # 4. Define Initial Conditions
    np.random.seed(24)  # For reproducible phase perturbation
    
    # Phase angles initialized near 0 with small random perturbations (within [-0.25, 0.25] radians)
    theta_0 = np.random.uniform(-0.25, 0.25, num_nodes)
    
    # Average frequency is perfectly synchronous (50Hz deviation is 0.0)
    omega_0 = np.zeros(num_nodes)
    
    y0 = np.concatenate([theta_0, omega_0])
    
    # 5. Integrate coupled 2nd-order ODE swing system
    # Time span: 40 seconds to show dynamic transition and complete locking
    t_span = (0, 40)
    t_eval = np.linspace(0, 40, 1000)
    
    print("\nIntegrating coupled swing equations over 40-second horizon...")
    sol = solve_ivp(
        swing_ode,
        t_span,
        y0,
        args=(num_nodes, M, D, P, K_matrix),
        t_eval=t_eval,
        method='RK45'
    )
    
    # 6. Post-processing & Verification
    phases = sol.y[:num_nodes, :]
    frequencies = sol.y[num_nodes:, :]
    
    # Calculate phase cohesion parameter (Kuramoto order parameter R)
    # R(t) = (1/N) * |sum_j e^(i * theta_j(t))|
    order_parameter = []
    for t_idx in range(len(sol.t)):
        phases_t = phases[:, t_idx]
        complex_sum = np.sum(np.exp(1j * phases_t))
        r = np.abs(complex_sum) / num_nodes
        order_parameter.append(r)
        
    print(f"\nFinal State at t = 40.0s:")
    print(f"  Kuramoto Order Parameter (R) = {order_parameter[-1]:.4f} (R -> 1 implies complete synchronization)")
    print(f"  Maximum Frequency Deviation  = {np.max(np.abs(frequencies[:, -1])):.6f} rad/s")
    
    # 7. Generate Visualization
    print("\nPlotting coupled swing frequency and phase locking curves...")
    fig, axes = plt.subplots(3, 1, figsize=(12, 12))
    
    # Panel 1: Frequency Deviations
    for i, name in enumerate(node_names):
        color = nodes_config[name]["type"]
        if color == "Thermal_Gen":
            axes[0].plot(sol.t, frequencies[i, :], 'r-', alpha=0.6)
        elif color == "Solar_Gen":
            axes[0].plot(sol.t, frequencies[i, :], 'g-', alpha=0.7)
        else:
            axes[0].plot(sol.t, frequencies[i, :], 'b-', alpha=0.3)
            
    axes[0].set_title("Frequency Deviations $\\omega_i(t)$ Across GSS Nodes", fontsize=12, fontweight='bold')
    axes[0].set_xlabel("Time (s)")
    axes[0].set_ylabel("Frequency Deviation (rad/s)")
    axes[0].grid(True)
    
    # Panel 2: Relative Phase Angles
    # We reference the phase relative to the average grid phase to remove global drift
    mean_phase = np.mean(phases, axis=0)
    for i, name in enumerate(node_names):
        color = nodes_config[name]["type"]
        relative_phase = phases[i, :] - mean_phase
        if color == "Thermal_Gen":
            axes[1].plot(sol.t, relative_phase, 'r-', alpha=0.6)
        elif color == "Solar_Gen":
            axes[1].plot(sol.t, relative_phase, 'g-', alpha=0.7)
        else:
            axes[1].plot(sol.t, relative_phase, 'b-', alpha=0.3)
            
    axes[1].set_title("Relative Phase Angles $\\theta_i(t) - \\bar{\\theta}(t)$ showing Phase-Locking", fontsize=12, fontweight='bold')
    axes[1].set_xlabel("Time (s)")
    axes[1].set_ylabel("Relative Phase (rad)")
    axes[1].grid(True)
    
    # Panel 3: Order Parameter
    axes[2].plot(sol.t, order_parameter, 'k-', linewidth=2.0)
    axes[2].set_title("Kuramoto Grid Synchronization Order Parameter $R(t)$", fontsize=12, fontweight='bold')
    axes[2].set_xlabel("Time (s)")
    axes[2].set_ylabel("Synchronization Cohesion $R(t)$")
    axes[2].set_ylim([0.8, 1.02])
    axes[2].grid(True)
    
    # Add custom legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color='red', label='NTPC Thermal Generators (High Inertia)'),
        Line2D([0], [0], color='green', label='Solar Renewable Nodes (Low Inertia)'),
        Line2D([0], [0], color='blue', alpha=0.5, label='Grid Substations (Load GSS)')
    ]
    axes[0].legend(handles=legend_elements, loc='upper right')
    
    plt.tight_layout()
    
    # Save figure
    fig_file_png = OUTPUTS_FIGURES / "bihar_grid_sync_simulation.png"
    fig_file_pdf = OUTPUTS_FIGURES / "bihar_grid_sync_simulation.pdf"
    plt.savefig(fig_file_png, dpi=300)
    plt.savefig(fig_file_pdf, bbox_inches='tight')

    plt.close()
    
    print(f"\nSuccessfully generated and saved simulation curves at: {fig_file_png}")

if __name__ == "__main__":
    simulate_grid_synchronization()
