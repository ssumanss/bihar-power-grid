from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_PROCESSED = REPO_ROOT / "data"
OUTPUTS_FIGURES = REPO_ROOT / "figures"
LATEX_PAPER2_FIGURES = REPO_ROOT / "figures"
LATEX_THESIS_FIGURES = REPO_ROOT / "figures"

DATA_PROCESSED.mkdir(parents=True, exist_ok=True)
OUTPUTS_FIGURES.mkdir(parents=True, exist_ok=True)
LATEX_PAPER2_FIGURES.mkdir(parents=True, exist_ok=True)
LATEX_THESIS_FIGURES.mkdir(parents=True, exist_ok=True)


def swing_ode(
    t: float,
    y: np.ndarray,
    num_nodes: int,
    M: np.ndarray,
    D: np.ndarray,
    P: np.ndarray,
    K_matrix: np.ndarray,
    coupling_multiplier: float
) -> np.ndarray:
    """Coupled Kuramoto-swing ODE system computing phase angle and frequency derivatives.

    Parameters
    ----------
    t : float
        Current simulation time.
    y : np.ndarray
        Flat state vector of length 2 * num_nodes (first half phases theta, second half frequencies omega).
    num_nodes : int
        Number of grid buses.
    M : np.ndarray
        Inertia coefficients vector.
    D : np.ndarray
        Damping coefficients vector.
    P : np.ndarray
        Net active power injections vector.
    K_matrix : np.ndarray
        Base susceptance / coupling matrix.
    coupling_multiplier : float
        Uniform coupling scale factor sigma.

    Returns
    -------
    np.ndarray
        Concatenated time derivatives [dtheta/dt, domega/dt].
    """
    theta = y[:num_nodes]
    omega = y[num_nodes:]
    
    dtheta_dt = omega
    domega_dt = np.zeros(num_nodes)
    
    # Scale the coupling matrix by the multiplier
    scaled_K = coupling_multiplier * K_matrix
    
    for i in range(num_nodes):
        # Coupling term sum_j K_ij * sin(theta_i - theta_j)
        coupling_sum = 0.0
        for j in range(num_nodes):
            if scaled_K[i, j] > 0:
                coupling_sum += scaled_K[i, j] * np.sin(theta[i] - theta[j])
                
        # Swing Equation
        domega_dt[i] = (1.0 / M[i]) * (P[i] - D[i] * omega[i] - coupling_sum)
        
    return np.concatenate([dtheta_dt, domega_dt])


def evaluate_synchronization_threshold() -> pd.DataFrame:
    """Evaluate grid synchronization threshold across coupling multipliers.

    Returns
    -------
    pd.DataFrame
        DataFrame of sweep results across coupling multipliers.
    """
    print("================================================================================")
    print("    BSPTCL Grid: Numerical Sweep of Synchronization Threshold vs. Fiedler      ")
    print("================================================================================")
    
    # 1. Load Node Configuration
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
    
    # Extract parameter vectors
    M = np.array([nodes_config[name]["inertia"] for name in node_names])
    D = np.array([nodes_config[name]["damping"] for name in node_names])
    P = np.array([nodes_config[name]["power"] for name in node_names])
    
    # Power Balance adjustment
    load_indices = [i for i, name in enumerate(node_names) if nodes_config[name]["type"] == "Load_GSS"]
    gen_power = np.sum([P[i] for i in range(num_nodes) if P[i] > 0])
    total_baseline_load = np.sum([P[i] for i in load_indices])
    scale_factor = -gen_power / total_baseline_load
    for idx in load_indices:
        P[idx] *= scale_factor
        
    print(f"Nodes: {num_nodes}, Adjusted Power Sum: {np.sum(P):.6f}")
    
    # Read edges from CSV
    edges_csv = DATA_PROCESSED / "bihar_grid_edges.csv"
    df_edges = pd.read_csv(edges_csv)
    
    # Construct base susceptance matrix K_matrix (unscaled)
    K_matrix = np.zeros((num_nodes, num_nodes))
    for _, row in df_edges.iterrows():
        u = row["From"]
        v = row["To"]
        k_val = row["Coupling_K"]
        idx_u = node_names.index(u)
        idx_v = node_names.index(v)
        K_matrix[idx_u, idx_v] = k_val
        K_matrix[idx_v, idx_u] = k_val
        
    # Calculate Fiedler Eigenvalue of the base K_matrix
    A_m = K_matrix
    D_m = np.diag(np.sum(A_m, axis=1))
    L_m = D_m - A_m
    eigenvalues = np.sort(np.linalg.eigvalsh(L_m))
    fiedler_base = eigenvalues[1]
    print(f"Base Susceptance Graph Fiedler Eigenvalue (lambda_2_base) = {fiedler_base:.6f}")
    
    # Define Initial conditions (small random perturbations)
    np.random.seed(24)
    theta_0 = np.random.uniform(-0.1, 0.1, num_nodes)
    omega_0 = np.zeros(num_nodes)
    y0 = np.concatenate([theta_0, omega_0])
    
    # Sweep coupling multiplier from 0.1 to 6.0 (total 40 steps)
    multipliers = np.linspace(0.1, 6.0, 40)
    
    sweep_results = []
    
    t_span = (0, 30)
    t_eval = [30.0]  # Only evaluate final state to save time
    
    print("\nSweeping coupling multiplier and simulating swing equations (30 seconds per step)...")
    for idx, mult in enumerate(multipliers):
        # Run integration
        sol = solve_ivp(
            swing_ode,
            t_span,
            y0,
            args=(num_nodes, M, D, P, K_matrix, mult),
            t_eval=t_eval,
            method='RK45'
        )
        
        final_y = sol.y[:, -1]
        final_theta = final_y[:num_nodes]
        final_omega = final_y[num_nodes:]
        
        # Calculate Order Parameter R
        complex_sum = np.sum(np.exp(1j * final_theta))
        r = np.abs(complex_sum) / num_nodes
        
        # Calculate frequency standard deviation (measure of desynchronization)
        freq_std = np.std(final_omega)
        
        # Calculate Fiedler eigenvalue for this coupling strength
        fiedler_scaled = mult * fiedler_base
        
        status = "Locked" if (r > 0.95 and freq_std < 0.05) else "Drifting"
        
        sweep_results.append({
            "Coupling_Multiplier": mult,
            "Fiedler_lambda_2": fiedler_scaled,
            "Order_Parameter_R": r,
            "Freq_Std_Dev": freq_std,
            "Status": status
        })
        
        if idx % 5 == 0 or idx == len(multipliers)-1:
            print(f"  Step {idx:2d}/40: Multiplier = {mult:.2f} | Scaled lambda_2 = {fiedler_scaled:.4f} | R = {r:.4f} | Freq Std = {freq_std:.6f} | {status}")
            
    df_sweep = pd.DataFrame(sweep_results)
    
    # Save to CSV
    sweep_csv = DATA_PROCESSED / "grid_synchronization_threshold_sweep.csv"
    df_sweep.to_csv(sweep_csv, index=False)
    print(f"\nSaved synchronization threshold results to: {sweep_csv}")
    
    # Determine the critical coupling multiplier
    locked_df = df_sweep[df_sweep['Status'] == 'Locked']
    if len(locked_df) > 0:
        critical_mult = locked_df['Coupling_Multiplier'].min()
        critical_fiedler = locked_df['Fiedler_lambda_2'].min()
        print("\nCRITICAL COUPLING THRESHOLD IDENTIFIED:")
        print(f"  Critical Coupling Multiplier (sigma_c) = {critical_mult:.4f}")
        print(f"  Critical Fiedler Connectivity Margin   = {critical_fiedler:.4f}")
    else:
        print("\nWARNING: No synchronization achieved in the swept range!")
        critical_mult = None
        
    # Generate Phase-Locking Transition Plot
    print("\nGenerating phase-locking transition plot...")
    plt.figure(figsize=(10, 6))
    
    plt.plot(df_sweep['Coupling_Multiplier'], df_sweep['Order_Parameter_R'], 'b-o', label='Order Parameter $R$ (Cohesion)', linewidth=2.0)
    plt.plot(df_sweep['Coupling_Multiplier'], df_sweep['Freq_Std_Dev'], 'r-s', label='Frequency Dispersion $\\sigma_\\omega$ (rad/s)', linewidth=1.5)
    
    if critical_mult is not None:
        plt.axvline(x=critical_mult, color='darkgreen', linestyle='--', linewidth=1.5, 
                    label=f'Critical Threshold $\\sigma_c = {critical_mult:.2f}$')
        plt.fill_between(df_sweep['Coupling_Multiplier'], 0, 1, 
                         where=(df_sweep['Coupling_Multiplier'] < critical_mult), 
                         color='red', alpha=0.08, label='Desynchronized Drift Regime')
        plt.fill_between(df_sweep['Coupling_Multiplier'], 0, 1, 
                         where=(df_sweep['Coupling_Multiplier'] >= critical_mult), 
                         color='green', alpha=0.08, label='Synchronized Locking Regime')
                         
    plt.title("Synchronization Bifurcation in Coupled swing equations (BSPTCL 24-Node Grid)", fontsize=12, fontweight='bold')
    plt.xlabel("Coupling Strength Scale Multiplier ($\\sigma$)", fontsize=11, fontweight='bold')
    plt.ylabel("Synchronization Performance Metrics", fontsize=11, fontweight='bold')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.ylim([-0.05, 1.05])
    plt.legend(loc="center right")
    
    plt.tight_layout()
    
    fig_file = OUTPUTS_FIGURES / "grid_synchronization_bifurcation.png"
    plt.savefig(fig_file, dpi=300)
    fig_file_thesis = LATEX_THESIS_FIGURES / "grid_synchronization_bifurcation.png"
    plt.savefig(fig_file_thesis, dpi=300)
    plt.close()
    
    print(f"Successfully saved synchronization bifurcation plot at: {fig_file} and {fig_file_thesis}")
    return df_sweep


if __name__ == "__main__":
    evaluate_synchronization_threshold()
