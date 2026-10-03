"""Reproducibility tests for Bihar Power Grid numerical simulations.

Verifies grid topology integrity, Laplacian spectral properties, and coupling
sensitivity invariance under LAPACK re-diagonalization.
"""

from pathlib import Path
import sys
import networkx as nx
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = REPO_ROOT / "data"
SCRIPTS_DIR = REPO_ROOT / "scripts"


def test_bihar_grid_edges_exists() -> None:
    """Verify that bihar_grid_edges.csv exists and contains exactly 31 transmission corridors."""
    csv_file = DATA_DIR / "bihar_grid_edges.csv"
    assert csv_file.exists(), f"Missing {csv_file}"
    df = pd.read_csv(csv_file)
    assert len(df) == 31, f"Expected 31 lines, got {len(df)}"


def test_vulnerability_invariance_across_kappa() -> None:
    """Verify algebraic connectivity vulnerability drop invariance across coupling scale kappa.

    Re-diagonalizes the scaled Laplacian matrices via LAPACK (np.linalg.eigvalsh) at each kappa
    to verify numerical stability without drift or mode-crossing.
    """
    csv_file = DATA_DIR / "bihar_grid_edges.csv"
    df = pd.read_csv(csv_file)
    G = nx.Graph()
    for _, r in df.iterrows():
        G.add_edge(r['From'], r['To'], weight=r['Coupling_K'])
    nodes = list(G.nodes())
    A = nx.adjacency_matrix(G, nodelist=nodes, weight='weight').toarray()
    D = np.diag(np.sum(A, axis=1))
    L = D - A

    G_out = G.copy()
    G_out.remove_edge("Nawada_GSS", "Biharsharif_GSS")
    A_out = nx.adjacency_matrix(G_out, nodelist=nodes, weight='weight').toarray()
    D_out = np.diag(np.sum(A_out, axis=1))
    L_out = D_out - A_out

    for kappa in [2.0, 3.0, 4.0, 5.0, 6.0]:
        l2_b = np.sort(np.linalg.eigvalsh(kappa * L))[1]
        l2_o = np.sort(np.linalg.eigvalsh(kappa * L_out))[1]
        vuln = (l2_b - l2_o) / l2_b
        assert np.isclose(vuln, 0.6544, atol=1e-3), f"Failed for kappa={kappa}: {vuln}"


def test_sweep_coupling_sensitivity_function() -> None:
    """Verify sweep_coupling_sensitivity function returns consistent 65.44% vulnerability drop."""
    if str(SCRIPTS_DIR) not in sys.path:
        sys.path.insert(0, str(SCRIPTS_DIR))
    from vulnerability_analysis import sweep_coupling_sensitivity
    results = sweep_coupling_sensitivity()
    assert len(results) == 5
    for r in results:
        assert np.isclose(r["vulnerability_index"], 0.6544, atol=1e-3)
