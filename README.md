# Synchronisation and Vulnerability Analysis of the Bihar Power Grid under Solar Integration

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![Formal Verification: Lean 4](https://img.shields.io/badge/Lean_4-Mathlib_v4.16.0-green.svg)](proofs/Synchronization.lean)

This repository contains the complete open-access replication package, topological transmission dataset, numerical simulation suite, and formal Lean 4 verification proofs for the research paper:

> **"Synchronisation and Vulnerability Analysis of the Bihar Power Grid under Solar Integration"**  
> *Raj Kumar and Sandeep Suman*  
> University Department of Mathematics, Tilka Manjhi Bhagalpur University, Bhagalpur 812007, Bihar, India.  
> Submitted to: *Chaos: An Interdisciplinary Journal of Nonlinear Science* (AIP Publishing).

---

## Repository Structure

```
bihar-power-grid/
├── README.md                                    # Replication overview and guide
├── LICENSE                                      # MIT Open Source License
├── pyproject.toml                               # Python project and dependency specification
├── .gitignore
├── data/                                        # Processed topological and simulation datasets
│   ├── bihar_grid_edges.csv                     # 24-node, 31-corridor reactances, distances & susceptances
│   ├── grid_vulnerability_index.csv             # Full N-1 contingency sweep & lambda_2 degradation rankings
│   ├── solar_penetration_sweep_results.csv      # 0% to 100% solar penetration transient frequency data
│   └── grid_synchronization_threshold_sweep.csv # Coupling sensitivity invariance across kappa in [2.0, 6.0]
├── figures/                                     # High-resolution vector PDF and PNG publication figures
│   ├── bihar_grid_topology_tikz.pdf             # TikZ vector map of 24 substations and transmission lines
│   ├── bihar_grid_sync_simulation.pdf           # 40-second Kuramoto swing equation dynamics and order parameter
│   ├── solar_penetration_sweep.pdf              # RoCoF and frequency deviance vs solar penetration
│   └── bihar_grid_vulnerability_centrality.pdf  # Outage vulnerability vs betweenness centrality ranking
├── proofs/                                      # Interactive theorem prover verification files
│   └── Synchronization.lean                     # Lean 4 / Mathlib formal proof of Rayleigh quotient outage sensitivity
├── scripts/                                     # Standalone Python simulation and analysis suite
│   ├── grid_laplacian.py                        # Graph Laplacian construction, GPS Haversine and Fiedler value
│   ├── swing_equation.py                        # RK45 numerical integration of non-uniform swing equations
│   ├── solar_penetration_sweep.py               # Inertia decommissioning parameter sweep
│   ├── vulnerability_analysis.py                # Reactance-weighted betweenness centrality and N-1 sweep
│   └── synchronization_threshold.py             # Analytical synchronization threshold validation
└── tests/                                       # Automated verification and reproducibility tests
    └── test_reproducibility.py                  # Pytest asserting exact numerical invariance
```

---

## Key Findings & Data Description

1. **Topological Guaranteed Synchrony:**  
   The 24-node, 31-corridor 132/220 kV high-voltage transmission network of Bihar State Power Transmission Company Limited (BSPTCL) has an algebraic connectivity of $\lambda_2 = 0.264923 > 0$. By the Master Stability Framework, complete frequency synchronisation is structurally guaranteed.

2. **Solar Penetration & Transient Volatility:**  
   As conventional synchronous generation (NTPC Barh, Kahalgaon, Nabinagar, Barauni, Kanti) is displaced by zero-inertia utility solar PV (0% to 100% sweep), transient frequency deviations surge by **66.2%** to $0.555\text{ rad/s}$ ($\approx 0.088\text{ Hz}$), approaching the statutory $0.1\text{ Hz}$ safety margin established by the Central Electricity Authority (CEA) of India.

3. **Critical Non-Islanding Bottleneck:**  
   An $N-1$ single-line outage sweep combined with reactance-weighted betweenness centrality identifies the **Biharsharif GSS $\longleftrightarrow$ Nawada GSS** corridor as the state's most critical non-islanding bottleneck, whose tripping causes a **65.4% collapse** in algebraic connectivity.

4. **Formal Verification in Lean 4:**  
   The first-order Rayleigh quotient perturbation formula:
   $$\Delta\lambda_2 \approx -K_{ij} (v_{2,i} - v_{2,j})^2$$
   is formally proven in Lean 4 (Mathlib v4.16.0) with zero axioms (`sorry`-free) in [`proofs/Synchronization.lean`](proofs/Synchronization.lean).

---

## Quick Start & Reproduction

### Prerequisites
- Python $\ge$ 3.10
- Recommended: [uv](https://docs.astral.sh/uv/) (fast Python package manager)

### Installation
Clone the repository:
```bash
git clone https://github.com/ssumanss/bihar-power-grid.git
cd bihar-power-grid
```

Install dependencies:
```bash
uv pip install -e ".[dev]"
# Or via standard pip:
pip install -e ".[dev]"
```

### Running the Simulations
Generate all network datasets and figures:
```bash
# 1. Graph Laplacian spectrum and GPS topology:
uv run python scripts/grid_laplacian.py

# 2. Non-uniform Kuramoto swing equation dynamics:
uv run python scripts/swing_equation.py

# 3. Solar penetration inertia decommissioning sweep:
uv run python scripts/solar_penetration_sweep.py

# 4. Reactance-weighted betweenness and N-1 contingency sweep:
uv run python scripts/vulnerability_analysis.py
```

### Running Automated Reproducibility Tests
```bash
uv run pytest tests/ -v
```

---

## Formal Proofs in Lean 4

To verify the mathematical proofs in Lean 4:
```bash
# Requires Lean 4 and Lake
cd proofs/
lake build
```

---

## Citation

If you use this dataset, code, or formal verification proofs in your research, please cite:

```bibtex
@article{kumar2026synchronisation,
  title={Synchronisation and Vulnerability Analysis of the Bihar Power Grid under Solar Integration},
  author={Kumar, Raj and Suman, Sandeep},
  journal={Chaos: An Interdisciplinary Journal of Nonlinear Science},
  year={2026},
  publisher={AIP Publishing},
  url={https://github.com/ssumanss/bihar-power-grid}
}
```

---

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.
