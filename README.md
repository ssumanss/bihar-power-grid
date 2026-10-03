# Bihar Power Grid: Synchronization, Vulnerability & Solar Integration

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Format: CSV](https://img.shields.io/badge/Format-CSV-blue.svg)](data/)
[![Lean 4: v4.16.0](https://img.shields.io/badge/Lean_4-v4.16.0-green.svg)](proofs/)
[![Tests: Pytest](https://img.shields.io/badge/Tests-Passing-brightgreen.svg)](tests/)

Open-access replication codebase, formal mathematical proofs, and research datasets accompanying the manuscript:

> **"Synchronisation and Vulnerability Analysis of the Bihar Power Grid under Solar Integration"**  
> *Raj Kumar Raj and Sandeep Suman*  
> University Department of Mathematics, Tilka Manjhi Bhagalpur University, Bhagalpur 812007, Bihar, India.  
> *Chaos: An Interdisciplinary Journal of Nonlinear Science (AIP Publishing), Submitted (2026).*

---

## Repository Overview

This repository provides full end-to-end reproducibility for the empirical transmission topology, spectral graph Laplacians, coupled Kuramoto-swing ODE simulations, $N-1$ contingency vulnerability rankings, and interactive formal proofs in Lean 4 for the 24-bus transmission network of the Bihar State Power Transmission Company Limited (BSPTCL).

```
bihar-power-grid/
├── README.md                                    # Comprehensive replication guide & documentation
├── LICENSE                                      # MIT Open Source License
├── .gitignore
├── data/                                        # Empirical and simulated research datasets
│   ├── bihar_grid_edges.csv                     # 24-bus, 31-corridor reactances, distances & susceptances
│   ├── grid_vulnerability_index.csv             # Full N-1 contingency sweep & lambda_2 degradation rankings
│   ├── solar_penetration_sweep_results.csv      # 0% to 100% solar penetration transient frequency dynamics
│   └── grid_synchronization_threshold_sweep.csv # Dynamic coupling multiplier sensitivity sweep
├── scripts/                                     # Standalone numerical simulation & analysis pipeline
│   ├── grid_laplacian.py                        # Network graph formulation, Laplacian spectrum & Fiedler vector
│   ├── swing_equation.py                        # Coupled second-order Kuramoto-swing ODE integrator (RK45)
│   ├── solar_penetration_sweep.py               # Inertia reduction sweep (0-100% solar PV) & RoCoF dynamics
│   ├── vulnerability_analysis.py                # Exhaustive N-1 line outages, islanding flags & kappa invariance
│   └── synchronization_threshold.py             # Phase-locking bifurcation & critical coupling margin
├── proofs/                                      # Formal interactive theorem proving in Lean 4
│   ├── Synchronization.lean                     # Zero-axiom formal proof of Rayleigh quotient outage sensitivity
│   ├── lakefile.lean                            # Lake package configuration (Mathlib v4.16.0)
│   └── lean-toolchain                           # Lean toolchain pin (v4.16.0)
├── tests/                                       # Automated regression and reproducibility tests
│   └── test_reproducibility.py                  # Pytest suite asserting topology, spectrum & invariance
└── figures/                                     # Publication-grade vector (PDF) and raster (PNG) figures
    ├── bihar_grid_topology.png                  # Transmission grid geographic topology map
    ├── bihar_grid_sync_simulation.pdf / .png    # Swing dynamics, phase locking & frequency trajectories
    ├── solar_penetration_sweep.pdf / .png       # Order parameter R, peak deviance & RoCoF vs solar %
    ├── bihar_grid_vulnerability_centrality.pdf  # N-1 lambda_2 drop vs reactance-weighted betweenness
    └── grid_synchronization_bifurcation.png     # Order parameter & frequency dispersion bifurcation curve
```

---

## Datasets (`data/`)

### 1. `data/bihar_grid_edges.csv`
Geographic and physical corridor parameters for the 24-node, 31-line transmission network:
- **`From`**: Origin Grid Substation (GSS) or generation hub name
- **`To`**: Destination GSS or generation hub name
- **`Distance_km`**: Geodesic line length calculated via Haversine formula from GPS coordinates
- **`Voltage_kV`**: Operating transmission voltage ($132\text{ kV}$ Panther or $220\text{ kV}$ Zebra)
- **`Conductor`**: ACSR conductor specification (`Zebra` for $220\text{ kV}$, `Panther` for $132\text{ kV}$)
- **`Reactance_Ohm`**: Line reactance $X_{ij} = X_0 L_{ij}$ ($X_0 \approx 0.32\,\Omega/\text{km}$)
- **`Coupling_K`**: Dimensionless coupling susceptance $K_{ij} = \gamma_0 / X_{ij}$ normalized to reference impedance $\gamma_0 = 50.0\,\Omega$

### 2. `data/grid_vulnerability_index.csv`
Comprehensive $N-1$ contingency sweep assessing structural network capacity under single transmission line outages:
- **`Line`**: Line identifier
- **`From_Node`**, **`To_Node`**: Endpoints of the tripped corridor
- **`Base_Lambda2`**: Baseline algebraic connectivity ($\lambda_2 = 0.264923$)
- **`Outage_Lambda2`**: Post-contingency algebraic connectivity $\lambda_{2,\text{outage}}$
- **`Vulnerability_Drop`**: Fractional algebraic connectivity drop $(\lambda_{2,\text{base}} - \lambda_{2,\text{outage}}) / \lambda_{2,\text{base}}$
- **`Islanding`**: Boolean flag (`True`/`False`) indicating topological network partitioning

### 3. `data/solar_penetration_sweep_results.csv`
System swing dynamics metrics across six discrete solar penetration stages ($0\%$ to $100\%$):
- **`Solar_Penetration_Pct`**: Percentage of generation from utility solar PV inverters ($0\%$, $20\%$, $40\%$, $60\%$, $80\%$, $100\%$)
- **`Grid_Cohesion_R`**: Steady-state Kuramoto phase order parameter $R(t_{\text{end}})$
- **`Max_Freq_Dev_rad_s`**: Peak transient frequency deviance from $50\text{ Hz}$ reference ($\text{rad/s}$)
- **`RoCoF_rad_s2`**: Peak transient Rate of Change of Frequency ($\text{rad/s}^2$)

### 4. `data/grid_synchronization_threshold_sweep.csv`
Dynamic coupling multiplier sensitivity sweep confirming the exact invariance of the Nawada--Biharsharif bottleneck vulnerability drop ($65.44\%$) across scaling factors $\kappa \in [2.0, 6.0]$ under LAPACK re-diagonalization.

---

## Python Simulation Pipeline (`scripts/`)

All scripts run standalone with standard scientific Python dependencies (`numpy`, `scipy`, `pandas`, `networkx`, `matplotlib`):

```bash
# 1. Compute Graph Laplacian, Fiedler eigenvalue (lambda_2 = 0.2649), and spectral properties
python scripts/grid_laplacian.py

# 2. Integrate coupled 2nd-order Kuramoto-swing ODEs (RK45 Dormand-Prince, 40 s horizon)
python scripts/swing_equation.py

# 3. Execute the 6-stage solar penetration inertia degradation sweep (0% to 100%)
python scripts/solar_penetration_sweep.py

# 4. Perform exhaustive N-1 line contingency analysis & coupling sensitivity invariance sweep
python scripts/vulnerability_analysis.py

# 5. Evaluate phase-locking bifurcation and identify critical coupling threshold
python scripts/synchronization_threshold.py
```

Generated publication plots are saved automatically to `figures/`.

---

## Automated Reproducibility Tests (`tests/`)

Automated tests verify topological completeness, spectral bounds, and numerical invariance:

```bash
# Run test suite with pytest
pytest tests/ -v
```

Tests assert:
1. `test_bihar_grid_edges_exists`: Exactly 31 physical transmission lines across 24 substations.
2. `test_vulnerability_invariance_across_kappa`: Nawada--Biharsharif corridor outage causes an invariant $65.44\%$ collapse in algebraic connectivity under LAPACK re-diagonalization across $\kappa \in [2.0, 6.0]$.
3. `test_sweep_coupling_sensitivity_function`: Script-level sensitivity functions match analytical expectations.

---

## Formal Verification in Lean 4 (`proofs/`)

The theoretical foundation of the paper includes a zero-axiom formal verification of the first-order Rayleigh quotient perturbation theorem:
$$\Delta \lambda_2 \approx -K_{ij} (v_i - v_j)^2$$
for arbitrary line outages in symmetric Laplacian networks.

### Building Proofs:
Prerequisites: Lean 4 (v4.16.0) and Lake (included with standard `elan` toolchain).

```bash
cd proofs
lake update
lake build
```

Theorem `fiedler_outage_sensitivity` in `proofs/Synchronization.lean` verifies cleanly with **zero `sorry` axioms**, establishing machine-checked mathematical certainty.

---

## Citation

```bibtex
@article{raj2026bihar_grid,
  author  = {Raj, Raj Kumar and Suman, Sandeep},
  title   = {Synchronisation and Vulnerability Analysis of the Bihar Power Grid under Solar Integration},
  journal = {Chaos: An Interdisciplinary Journal of Nonlinear Science},
  year    = {2026},
  note    = {Submitted. Open data and replication repository: \url{https://github.com/ssumanss/bihar-power-grid}}
}
```

---

## License

This repository is distributed under the open-source **MIT License** — see [LICENSE](LICENSE) for details.
