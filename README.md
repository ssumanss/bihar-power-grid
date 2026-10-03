# Bihar Power Grid: Transmission Topology and Vulnerability Dataset

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Format: CSV](https://img.shields.io/badge/Format-CSV-blue.svg)](data/)

Open-access research dataset for the research paper:

> **"Synchronisation and Vulnerability Analysis of the Bihar Power Grid under Solar Integration"**  
> *Raj Kumar and Sandeep Suman*  
> University Department of Mathematics, Tilka Manjhi Bhagalpur University, Bhagalpur 812007, Bihar, India.  
> *Manuscript submitted for publication (2026).*

---

## Dataset Overview

This repository provides the essential open datasets required for evaluating the network topology, Laplacian spectrum, solar penetration dynamics, and $N-1$ line outage vulnerabilities of the high-voltage transmission system of the Bihar State Power Transmission Company Limited (BSPTCL).

```
bihar-power-grid/
├── README.md                                    # Dataset documentation
├── LICENSE                                      # MIT Open Source License
├── .gitignore
└── data/                                        # Essential research datasets
    ├── bihar_grid_edges.csv                     # 24-bus, 31-corridor reactances, distances & susceptances
    ├── grid_vulnerability_index.csv             # Full N-1 contingency sweep & lambda_2 degradation rankings
    ├── solar_penetration_sweep_results.csv      # 0% to 100% solar penetration transient frequency data
    └── grid_synchronization_threshold_sweep.csv # Coupling sensitivity invariance across kappa in [2.0, 6.0]
```

---

## Data Files & Field Descriptions

### 1. `data/bihar_grid_edges.csv`
Geographic and physical corridor parameters for the 24-node, 31-line transmission network:
- **`From`**: Origin substation / generation hub name
- **`To`**: Destination substation / generation hub name
- **`Distance_km`**: Geodesic distance in kilometers calculated via Haversine formula from GPS coordinates
- **`Voltage_kV`**: Operating transmission line voltage ($132\text{ kV}$ or $220\text{ kV}$)
- **`Conductor`**: ACSR conductor type (`Zebra` for $220\text{ kV}$, `Panther` for $132\text{ kV}$)
- **`Reactance_Ohm`**: Line reactance $X_{ij} = X_0 L_{ij}$ with $X_0 \approx 0.32\,\Omega/\text{km}$
- **`Coupling_K`**: Dimensionless coupling susceptance $K_{ij} = \gamma_0 / X_{ij}$ normalized with $\gamma_0 = 50.0\,\Omega$

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
Dynamic coupling multiplier sensitivity sweep:
- Confirms the exact invariance of the Nawada--Biharsharif bottleneck vulnerability drop ($65.44\%$) across scaling factors $\kappa \in [2.0, 6.0]$ under LAPACK re-diagonalization.

---

## Citation

This repository accompanies the manuscript:

> Raj Kumar and Sandeep Suman, *"Synchronisation and Vulnerability Analysis of the Bihar Power Grid under Solar Integration"*, submitted for publication (2026).

If using this dataset in research prior to journal publication, please cite:

```bibtex
@misc{bihar_power_grid_2026,
  author       = {Kumar, Raj and Suman, Sandeep},
  title        = {Bihar Power Grid: Transmission Topology and Vulnerability Dataset},
  year         = {2026},
  howpublished = {\url{https://github.com/ssumanss/bihar-power-grid}},
  note         = {Data repository accompanying manuscript submitted for publication}
}
```

*(Full citation with journal DOI and volume details will be updated here upon publication).*

---

## License

This dataset is distributed under the open-source **MIT License** — see [LICENSE](LICENSE) for details.
