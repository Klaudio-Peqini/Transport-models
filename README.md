# hospital-wwtp

A Python package for **dynamic simulation of hospital wastewater treatment** using an explicitly staged process train:

> **Influent → Equalization (EQ) → Membrane Bioreactor (MBR) → Advanced Oxidation Process (AOP) → Granular Activated Carbon (GAC)**

This package is a research-friendly and GitHub-ready rewrite of the earlier MATLAB prototype, with several important realism upgrades:

- explicit **EQ tank dynamics** instead of a conceptual buffer only
- **non-sinusoidal hospital-style influent** by default, with daily structure, random variability, and optional shock loads
- richer multi-component state vector inspired by the parameter sheet you provided
- **split COD** into biodegradable and inert fractions
- explicit **NH4 → NO3** conversion and simple denitrification
- simple **inhibition law** from disinfectants and toxic compounds
- empirical **membrane fouling indicator**
- **species-specific AOP and GAC** behaviour
- terminal control through **argparse / CLI**

The model is still intentionally mid-complexity. It is designed for:

- thesis and teaching work
- engineering sensitivity studies
- figure generation for presentations and reports
- GitHub publication and extension

It is **not** intended to replace a full ASM/CFD/industrial design simulator.

---

## 1. Why this version is more realistic than the preliminary one

The original prototype was useful, but it had several simplifications that could produce visually unrealistic behaviour, especially smooth sinusoidal time series. This version addresses the main issues:

### Influent realism
Instead of a pure sinusoid, the default influent uses:
- a hospital-like daily schedule
- meal-related organic peaks
- cleaning/disinfection peaks
- random multiplicative variability
- optional shock events for toxicants / pharmaceuticals

### Process realism
This package now includes:
- **EQ** as a real dynamic mixing volume
- **MBR** with reaction + transport + inhibition + fouling
- **AOP** with species-specific pseudo-first-order oxidation
- **GAC** with species-specific sigmoidal breakthrough behaviour

### State realism
The model now tracks or derives the following quantities:

- BOD
- biodegradable COD
- inert COD
- total COD (derived)
- TOC
- NH4
- NO3
- PO4
- suspended solids (SS)
- alkalinity
- conductivity
- phenol
- formaldehyde
- glutaraldehyde
- anionic detergent
- oils / fats
- Pb
- carbamazepine (CBZ)
- diclofenac (DCF)
- pH (derived)

These additions were chosen to better reflect the ranges you uploaded, including BOD, COD, pH, phosphate, nitrate, NH4, suspended solids, conductivity, TOC, phenols, formaldehyde, glutaraldehyde, detergents, oils, and lead.

---

## 2. Mathematical model

### 2.1 EQ tank
For every state variable $(C_i)$:

$rac{dC_{EQ,i}}{dt} = rac{Q}{V_{EQ}} \left( C_{in,i} - C_{EQ,i} \right)$

This damps short-timescale fluctuations before they reach the bioreactor.

### 2.2 MBR block
For the MBR, the general form is:

$rac{dC_{MBR,i}}{dt} = rac{Q}{V_{MBR}} \left( C_{EQ,i} - C_{MBR,i} \right) - r_i$

where reaction terms $(r_i)$ are empirical and species-specific.

#### Organics
- BOD and biodegradable COD decay faster than inert COD.
- TOC removal is moderate.
- SS and oils are reduced more strongly due to membrane retention / empirical removal.

#### Nitrogen
- NH4 is converted to NO3 by nitrification.
- NO3 is partly consumed by denitrification.
- nitrification is pH-sensitive.

#### Phosphorus
- PO4 is removed by a first-order term plus a coupling to SS removal.

#### Toxic compounds and pharmaceuticals
- phenol, formaldehyde, glutaraldehyde, detergents, Pb, CBZ and DCF all have specific apparent removal rates.

### 2.3 Inhibition law
Biological performance is reduced by disinfectants / toxicants:

$I = rac{1}{1 + a_{ph}C_{ph} + a_{fo}C_{fo} + a_{gl}C_{gl} + a_{det}C_{det} + a_{Pb}C_{Pb}}$

This captures the idea that hospital disinfectants can suppress biological treatment.

### 2.4 Fouling indicator
A simple empirical fouling state \(F\) evolves as:

$rac{dF}{dt} = a_f\,f(SS, Oils, Detergent) - b_f F$

and MBR performance is multiplied by a fouling factor:

$f_{foul} = rac{1}{1 + \gamma F}$

This is not a mechanistic membrane model, but it improves realism.

### 2.5 AOP block
The AOP block applies pseudo-first-order polishing:

\[
C_{out,i} = C_{in,i}\exp(-k_{AOP,i}	au_{AOP})
\]

with species-specific \(k_{AOP,i}\). Oxidation strength can be increased from the CLI.

### 2.6 GAC block
Selected compounds undergo breakthrough-type polishing:

\[
F_i(t)=rac{1}{1+\exp(k_{Th,i}(	au_i^* - t))}
\]

\[
C_{out,i}=C_{in,i}F_i(t)
\]

This is still a reduced adsorption model, but it is more informative than a static efficiency.

---

## 3. Repository structure

```text
hospital_wwtp_python_package_v2/
├── src/hospital_wwtp/
│   ├── __init__.py
│   ├── config.py
│   ├── influent.py
│   ├── simulation.py
│   ├── metrics.py
│   ├── plotting.py
│   ├── io_utils.py
│   └── cli.py
├── examples/
│   └── run_example.py
├── outputs/
├── pyproject.toml
├── requirements.txt
├── LICENSE
└── README.md
```

---

## 4. Installation

### Option A — run without installing the package
From the project root:

```bash
PYTHONPATH=src python examples/run_example.py --scenario nominal
```

or

```bash
PYTHONPATH=src python -m hospital_wwtp.cli simulate --scenario nominal --output-dir outputs
```

### Option B — editable install

```bash
pip install --user -e .
```

then run:

```bash
hospital-wwtp simulate --scenario nominal --output-dir outputs
```

---

## 5. CLI usage

### Single scenario

```bash
PYTHONPATH=src python -m hospital_wwtp.cli simulate   --scenario nominal   --duration-h 96   --dt-minutes 5   --flow-m3-day 500   --output-dir outputs
```

### High-load case with stronger oxidation and larger GAC bed

```bash
PYTHONPATH=src python -m hospital_wwtp.cli simulate   --scenario high   --oxidation-scale 1.3   --adsorption-scale 1.2   --V-GAC-bed-m3 8   --output-dir outputs/high_case
```

### Shock scenario

```bash
PYTHONPATH=src python -m hospital_wwtp.cli simulate   --scenario shock   --shock-hour 42   --shock-width-h 1.0   --shock-multiplier 3.5   --output-dir outputs/shock_case
```

### Run all scenarios (low / nominal / high / shock)

```bash
PYTHONPATH=src python -m hospital_wwtp.cli batch --output-dir outputs/batch
```

---

## 6. Available command-line arguments

### Global simulation control
- `--scenario {low,nominal,high,shock}`
- `--influent-mode {hospital,sinusoidal}`
- `--duration-h`
- `--dt-minutes`
- `--flow-m3-day`
- `--seed`
- `--noise-sigma`

### Shock-event control
- `--shock-hour`
- `--shock-width-h`
- `--shock-multiplier`

### Process-strength control
- `--oxidation-scale`
- `--adsorption-scale`
- `--inhibition-scale`

### Reactor sizes
- `--V-EQ-m3`
- `--V-MBR-m3`
- `--V-AOP-m3`
- `--V-GAC-bed-m3`

### Output control
- `--output-dir`
- `--no-plots`

---

## 7. Example script

A simple example is included:

```bash
PYTHONPATH=src python examples/run_example.py --scenario nominal --output-dir outputs
```

This generates all main tables and figures automatically.

---

## 8. Generated outputs

The package writes:

### Figures
- `results_trends.png`
- `micropollutants_trends.png`
- `effluent_24h_avg.png`
- `gac_breakthrough.png`
- `pH_and_fouling.png`

### Tables / CSV files
- `simulation_timeseries.csv`
- `stage_24h_average.csv`
- `removal_summary.csv`
- `compliance_summary.csv`
- `risk_summary.csv`
- `breakthrough_summary.csv`
- `config_used.json`

These are intended to be presentation-ready and report-friendly.

---

## 9. Interpreting the outputs

### `results_trends.png`
Shows how key bulk quantities evolve through the plant over time.

### `micropollutants_trends.png`
Shows behaviour of toxicants and pharmaceuticals, especially the importance of AOP + GAC polishing.

### `effluent_24h_avg.png`
Compares final 24-hour average effluent values to target values.

### `gac_breakthrough.png`
Illustrates when selected adsorbable compounds begin to appear more strongly at the final outlet.

### `pH_and_fouling.png`
Shows whether biological operation stays in a plausible pH window and whether fouling is trending upward.

---

## 10. Practical modelling notes

### What this package is good for
- scenario comparison
- sensitivity studies
- educational demonstrations
- thesis plots and GitHub publication
- testing whether EQ/MBR/AOP/GAC process choices are qualitatively consistent

### What this package is not yet
- a full activated sludge model
- a membrane-fouling design tool
- a rigorous adsorption bed simulator
- a regulatory compliance calculator with field calibration

---

## 11. Suggested next upgrades

If you continue developing the repository, the most meaningful upgrades would be:

1. calibration to measured hospital wastewater data
2. explicit biomass / oxygen states in the MBR
3. separate weekday/weekend occupancy schedules
4. compound-specific AOP chemistry based on pH / oxidant dose
5. better GAC mass-balance and capacity tracking
6. uncertainty propagation or Monte Carlo ensembles
7. comparison of multiple treatment trains, not only EQ→MBR→AOP→GAC

---

## 12. Citation / reuse note

If you upload this to GitHub, it is a good idea to add:
- `CITATION.cff`
- `.gitignore`
- a short discussion in the README explaining that this is a **mid-complexity research prototype**, not a certified design package.

---

## 13. Minimal workflow summary

```bash
unzip hospital_wwtp_python_package_v2.zip
cd hospital_wwtp_python_package_v2
pip install --user -r requirements.txt
PYTHONPATH=src python examples/run_example.py --scenario nominal --output-dir outputs
```

or, after editable install:

```bash
pip install --user -e .
hospital-wwtp simulate --scenario nominal --output-dir outputs
```

---

## 14. Scientific positioning

This repository sits between:
- a very simple first-order teaching model, and
- a full-scale mechanistic wastewater treatment simulator.

That is deliberate.

It is complex enough to produce realistic non-sinusoidal transient figures, stagewise removals, risk indicators, and sensitivity studies, while remaining simple enough to be readable, teachable, and extendable.
