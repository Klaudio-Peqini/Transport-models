# Transport Models Suite: Hospital Wastewater + Analytical Buckley–Leverett Polymer Flooding

This repository brings together **two complementary modelling layers**:

1. **`hospital_wwtp`**  
   A dynamic multi-component treatment-train simulator for hospital wastewater:
   $$\mathrm{EQ} \rightarrow \mathrm{MBR} \rightarrow \mathrm{AOP} \rightarrow \mathrm{GAC}$$
   with realistic diurnal influent structure, toxic-load inhibition, stage-specific polishing, and CLI control.

2. **`polymer_bl`**  
   A corrected and publication-oriented **analytical Buckley–Leverett implementation** for 1D polymer flooding in porous media using:
   - Corey relative permeabilities
   - a physically consistent fractional-flow formulation
   - the **method of characteristics**
   - the **Welge tangent condition** for shock / breakthrough
   - waterflood vs polymer-flood comparison
   - publication-quality plots

The point of collecting both models in one repository is not accidental. They are connected by the same underlying mathematical language:
**nonlinear transport, front propagation, breakthrough, and cumulative recovery / removal**.

---

## 1. Why this repository exists

Earlier work in this project focused on a **hospital wastewater treatment model** with explicit process units.  
The new request from the professor adds a second, more theoretical layer:

> implement the **analytical Buckley–Leverett solution** for polymer flooding in python3 and relate it conceptually to the transport ideas already used in the wastewater model.

This repository therefore serves two goals at once:

- provide an operational environmental-engineering simulator
- provide a compact nonlinear-transport reference implementation in porous media

---

## 2. Folder structure

```text
transport_models_suite/
├── src/
│   ├── hospital_wwtp/
│   │   ├── __init__.py
│   │   ├── config.py
│   │   ├── influent.py
│   │   ├── simulation.py
│   │   ├── metrics.py
│   │   ├── plotting.py
│   │   ├── io_utils.py
│   │   └── cli.py
│   └── polymer_bl/
│       ├── __init__.py
│       ├── model.py
│       ├── plotting.py
│       └── cli.py
├── examples/
│   ├── run_example.py
│   ├── run_polymer_example.py
│   └── run_transport_bridge_example.py
├── outputs/
├── outputs_polymer/
├── outputs_bridge/
├── pyproject.toml
├── requirements.txt
├── LICENSE
└── README.md
```

---

## 3. Analytical Buckley–Leverett model: what is implemented

### 3.1 Governing equation

For 1D incompressible displacement with constant total Darcy velocity $u$, no gravity, no dispersion, no adsorption, no degradation, and constant polymer-modified water viscosity, the water saturation satisfies

$$
\phi \frac{\partial S_w}{\partial t}
+
\frac{\partial}{\partial x}\Big(u f_w(S_w)\Big)=0.
$$

Here:

- $S_w$ is water saturation,
- $\phi$ is porosity,
- $u$ is Darcy velocity,
- $f_w$ is the water fractional-flow function.

### 3.2 Corey relative permeability model

The implementation uses the effective saturation

$$
S_e=\frac{S_w-S_{wi}}{1-S_{wi}-S_{or}},
\qquad 0 \le S_e \le 1
$$

and Corey laws

$$
k_{rw}=k_{rw0} S_e^{n_w}, \qquad
k_{ro}=k_{ro0}(1-S_e)^{n_o}.
$$

Mobilities are

$$
\lambda_w=\frac{k_{rw}}{\mu_w^{(*)}}, \qquad
\lambda_o=\frac{k_{ro}}{\mu_o},
$$

where $\mu_w^{(*)}$ is either the waterflood viscosity $\mu_w$ or the polymer-modified viscosity $\mu_w^p$.

Then the fractional flow is

$$
f_w(S_w)=\frac{\lambda_w}{\lambda_w+\lambda_o}.
$$

### 3.3 Method of characteristics

The nonlinear PDE is solved analytically by characteristics:

$$
\frac{x}{t}=\frac{u}{\phi}\frac{df_w}{dS_w}.
$$

Different saturation states move with different characteristic speeds.  
Because those speeds are nonlinear functions of $S_w$, a **shock front** forms.

### 3.4 Welge tangent / shock condition

The shock saturation $S_{wf}$ is identified numerically from the tangency condition

$$
\left.\frac{df_w}{dS_w}\right|_{S_{wf}} = \frac{f_w(S_{wf})-f_w(S_{wi})}{S_{wf}-S_{wi}}.
$$

This gives the dimensionless shock speed

$$
v_D^{\text{shock}} = \frac{f_w(S_{wf})-f_w(S_{wi})}{S_{wf}-S_{wi}},
$$

and the breakthrough pore volume injected (PVI)

$$
\mathrm{PVI}_{bt}=\frac{1}{v_D^{\text{shock}}}.
$$

### 3.5 Production history

The production-side water cut is obtained from the outlet saturation after breakthrough through the characteristic relation

$$
\frac{df_w}{dS_w}=\frac{1}{t_D}, \qquad t_D=\frac{ut}{\phi L}.
$$

The recovery factor is computed by numerical quadrature of the oil-cut history.

---

## 4. What has been improved?

### 4.1 Correct Corey relative permeabilities
The original code normalized Corey powers by a sum,
which is not the standard Corey model.  
This version uses:

$$
k_{rw}=k_{rw0} S_e^{n_w}, \qquad
k_{ro}=k_{ro0}(1-S_e)^{n_o}.
$$

### 4.2 Proper Welge shock criterion
The shock is not identified by `argmax(dfw * Sw / fw)` in a robust way.  
Instead, this package evaluates the **difference between derivative and secant slope** and selects the best tangency point.

### 4.3 Indexing bug removed
The original sketch used expressions like `fw[S_wf]`, where `S_wf` is a floating-point saturation value, not an array index.  
This package consistently uses interpolation in saturation space.

### 4.4 Saturation profiles built from characteristic inversion
The original loop for $S_w(x,t)$ was not fully consistent with the rarefaction–shock structure.  
The new implementation:
- constructs the rarefaction branch explicitly,
- sorts the characteristic-speed relation,
- interpolates the physically admissible profile,
- then imposes the shock connection to the initial state.

### 4.5 Recovery calculation made dimensionally consistent
The original `welge_recovery` block mixed geometric assumptions and dimensionless quantities in an unrobust way.  
This repository computes:
- water cut,
- oil cut,
- cumulative oil,
- recovery factor,
using the standard dimensionless time $t_D=\mathrm{PVI}$.

---

## 5. Publication-quality outputs produced by `polymer_bl`

Running the polymer examples generates:

- `fractional_flow_comparison.png`
- `saturation_profiles_comparison.png`
- `production_curves_comparison.png`
- `polymer_summary.csv`
- `fractional_flow_table.csv`
- `polymer_production_curves.csv`
- `polymer_profiles.csv`

These plots are intended for:
- reports,
- presentations,
- thesis chapters,
- side-by-side comparison with the hospital-wastewater breakthrough logic.

---

## 6. How the Buckley–Leverett model connects to the wastewater model

This is the most important conceptual bridge in the repository.

### 6.1 Common mathematical structure

The Buckley–Leverett equation is a nonlinear hyperbolic conservation law:

$$
\frac{\partial U}{\partial t} + \frac{\partial F(U)}{\partial x}=0.
$$

A simplified pollutant-transport equation without diffusion/reaction has the same structure:

$$
\phi\frac{\partial C}{\partial t}+u\frac{\partial C}{\partial x}=0.
$$

Once adsorption, retardation, or nonlinear flux laws are introduced, the wastewater system also becomes a nonlinear transport problem.

### 6.2 Physical analogy

| Porous-media displacement | Wastewater system |
|---|---|
| saturation front | contamination / polishing front |
| water cut at producer | breakthrough at GAC outlet |
| shock formation | sharp concentration breakthrough |
| cumulative oil recovery | cumulative contaminant removal |
| method of characteristics | transport-front interpretation |

### 6.3 Practical bridge implemented here

The script `examples/run_transport_bridge_example.py` runs:

1. the **Buckley–Leverett waterflood vs polymer-flood comparison**  
2. the **hospital wastewater model**  
3. a bridge plot, `transport_bridge.png`, showing normalized breakthrough curves on the same dimensionless axis.

This is not claiming the two systems are physically identical.  
It demonstrates that both systems share:
- front-propagation logic,
- breakthrough behaviour,
- nonlinear transport interpretation.

---

## 7. Running the code from terminal

No virtual environment is required if you use `python3PATH=src`.

### 7.1 Hospital wastewater model

```bash
python3PATH=src python3 -m hospital_wwtp.cli simulate --scenario nominal --output-dir outputs
```

### 7.2 Analytical Buckley–Leverett model

```bash
python3PATH=src python3 -m polymer_bl.cli --mu-w 1 --mu-wp 12 --mu-o 5 --output-dir outputs_polymer
```

### 7.3 Direct example script for Buckley–Leverett

```bash
python3PATH=src python3 examples/run_polymer_example.py
```

### 7.4 Direct bridge example

```bash
python3PATH=src python3 examples/run_transport_bridge_example.py
```

---

## 8. Useful command-line flags for the polymer model

```bash
python3PATH=src python3 -m polymer_bl.cli \
  --phi 0.20 \
  --L 100 \
  --u 1e-5 \
  --Swi 0.20 \
  --Sor 0.0 \
  --nw 2 --no 2 \
  --krw0 0.8 --kro0 1.0 \
  --mu-w 1.0 \
  --mu-wp 15.0 \
  --mu-o 5.0 \
  --pvi-max 2.5 \
  --profile-pvi 0.2 0.5 1.0 1.5 2.0 \
  --output-dir outputs_polymer/high_polymer_case
```

Meaning of the most important flags:

- `--mu-w`: baseline water viscosity in cP
- `--mu-wp`: polymer-modified water viscosity in cP
- `--mu-o`: oil viscosity in cP
- `--profile-pvi`: PVI values at which saturation profiles are plotted
- `--pvi-max`: maximum injected pore volumes for the production curves

---

## 9. Interpreting the outputs

### 9.1 `fractional_flow_comparison.png`
Shows how polymer thickening modifies the fractional-flow curve and shifts the Welge tangent.  
A higher polymer viscosity generally:
- reduces the mobility ratio,
- delays breakthrough,
- improves sweep,
- increases cumulative recovery.

### 9.2 `saturation_profiles_comparison.png`
Shows the saturation front in dimensionless space $x/L$.  
Compared with waterflooding, polymer flooding generally displays a more favourable displacement profile.

### 9.3 `production_curves_comparison.png`
Contains:
- **water cut vs PVI**
- **recovery factor vs PVI**

These are the clearest engineering outputs for discussing improved displacement.

### 9.4 `transport_bridge.png`
Connects the porous-media transport front to the wastewater breakthrough logic by comparing normalized front-like curves.

---

## 10. Recommended scientific workflow

A useful workflow for presentations and thesis writing is:

1. run the wastewater model for the nominal case  
2. run the polymer analytical model for baseline vs polymer  
3. generate the bridge plot  
4. discuss:
   - the same transport mathematics,
   - different physical contexts,
   - why breakthrough curves matter in both.

---

## 11. Installation

From the project root:

```bash
pip install --user -e .
```

Then you may run:

```bash
hospital-wwtp simulate --scenario nominal --output-dir outputs
polymer-bl --mu-w 1 --mu-wp 12 --mu-o 5 --output-dir outputs_polymer
```

---

## 12. Minimal dependencies

- `numpy`
- `scipy`
- `pandas`
- `matplotlib`

They are listed in `requirements.txt` and `pyproject.toml`.

---

## 13. Suggested next extensions

### For the porous-media side
- add adsorption / retardation
- add dispersion and compare analytical vs numerical fronts
- add polymer slug followed by chase water
- add inaccessible pore volume (IPV)
- add resistance factor / residual resistance factor (RF/RRF)

### For the wastewater side
- calibrate parameters against measurements
- add uncertainty propagation
- add stronger compound-specific AOP/GAC kinetics
- link pH explicitly to reaction rates
- extend compliance/risk outputs

---

## 14. Final note

This repository intentionally places **environmental process modelling** and **porous-media transport theory** in the same framework.  
That is not a software convenience alone; it is a modelling statement:

> the same mathematical structures — conservation laws, fronts, breakthrough, and cumulative response — recur across apparently different physical systems.

That is the main scientific message tying this work together.
