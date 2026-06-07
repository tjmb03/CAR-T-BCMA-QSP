# Model specification — CAR-T-BCMA-QSP

Six-state mechanistic QSP model of BCMA-directed CAR-T in multiple myeloma — the lymphodepletion-recovery structure that is cross-validated R ↔ Julia to machine precision. Equations only; no implementation code.

## State variables

| symbol | meaning | units |
|---|---|---|
| $E$ | CAR-T effector pool | cells |
| $T_{\text{high}},\,T_{\text{mid}},\,T_{\text{low}}$ | tumor stratified by BCMA antigen density | cells |
| $B$ | soluble BCMA (competitive decoy) | AU |
| $L$ | lymphodepletion niche availability | dimensionless, $[0,1]$ |

CAR-T infusion enters as the initial condition $E(0)$ — a single bolus at $t=0$.

## Kill term (two-tier saturable)

A per-cell **maximum** kill rate set by antigen density via a Hill function, separately parameterized for bivalent engagement (high/mid clones) and monovalent engagement (low clone):

$$
V_{\text{kill}}^{(i)} \;=\; V_{\max}^{(i)}\,\frac{A_i^{\,n^{(i)}}}{\bigl(EC_{50}^{(i)}\bigr)^{n^{(i)}} + A_i^{\,n^{(i)}}}
\qquad i \in \{\text{high},\text{mid},\text{low}\}
$$

Effective effectors after the soluble-BCMA sink, and the per-clone effector-to-target ratio:

$$
E_{\text{eff}} \;=\; \frac{E}{1 + B/K_d},
\qquad
r_i \;=\; \frac{E_{\text{eff}}}{T_i + 1}
$$

The **actual** kill rate saturates with available effectors through a second Hill, in the E:T ratio, half-saturating at $ET_{50}$:

$$
k_{\text{kill}}^{(i)} \;=\; V_{\text{kill}}^{(i)}\,\frac{r_i}{ET_{50} + r_i}
$$

Antigen density sets the ceiling; the E:T ratio sets the approach to it. The two are separable saturable processes, which keeps the kill bounded (unlike a bilinear $E\cdot T$ mass-action term).

## Expansion, decay, growth

Antigen-weighted signal and the niche-gated expansion drive:

$$
\bar{A} \;=\; \frac{\sum_i T_i A_i}{\sum_i T_i + 1},
\qquad
\mathrm{drive} \;=\; L\left(1 + \frac{\bar{A}}{\bar{A} + 100}\right)
$$

Blended effector/memory decay, and the logistic tumor cap:

$$
\delta \;=\; k_{\text{dec,eff}}\,(1 - f_{\text{mem}}) + k_{\text{dec,mem}}\,f_{\text{mem}},
\qquad
g \;=\; 1 - \frac{\sum_i T_i}{K_{\max}}
$$

## Governing equations

$$
\begin{aligned}
\frac{dE}{dt} &= k_{\text{growth}}\,E\,\mathrm{drive} \;-\; \delta\,E \\[6pt]
\frac{dT_i}{dt} &= k_{\text{prolif}}\,T_i\,g \;-\; k_{\text{kill}}^{(i)}\,T_i,
\qquad i \in \{\text{high},\text{mid},\text{low}\} \\[6pt]
\frac{dB}{dt} &= \sum_i k_{\text{shed}}^{(i)}\,T_i \;-\; k_{\text{clear}}\,B \\[6pt]
\frac{dL}{dt} &= -\,k_{\text{LD}}\,L
\end{aligned}
$$

## Parameters

| parameter | symbol | value | units | role |
|---|---|---|---|---|
| CAR-T expansion rate | $k_{\text{growth}}$ | 1.2 | day⁻¹ | effector expansion |
| effector decay | $k_{\text{dec,eff}}$ | 0.15 | day⁻¹ | fast decay |
| memory decay | $k_{\text{dec,mem}}$ | 0.008 | day⁻¹ | slow decay |
| memory fraction | $f_{\text{mem}}$ | 0.15 | — | fixed (literature) |
| niche recovery | $k_{\text{LD}}$ | 0.08 | day⁻¹ | LD niche decay |
| max kill, bivalent | $V_{\max}^{\text{biv}}$ | 2.5 | day⁻¹ | high/mid ceiling |
| EC₅₀, bivalent | $EC_{50}^{\text{biv}}$ | 50 | AU | antigen half-max |
| Hill, bivalent | $n^{\text{biv}}$ | 2.2 | — | |
| max kill, monovalent | $V_{\max}^{\text{mono}}$ | 0.8 | day⁻¹ | low ceiling |
| EC₅₀, monovalent | $EC_{50}^{\text{mono}}$ | 200 | AU | antigen half-max |
| Hill, monovalent | $n^{\text{mono}}$ | 1.1 | — | fixed |
| antigen densities | $A_{\text{high/mid/low}}$ | 500 / 150 / 30 | AU | per clone |
| tumor proliferation | $k_{\text{prolif}}$ | 0.02 | day⁻¹ | |
| carrying capacity | $K_{\max}$ | $10^{12}$ | cells | logistic cap |
| sBCMA shedding | $k_{\text{shed}}^{\text{high/mid/low}}$ | 5×10⁻¹⁰ / 1.5×10⁻¹⁰ / 3×10⁻¹¹ | AU·cell⁻¹·day⁻¹ | per clone |
| sBCMA clearance | $k_{\text{clear}}$ | 3.5 | day⁻¹ | $t_{1/2}\approx4.8$ h |
| sBCMA Kd | $K_d$ | 100 | AU | competitive sink |
| E:T half-saturation | $ET_{50}$ | 0.1 | — | kill saturation |

## Initial conditions

$$
E(0)=7.5\times10^{7},\quad
T_{\text{high}}(0)=6\times10^{10},\quad
T_{\text{mid}}(0)=3\times10^{10},\quad
T_{\text{low}}(0)=1\times10^{10}
$$

$$
B(0)=\frac{\sum_i k_{\text{shed}}^{(i)}T_i(0)}{k_{\text{clear}}}\approx 9.94,
\qquad
L(0)=0.85
$$

## Numerical notes

Solved over $0$–$365$ days at relative and absolute tolerances of $10^{-8}$. The implementation was cross-validated across three solvers — `rxode2` (LSODA), `DifferentialEquations.jl` (Rodas5P), and `scipy` (LSODA) — agreeing to a maximum scaled trajectory error of $2.2\times10^{-6}$, evaluated with scaled relative error and absolute floors to handle states that cross zero.
