# SystemAnalysis
SystemAnalysis is a Mathematica package for preprocessing in a fully automatic way first order systems of differential equations (developed as part of a PhD thesis at the Research Institute for Symbolic Computations, Johannes Kepler University, in Linz). More precisely, it's a package who's aim is to **analyse first order coupled differential systems depending on a parameter ε** and to determine, before any solving is attempted, the *minimal* ε-orders to which the inhomogeneous part has to be expanded.

## Aim

Systems of the form

$$\alpha_i(x,d)\frac{d}{dx}I_i(x) = \sum_j \beta_{i,j}(x,d)\,I_j(x) + \sum_k \gamma_{i,k}(x,d)B_k(x)$$

arise typically in particle physics, where the $I_i$ are *master integrals* (MIs) and the $B_k$ are *base case integrals* (BIs) whose ε-expansions must be supplied by the user. Since these systems can only be solved up to a finite order in ε, computing the BIs to an unnecessarily high order is the main bottleneck.

`SystemAnalysis` answers the question: *given a target order for the master integrals, what is the cheapest set of inputs and equations that gets me there?* It does the following:

1. **Triangularizes** the system, splitting it into irreducible subsystems,
2. **Uncouples** each subsystem in several schemes (Gauß, Zürcher), comparing the resulting ε-orders as well as the orders of the associated ODEs/recurrences,
3. **Recombines** the local corrections into a global answer for the full system.

The output is the list of MIs with the ε-order and the order of the higher-order ODE (HODE) needed for each, the list of BIs with the ε-order to which each must be provided, and the corresponding set of equations.

## Requirements

`SystemAnalysis` is built on top of several RISC packages, which must be loaded **before** it:

```mathematica
<< Sigma.m;               (* C. Schneider *)
<< HarmonicSums.m;        (* J. Ablinger *)
<< SumProduction.m;       (* C. Schneider *)
<< EvaluateMultiSums.m;   (* C. Schneider *)
<< OreSysG.m;             (* S. Gerhold *)
<< SolveCoupledSystem.m;  (* C. Schneider *)

<< SystemAnalysis.m;
```

The packages can be obtained from:

- [**Sigma**](https://www.risc.jku.at/research/combinat/software/Sigma/) — a summation package by Carsten Schneider
- [**HarmonicSums**](https://www.risc.jku.at/research/combinat/software/HarmonicSums/) — a package for allows to deal with nested sums such as harmonic sums, S-sums, cyclotomic sums and cyclotmic S-sums as well as iterated integrals such as harmonic polylogarithms, multiple polylogarithms and cyclotomic polylogarithms in an algorithmic fashion by Jakob Ablinger
- [**SumProduction and EvaluateMultiSums**](https://www3.risc.jku.at/research/combinat/software/EvaluateMultiSums/index.php) — packages for the simplification of definite multi-sums and large multi-sums in terms of indefinite nested sums and products by Carsten Schneider
- [**SolveCoupledSystems**](https://www3.risc.jku.at/research/QFT/software.html) — a package for solving coupled systems of differential and difference equations in terms of nested sums and products by Carsten Schneider
- [**OreSysG**](https://www3.risc.jku.at/research/combinat/software/ergosum/RISC/OreSys.html) — a package for uncoupling systems of linear Ore operator equations by Stefan Gerhold

## Basic commands

### `getMI[sys, ord, x]`

Helper returning the list of `{master integral, order}` pairs, i.e. all the MIs of `sys` requested up to the same ε-order `ord`.

### `analysisSystem[sys, x, H, J, W, R, ord, d, Ep, opts]`

The single, centralised entry point performing the full analysis.

| Argument | Meaning |
|---|---|
| `sys`  | the system to analyse, as a list of equations |
| `x`    | the differentiation variable |
| `H`    | head of the MIs/BIs (e.g. `J`, or `HFORM3l1a` for physics input) |
| `J`    | head used to rename the integrals internally |
| `W`    | head used for the dummy variables during the uncoupling |
| `R`    | head used for the abstract right-hand side vector |
| `ord`  | list of `{MI, target ε-order}` pairs, typically built with `getMI` |
| `d`    | symbolic spacetime dimension |
| `Ep`   | the symbol ε, used through the relation `d = 4 - Ep` |

Options:

| Option | Values | Effect |
|---|---|---|
| `"Speed"` | `0`, `1`, `2` | `0` performs all possible uncouplings; `1` uses the time of the first uncouplings as a time limit for the remaining ones; `2` launches Gauß and Zürcher in parallel and keeps the first to finish |
| `"Scheme"` | `"Default"`, `"EpsOrder"`, `"DiffEqn"`, `"RecEqn"` | selection criterion: combined cost $-n_\varepsilon + n_d + n_r$, or minimisation of the ε-order, of the ODE order, or of the recurrence order alone |
| `"Save"` | `True`, `False` | writes intermediate save states (`Backup_sys/`) and final solutions (`Output_sys/`); strongly recommended for long computations |
| `"Parallelization"` | `True`, `False` | tries to parallelise the uncouplings |

Output:

```mathematica
{"MasterIntegralsUpdated" -> {..., {I[i], epsOrder, HODEorder}, ...},
 "BaseIntegralsUpdated"   -> {..., {J[j], epsOrder}, ...},
 "Equations"              -> {..., {I[i], lhs == rhs}, ...}}
```

where `HODEorder = 0` means the master integral is obtained from a linear relation rather than from a differential equation.

Two lower-level utilities are also exported: `systemClustering[sys, x, varRoot, choice]`, which splits a system into its irreducible subsystems, and `ToMatrixForm[sys, x, varRoot]`, which returns the matrix representation of a system.
