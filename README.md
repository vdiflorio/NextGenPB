# NextGenPB 
-----------  
Copyright (C) 2021-2025 Vincenzo Di Florio

Copyright (C) 2019-2025 Carlo de Falco

Copyright (C) 2020-2021 Martina Politi

This software is distributed under the terms
the terms of the GNU/GPL licence v3

# Overview
----------

**NextGenPB** is a high-performance solver for the linearized Poisson–Boltzmann equation (PBE), built on an adaptive octree mesh.
It efficiently computes electrostatic potentials in heterogeneous dielectric media using a flexible, hierarchical discretization scheme.

The equation solved is:


$$
-\mathrm{div} \left( \varepsilon_0 \varepsilon_r \nabla \varphi \right) + \kappa^2 \varphi = \rho^f
$$

on a rectangular domain.

# Nonlinear Poisson--Boltzmann Extension
---

This repository also contains a course project extension that adds a
Newton-based solver for the nonlinear Poisson--Boltzmann equation.

The implementation is available in the `nonlinear-solver` branch.

The nonlinear formulation replaces the linear ionic contribution with a
hyperbolic sine term. At every Newton iteration, the Jacobian contains
the corresponding hyperbolic cosine contribution.

The solver mode is selected in the parameter file:

```text
linearized = 1   # original linearized solver
linearized = 0   # nonlinear Newton solver
```

The main modifications are located in:

```text
include/pb_class.h
src/pb_class.cpp
src/poisson_boltzmann.cpp
```

The nonlinear extension adds the following functions:

```text
assemble_newton_system
newton_solve
```


See `REPRODUCE.md` for step-by-step instructions to reproduce our
course project results (real-molecule test, clamping demonstration,
and weak/strong scalability tests).


# Documentation & Tutorials
---

Comprehensive installation instructions, examples, and usage guides are available here:

[NextGenPB Tutorial and Guide](https://vdiflorio.github.io/nextgenpb_tutorial/)

---

# Citation

If you use **NextGenPB** in your research, please cite the following article:

>Di Florio, V., Ansalone, P., Siryk, S. V., Decherchi, S., De Falco, C., & Rocchia, W. (2025). NextGenPB: An analytically-enabled super resolution tool for solving the Poisson-Boltzmann Equation featuring local (de) refinement. Computer Physics Communications, 109816.
