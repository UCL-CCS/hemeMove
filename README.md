# HemeLB: Haemodynamic simulation with lattice Boltzmann

![Main code status](https://github.com/hemelb-codes/hemelb/actions/workflows/main-app.yml/badge.svg)
![Python tools status](https://github.com/hemelb-codes/hemelb/actions/workflows/py-hemetools.yml/badge.svg) 
![Geometry tool status](https://github.com/hemelb-codes/hemelb/actions/workflows/gmy-tool.yml/badge.svg)

HemeLB uses the lattice Boltzmann method to simulate fluid flow in
complex geometries, such as a blood vessel network.

This software was started at University College London and has since
been developed by a large number of people (see AUTHORS). It is open
source under the LGPL license (see LICENSE).

Please see the [doc](doc) folder for more details.

Key features:

- highly scalable
- simulations in complex geometry
- multiple LB velocity sets (D3Q15, D3Q19, D3Q27)
- choice of fluid model (Newtonian via LBGK, MRT, Carreau-Yasuda,
  Casson, Truncated power law)
- various solid wall boundary conditions (SBB, BFL, GZS, JY)
- inlet/outlet boundary conditions by pressure or velocity
