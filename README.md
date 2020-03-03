# Open-source SHEMAT-Suite #

**SHEMAT-Suite (Simulator for HEat and MAss Transport)** is a
numerical code for computing flow, heat and species transport
equations in porous media. The governing equations of the code are the
groundwater flow equation, the heat transport equation and the species
transport equation.

SHEMAT-Suite includes parameter estimation and data assimilation
approaches, both stochastic (Monte Carlo, ensemble Kalman filter) and
deterministic (Bayesian inversion using automatic differentiation for
calculating derivatives).

SHEMAT-Suite is written in Fortran-95.

SHEMAT-Suite uses finite-difference discretization on a Cartesian grid
with x, y, z coordinates and block centered nodes.  The system of
equations can be solved explicitly, implicitly or semi-implicitly.

[More detailed information on the
equations](https://git.rwth-aachen.de/SHEMAT-Suite/SHEMAT-Suite-open/-/wikis/equations)


## Installation ##

Guide to installing SHEMAT-Suite:
[Installation](https://git.rwth-aachen.de/SHEMAT-Suite/SHEMAT-Suite-open/-/wikis/howtos/installation_requirements)

## Tutorial ##

If you want to run SHEMAT-Suite for the first time you should check
the [SHEMAT-Suite Tutorial](https://git.rwth-aachen.de/SHEMAT-Suite/SHEMAT-Suite-open/-/wikis/tutorial).

## For Developers ##

[Fortran-95 styleguide](https://git.rwth-aachen.de/SHEMAT-Suite/SHEMAT-Suite-open/-/wikis/howtos/styleguide) for SHEMAT-Suite.

[Doxygen documentation](https://git.rwth-aachen.de/SHEMAT-Suite/SHEMAT-Suite-open/-/wikis/howtos/doxygen).

## Most visited pages ##

[Input File](https://git.rwth-aachen.de/SHEMAT-Suite/SHEMAT-Suite-open/-/wikis/tutorial/input_file)

[Compilation](https://git.rwth-aachen.de/SHEMAT-Suite/SHEMAT-Suite-open/-/wikis/tutorial/compilation)

[Equations and Discretization](https://git.rwth-aachen.de/SHEMAT-Suite/SHEMAT-Suite-open/-/wikis/equations)

## Howtos ##
[howtos](https://git.rwth-aachen.de/SHEMAT-Suite/SHEMAT-Suite-open/-/wikis/howtos)

## Property Modules ##
[props](https://git.rwth-aachen.de/SHEMAT-Suite/SHEMAT-Suite-open/-/wikis/props)

## User Directories ##
[user](https://git.rwth-aachen.de/SHEMAT-Suite/SHEMAT-Suite-open/-/wikis/user)

## Branches ##
[branches](https://git.rwth-aachen.de/SHEMAT-Suite/SHEMAT-Suite-open/-/wikis/branches)

<!-- ![shemat_suite_branches](https://git.rwth-aachen.de/SHEMAT-Suite/SHEMAT-Suite-open/-/wikis/uploads/1d66b69ef413debec99f1d86e7aa09a2/shemat_suite_branches.png) -->

<img src="https://git.rwth-aachen.de/SHEMAT-Suite/SHEMAT-Suite-open/-/wikis/uploads/1d66b69ef413debec99f1d86e7aa09a2/shemat_suite_branches.png"  width="300">
<!-- <img src="https://git.rwth-aachen.de/SHEMAT-Suite/SHEMAT-Suite-open/-/wikis/uploads/1d66b69ef413debec99f1d86e7aa09a2/shemat_suite_branches.png"  width="120" height="120"> -->

## Outline of the directories in branch master ##

![shemat_suite_forward](https://git.rwth-aachen.de/SHEMAT-Suite/SHEMAT-Suite-open/-/wikis/uploads/6c0aef0549d4ffbd608dd7a077507c17/shemat_suite_forward.png)

  **Table:** Important code developments of
  SHEMAT-Suite.
  * ($`^{\mathrm{x}}`$) Functionalities not available in
	the open-source package.
  * (\*) Simplified functionality available in
	the open-source package.
  * (\*\*) SHEMAT-Suite functionality available
	open-source, additional software required.\

| **Newly implemented functionality**                                                                              | **Key Reference**                               |
| ---------------------------------------------------------------------------------------------------------------- | ----------------------------------------------- |
| Inverse parameter estimation based on automatic differentiation                                                  | Rath et al., 2006                               |
| Latent heat effects due to freezing and melting                                                                  | Mottaghy and Rath, 2006                         |
| Monte Carlo techniques for uncertainty quantification and reduction                                              | Vogt et al., 2010                               |
| Borehole heat exchanger module(\*)                                                                               | Mottaghy and Dijkshoorn, 2012                   |
| Shared-memory parallelization                                                                                    | Wolf, 2011                                      |
| Data assimilation based on the ensemble Kalman Filter                                                            | Vogt et al., 2012                               |
| Multi-phase flow module using automatic differentiation ($`^{\mathrm{x}}`$)                                      | Büsing et al., 2014                             |
| Distributed-memory parallelization ($`^{\mathrm{x}}`$)                                                           | Rostami and Bücker, 2014                        |
| Anisotropic flow module using the full permeability tensor ($`^{\mathrm{x}}`$)                                   | Chen et al., 2016                               |
| Supercritical water/steam module using automatic differentiation ($`^{\mathrm{x}}`$)                             | Büsing et al., 2017                             |
| Optimal borehole positioning with respect to reservoir characterization via optimal experimental design (\*\*)   | Seidler et al., 2016                            |
| Efficient two-phase flow in heterogeneous porous media using exact Jacobians ($`^{\mathrm{x}}`$)                 | Büsing, 2020                                    |

1. Rath, V., Wolf, A., & Bücker, H. M., Joint three-dimensional
   inversion of coupled groundwater flow and heat transfer based on
   automatic differentiation: sensitivity calculation, verification,
   and synthetic examples, Geophysical Journal International, 167(1),
   453–466 (2006).  http://dx.doi.org/10.1111/j.1365-246x.2006.03074.x
2. Mottaghy, D., & Rath, V., Latent heat effects in subsurface heat
   transport modelling and their impact on palaeotemperature
   reconstructions, Geophysical Journal International, 164(1), 236–245
   (2006).  http://dx.doi.org/10.1111/j.1365-246x.2005.02843.x
3. Vogt, C., Mottaghy, D., Wolf, A., Rath, V., Pechnig, R., & Clauser,
   C., Reducing temperature uncertainties by stochastic geothermal
   reservoir modelling, Geophysical Journal International, 181(1),
   321–333 (2010).  http://dx.doi.org/10.1111/j.1365-246x.2009.04498.x
4. Mottaghy, D., & Dijkshoorn, L., Implementing an effective finite
   difference formulation for borehole heat exchangers into a heat and
   mass transport code, Renewable Energy, 45(), 59–71 (2012).
   http://dx.doi.org/10.1016/j.renene.2012.02.013
5. Wolf, A., Ein softwarekonzept zur hierarchischen parallelisierung
   von stochastischen und deterministischen inversionsproblemen auf
   modernen ccnuma-plattformen unter nutzung automatischer
   programmtransformation. (Doctoral dissertation) (2011). RWTH Aachen
   University.
6. Vogt, C., Marquart, G., Kosack, C., Wolf, A., & Clauser, C.,
   Estimating the permeability distribution and its uncertainty at the
   egs demonstration reservoir soultz-sous-forêts using the ensemble
   Kalman filter, Water Resources Research, 48(8), (2012).
   http://dx.doi.org/10.1029/2011wr011673
7. Büsing, H., Willkomm, J., Bischof, C. H., & Clauser, C., Using
   exact jacobians in an implicit newton method for solving multiphase
   flow in porous media, International Journal of Computational
   Science and Engineering, 9(5/6), 499 (2014).
   http://dx.doi.org/10.1504/ijcse.2014.064535
8. Rostami, M. A., & H. M. B\"ucker, Preservation of non-uniform
   memory architecture characteristics when going from a nested OpenMP
   to a hybrid MPI/OpenMP approach, In M. S. Obaidat, J. Kacprzyk, &
   T. {\"O}ren, SIMULTECH 2014, Proceedings of the 4th International
   Conference on Simulation and Modeling Methodologies, Technologies
   and Applications, Vienna, Austria, August~28--30, 2014
   (pp. 286–291) (2014). : SciTePress.
9. Chen, T., Clauser, C., Marquart, G., Willbrand, K., & Büsing, H.,
   Modeling anisotropic flow and heat transport by using mimetic
   finite differences, Advances in Water Resources, 94(), 441–456
   (2016).  http://dx.doi.org/10.1016/j.advwatres.2016.06.006
10. Büsing, H., Vogt, C., Ebigbo, A., & Klitzsch, N., Numerical study
    on co2leakage detection using electrical streaming potential data,
    Water Resources Research, 53(1), 455–469 (2017).
    http://dx.doi.org/10.1002/2016wr019803
11. Seidler, R., Padalkina, K., Bücker, H. M., Ebigbo, A., Herty, M.,
    Marquart, G., & Niederau, J., Optimal experimental design for
    reservoir property estimates in geothermal exploration,
    Computational Geosciences, 20(2), 375–383 (2016).
    http://dx.doi.org/10.1007/s10596-016-9565-4
12. H. Büsing, Efficient solution techniques for two-phase flow in
    heterogeneous porous media using exact Jacobians, Computational
    Geosciences, In Review (2020).
