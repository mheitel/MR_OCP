# Manifold-Based Model Reduction of Optimal Control Problems

This project provides a MATLAB script for solving optimal control problems (OCPs). Manifold-Based Model Reduction Techniques reduce both stiffness and order of the underlying differential equations.

## Getting Started

Download the complete repository and run the enzyme example on your computer with MATLAB. Please see the following instructions for a correct use.

### Prerequisites

You need a proper installation of MATLABi, [CasADi](https://web.casadi.org/), and optionally [IPOPT](https://github.com/coin-or/Ipopt).  

### Installing

Just download this repository. Often, it is recommended to install the HSL routine MA27 for a faster solution of the underlying nonlinear program with IPOPT.

## Example

If you want to test your installation, please open MATLAB on you system and switch to the examples folder. Then just run the enzyme example with the following command.

```
enzyme_main
```
## Software 

* [CasADi](https://web.casadi.org/) - Tool for nonlinear optimization and algorithmic differentiation
* [IPOPT](https://github.com/coin-or/Ipopt) - Nonlinear Program Solver

## Authors

* **Marcus Heitel** - Institute for Numerical Mathematics, Ulm University, Ulm

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
