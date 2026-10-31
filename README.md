# Heat Exchanger Performance Simulator

**Thermofluids | Heat Transfer | MATLAB | Numerical Methods**

## Overview

A MATLAB-based numerical simulator for evaluating the thermal performance of different heat exchanger configurations. The program calculates heat transfer rates and outlet temperatures from user-specified operating conditions, supporting comparisons between exchanger designs and flow arrangements.

## Engineering Objective

Develop a flexible and physically consistent computational tool to predict heat exchanger performance under varying operating conditions, while solving the coupled temperature equations reliably.

## Supported Configurations

- Single-tube heat exchangers
- Multi-pass shell-and-tube heat exchangers
- Crossflow heat exchangers
- Parallel-flow and counterflow arrangements, where applicable

## Modeling and Numerical Approach

The simulator uses steady-state energy balances and the **log mean temperature difference (LMTD)** method to relate heat transfer to the temperature change of the hot and cold streams.

For a heat exchanger with no external heat losses, the governing energy balance is

$$
\dot Q = \dot m_h c_{p,h}(T_{h,in}-T_{h,out})
       = \dot m_c c_{p,c}(T_{c,out}-T_{c,in}).
$$

The LMTD formulation is

$$
\dot Q = U A F\,\Delta T_{lm},
$$

where $U$ is the overall heat transfer coefficient, $A$ is heat transfer area, and $F$ is a configuration-dependent correction factor when required. The log mean temperature difference is

$$
\Delta T_{lm} = \frac{\Delta T_1-\Delta T_2}{\ln(\Delta T_1/\Delta T_2)}.
$$

The appropriate terminal temperature differences depend on the flow arrangement. Coupled outlet-temperature equations are solved numerically using MATLAB's `fzero`, with physically admissible temperature bounds used to reject nonphysical solutions.

## Implementation

- Organized the program around a menu-driven interface for exchanger selection and operating-condition inputs.
- Implemented energy-balance and LMTD calculations for multiple exchanger types.
- Applied numerical root finding to solve nonlinear, temperature-coupled equations.
- Used physical constraints and convergence checks to improve solution reliability.
- Enabled parametric comparisons across operating conditions and configurations.

## Results

The simulator was used to calculate outlet temperatures, heat transfer rates, and thermal performance metrics for multiple exchanger types and flow arrangements. The numerical approach produced consistent, convergent solutions across the operating conditions investigated, supporting repeatable performance comparisons.

*Note: Specific validation datasets, numerical error metrics, and operating ranges are not included here because they were not provided in the project summary.*

## Tools and Skills

**MATLAB** · Heat Transfer · Thermofluids Modeling · LMTD Method · Energy Balances · Numerical Root Finding (`fzero`)
