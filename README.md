# dFBA: Dynamic Flux Balance Analysis of Pyruvate Routing in Yeast
dFBA: Dynamic Flux Balance Analysis of pyruvate routing in Saccharomyces cerevisiae under aerobic and anaerobic conditions with CDC19 (pyruvate kinase) constraints. MATLAB pipeline for simulating short-term metabolic responses to glucose pulses, with R scripts for visualization

# Overview

This repository contains the complete code and data for dynamic Flux Balance Analysis (dFBA) simulations of *Saccharomyces cerevisiae* central carbon metabolism under glucose pulse conditions. The study investigates oxygen-dependent pyruvate routing under CDC19 (pyruvate kinase) constraints, comparing single-pulse (SP) and repetitive-pulse (RP) conditions under both aerobic and anaerobic environments.

**Key features:**
- MATLAB pipeline for time-resolved dFBA simulations
- Four experimental conditions: Aerobic_SP, Aerobic_RP, Anaerobic_SP, Anaerobic_RP
- Pyruvate kinase (CDC19) constraints (50% anaerobic, 80% aerobic during feed phase)
- 12 metabolic parameters including glucose uptake, growth rate, CO₂, O₂, and fermentation products

## 🔧 Requirements

### MATLAB (R2023b or later)
- [COBRA Toolbox v3.0](https://opencobra.github.io/)
- [GLPK](https://www.gnu.org/software/glpk/) (or any supported LP solver)
- Statistics and Machine Learning Toolbox
- Parallel Computing Toolbox (optional, for speed)
