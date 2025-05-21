# Asymmetric_Rewiring
This repository contains the code and data necessary to plot the empirical and theoretical results presented in the article "Global Change Asymmetrically Rewires Ecosystems". In this synthesis, we examine how anthropogenic pressures asymmetrically rewire the spatial structure of food webs, with consequences for ecosystem functions and resilience. This repository contains both theoretical modelling code (Julia) and empirical synthesis code (R) used to generate the figures and results in the main text and supplementary materials of this article.

Repository Contents

**empirical_data_plots.R**

This R script generates empirical summary figures from our review of habitat coupling studies (Dataset S1 found in the Data folder of this repository).

The code includes:
- Visualization for Figure 2b and Figure 2c in the main text: Barplots summarizing study results by pressure category (e.g., climate change, nutrient pollution, etc.), mechanism of shift (e.g., changes in accessibility or resource density)
  
Data used in this analysis are available in DatasetS1.csv, which includes directional change (increase, decrease, no change) in habitat coupling for each study, the associated anthropogenic pressure, ecosystem type, and identified mechanism.

Reproducibility:
R version: 4.4
Key R packages: ggplot2, tidyverse


**asymmetric_foodweb_model.jl**

This Julia script contains a dynamical systems model of a generalist food web module, consisting of two basal resources (R1, R2), two intermediate consumers (C1, C2), and a top predator (P). It is used to simulate the ecological consequences of asymmetric rewiring driven by differential impacts from anthropogenic pressures on distinct habitats under both deterministic and stochastic scenarios.

The code includes:
- A five-species ODE system with flexible predator foraging preference (Ω) that can be fixed or density-dependent (ω).
- Simulations over a gradient of K1 values, holding K2 constant, to reflect differential change in habitat productivity.
- Calculation and plotting of: equilibrium densities for all species; predator:consumer biomass ratios; resource, consumer, and predator production; predator habitat coupling; local stability via maximum real eigenvalues of the Jacobian; predator population variability (CV) under Gaussian white noise

With this experiment, we illustrate how food webs are altered across structural (habitat coupling), functional (primary and secondary production), and dynamical (equilibrium and non-equilibrium stability) dimensions as a result of asymmetric rewiring. 

Reproducibility
Julia version: 1.11.2
Key Julia packages: Parameters.jl, LinearAlgebra.jl, ForwardDiff.jl, QuadGK.jl, NLsolve.jl, DifferentialEquations.jl, Plots.jl, Statistics.jl
