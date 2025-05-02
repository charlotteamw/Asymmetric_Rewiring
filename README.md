# Asymmetric_Rewiring
This repository accompanies the article "Global Change Asymmetrically Rewires Ecosystems". In this synthesis, we examine how anthropogenic pressures asymmetrically rewire the spatial structure of food webs, with consequences for ecosystem functions and resilience. This repository contains both theoretical modelling code (Julia) and empirical synthesis code (R) used to generate the figures and results in the main text and supplementary materials of this article.

Repository Contents

**empirical_data_plots.R**

This R script generates empirical summary figures from our review of habitat coupling studies (Dataset S1 found in the Data folder of this repository).

The code includes:
- Barplots summarizing study results by: pressure category (e.g., climate change, nutrient pollution, etc.), mechanism of shift (e.g., changes in accessibility or resource density)
- Visualization logic for Figure 2b and Figure 2c in the main text.

Data used in this analysis are available in DatasetS1.csv, which includes directional change (increase, decrease, no change) in habitat coupling for each study, the associated anthropogenic pressure, ecosystem type, and identified mechanism.

Reproducibility:
R version: ≥4.1
Key R packages: ggplot2, tidyverse, patchwork


**asymmetric_foodweb_model.jl**

This Julia script contains a dynamical model of a generalist food web module, consisting of two basal resources (R1, R2), two intermediate consumers (C1, C2), and a top predator (P). It is used to simulate the ecological consequences of asymmetric changes in habitat productivity (i.e., changes to K1, the carrying capacity of R1) under both deterministic and stochastic scenarios.

The code includes:
- A five-species ODE system with flexible predator foraging preference (Ω) that can be fixed or density-dependent (ω).
- Simulations over a gradient of K1 values, holding K2 constant, to reflect increasing asymmetry in habitat quality.
- Calculation and plotting of: equilibrium densities for all species; predator:consumer biomass ratios; primary, intermediate, and secondary (predator) production; predator habitat coupling; local stability via maximum real eigenvalues of the Jacobian; predator population variability (CV) under Gaussian white noise

The model is used to illustrate how food web responses shift across structural (coupling), functional (production), and dynamical (stability) dimensions in response to differential impacts of anthropogenic pressure of distinct habitats.

Reproducibility
Julia version: ≥1.9
Key Julia packages: DifferentialEquations.jl, Plots.jl, NLsolve.jl, ForwardDiff.jl
