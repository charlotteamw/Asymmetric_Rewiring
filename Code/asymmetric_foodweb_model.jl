
using Parameters: @with_kw, @unpack
using LinearAlgebra: eigvals
using ForwardDiff
using QuadGK: quadgk
using NLsolve
using DifferentialEquations
using Plots
using PyPlot
using Statistics

function adapt_pref(u, p, t)
    return p.ω * u[2] / (p.ω * u[2] + (1 - p.ω) * u[4])
end

function fixed_pref(u, p, t)
    return p.Ω
end

@with_kw mutable struct ModelPar
    a_R1C1 = 0.9
    h_R1C1 = 0.6
    e_R1C1 = 0.7
    a_R2C2 = 0.9
    h_R2C2 = 0.6
    e_R2C2 = 0.7
    a_PC1 = 1.2
    h_PC1 = 0.6
    e_PC1 = 0.7
    a_PC2 = 1.2
    h_PC2 = 0.6
    e_PC2 = 0.7
    m_P = 0.3
    r1 = 2.0
    K1 = 1.0
    r2 = 2.0
    K2 = 2.0
    m_C1 = 0.3
    m_C2 = 0.3
    pref::Function = adapt_pref
    Ω = 0.0
    ω = 0.6
    noise = 0.003
end

## functions for later

function jac(u, model, p)
    ForwardDiff.jacobian(u -> model(u, p, NaN), u)
end 

function habitat_coupling(u, p)
    _, C1, _, C2, P = u
    Ω = p.pref(u, p, 0.0)

    # Intake (functional response) includes handling time
    intake_C1 = (Ω * p.a_PC1 * C1) / (1 + Ω * p.a_PC1 * p.h_PC1 * C1 + (1 - Ω) * p.a_PC2 * p.h_PC2 * C2)
    intake_C2 = ((1 - Ω) * p.a_PC2 * C2) / (1 + Ω * p.a_PC1 * p.h_PC1 * C1 + (1 - Ω) * p.a_PC2 * p.h_PC2 * C2)

    total_intake = intake_C1 + intake_C2

    if total_intake == 0.0
        return 0.0
    end

    p1 = intake_C1 / total_intake
    return 0.5 - abs(0.5 - p1)
end

function basal_prod(u, p)
    R1, _, R2, _, _ = u
    prod_R1 = p.r1 * R1 * (1 - R1 / p.K1)
    prod_R2 = p.r2 * R2 * (1 - R2 / p.K2)
    return prod_R1 + prod_R2
end

function intermediate_prod(u, p)
    R1, C1, R2, C2, _ = u

    intake_C1 = (p.a_R1C1 * R1 * C1) / (1 + p.a_R1C1 * p.h_R1C1 * R1)
    intake_C2 = (p.a_R2C2 * R2 * C2) / (1 + p.a_R2C2 * p.h_R2C2 * R2)

    return p.e_R1C1 * intake_C1 + p.e_R2C2 * intake_C2
end

function top_pred_prod(u, p)
    _, C1, _, C2, P = u
    Ω = p.pref(u, p, 0.0)

    intake_C1 = (Ω * p.a_PC1 * C1 * P) / (1 + Ω * p.a_PC1 * p.h_PC1 * C1 + (1 - Ω) * p.a_PC2 * p.h_PC2 * C2)
    intake_C2 = ((1 - Ω) * p.a_PC2 * C2 * P) / (1 + Ω * p.a_PC1 * p.h_PC1 * C1 + (1 - Ω) * p.a_PC2 * p.h_PC2 * C2)

    return p.e_PC1 * intake_C1 + p.e_PC2 * intake_C2
end



### Food Web Model

function model!(du, u, p, t)
    @unpack r1, K1, r2, K2 = p
    @unpack a_R1C1, h_R1C1, e_R1C1, a_R2C2, h_R2C2, e_R2C2, m_C1, m_C2 = p
    @unpack a_PC1, h_PC1, e_PC1, a_PC2, h_PC2, e_PC2, m_P, Ω, ω  = p
    R1, C1, R2, C2, P = u

    Ω = p.pref(u, p, t)

    int_R1C1 = a_R1C1 * R1 * C1 / (1 + a_R1C1 * h_R1C1 * R1)
    int_R2C2 = a_R2C2 * R2 * C2 / (1 + a_R2C2 * h_R2C2 * R2)
    denom_PC1C2 = 1 + Ω * a_PC1 * h_PC1 * C1 + (1 - Ω) * a_PC2 * h_PC2 * C2
    num_PC1 = Ω * a_PC1 * C1 * P
    num_PC2 = (1 - Ω) * a_PC2 * C2 * P

    du[1] = r1 * R1 * (1 - R1 / K1) - int_R1C1
    du[2] = e_R1C1 * int_R1C1 - (num_PC1/denom_PC1C2) - m_C1 * C1
    du[3] = r2 * R2 * (1 - R2 / K2) - int_R2C2
    du[4] = e_R2C2 * int_R2C2 - (num_PC2/denom_PC1C2) - m_C2 * C2
    du[5] = (e_PC1 * num_PC1 + e_PC2 * num_PC2) / denom_PC1C2 - m_P * P

    return du
end

function model(u, ModelPar, t)
    du = similar(u)
    model!(du, u, ModelPar, t)
    return du
end

## time series for food web model with ModelPar
let
    u0 = [0.8, 0.4, 0.8, 0.4, 0.3]
    t_span = (0, 2000.0)
    p = ModelPar()
    prob = ODEProblem(model!, u0, t_span, p)
    sol = solve(prob)
    model_ts = figure()
    for i in 1:5
        PyPlot.plot(sol.t, getindex.(sol.u, i), label = ["R1", "C1", "R2", "C2", "P"][i])
    end
    xlabel("Time")
    ylabel("Density")
    legend()
    return model_ts
end

### Details for all the calculations that follow
K1_vals = 0.76:0.005:1.395
u0 = [0.8, 0.4, 0.8, 0.4, 0.3]
p = ModelPar()
ts = range(1000, 2000, length = 1001)  # Time steps
t_span = (0.0, 2000.0)  # Time span

## eq data frame
eq_hold = fill(0.0,length(K1_vals),6)
### Equilibrium densities
for i=1:length(K1_vals)
    p = ModelPar(K1 = K1_vals[i])
    u0 = [0.8, 0.4, 0.8, 0.4, 0.3]
    prob = ODEProblem(model!, u0, t_span, p)
    sol = solve(prob, reltol = 1e-8, abstol = 1e-8)
    grid = sol(ts)
    eq = mean(grid, dims=2)
    eq_hold[i,1] = K1_vals[i]
    eq_hold[i,2:end] = eq'
    println(eq_hold[i,:])
end

## plot equilibrium densities 
using Plots

eq_R1 = Plots.plot(eq_hold[:,1], eq_hold[:,2], legend = false, lw=2.0, colour = "black", xlabel = "", ylabel = "", xflip = true, xtickfont = font(12), ytickfont = font(12))

eq_C1 = Plots.plot(eq_hold[:,1], eq_hold[:,3], legend = false, lw=2.0, colour = "black", xlabel = "", ylabel = "", xflip = true, xtickfont = font(12), ytickfont = font(12))

eq_R2 = Plots.plot(eq_hold[:,1], eq_hold[:,4], legend = false, lw=2.0, linecolour = "darkorange", xlabel = "", ylabel = "", xflip = true, xtickfont = font(12), ytickfont = font(12))

eq_C2 = Plots.plot(eq_hold[:,1], eq_hold[:,5], legend = false, lw=2.0, linecolour = "green", xlabel = "", ylabel = "", xflip = true, xtickfont = font(12), ytickfont = font(12))

eq_P = Plots.plot(eq_hold[:,1], eq_hold[:,6], legend = false, grid = false, lw = 5.0, color = "black", xlabel = "", ylabel = "", xflip = true, xtickfont = font(14), ytickfont = font(14))

## plot P:C biomass ratio
PC_biomass = Plots.plot(eq_hold[:,1], eq_hold[:,6] ./ eq_hold[:,3], legend = false, lw = 5.0, color = "black", xlabel = "", ylabel = "", 
xflip = true, grid = false, xtickfont = font(14), ytickfont = font(14))

## calculate max real eigs 
maxeig_hold = fill(0.0,length(K1_vals),2)

for i=1:length(K1_vals)
    p = ModelPar(K1 = K1_vals[i])
    u0 = [0.8, 0.4, 0.8, 0.4, 0.3]
    prob = ODEProblem(model!, u0, t_span, p)
    sol = solve(prob, reltol = 1e-8, abstol = 1e-8)
    grid = sol(ts)
    eq = nlsolve((du, u) -> model!(du, u, p, 0.0), grid.u[end]).zero
    coup_jac = jac(eq, model, p)
    max_eig = maximum(real.(eigvals(coup_jac)))
    maxeig_hold[i,1] = K1_vals[i]
    maxeig_hold[i,2] = max_eig
    println(maxeig_hold[i,:])
end

## plot max real eig 
max_eig = Plots.plot(maxeig_hold[:,1], maxeig_hold[:,2], legend = false, lw = 5.0, colour = "black", xlabel = "", ylabel = "", xflip = true, 
yflip = true, grid = false, xtickfont = font(14), ytickfont = font(14))



### Degree coupling

coupling_hold= zeros(length(K1_vals), 2)


for i=1:length(K1_vals)
    p = ModelPar(K1 = K1_vals[i])
    u0 = [0.8, 0.4, 0.8, 0.4, 0.3]
    prob = ODEProblem(model!, u0, t_span, p)
    sol = solve(prob, reltol = 1e-8, abstol = 1e-8)
    grid = sol(ts)
    eq = mean(grid, dims=2)
    coupling_val = habitat_coupling(eq, p)
    coupling_hold[i,1] = K1_vals[i]
    coupling_hold[i,2] = coupling_val
    println(coupling_hold[i,:])
end

degree_coupling_plot = Plots.plot(coupling_hold[:,1], coupling_hold[:,2], legend = false, grid = false, lw = 5.0, color = "black", xlabel = "", ylabel = "", 
xflip = true, xtickfont = font(14), ytickfont = font(14))


### Secondary (predator) Production

pred_prod_hold= zeros(length(K1_vals), 2)

for i=1:length(K1_vals)
    p = ModelPar(K1 = K1_vals[i])
    u0 = [0.8, 0.4, 0.8, 0.4, 0.3]
    prob = ODEProblem(model!, u0, t_span, p)
    sol = solve(prob, reltol = 1e-8, abstol = 1e-8)
    grid = sol(ts)
    eq = mean(grid, dims=2)
    pred_prod_val = top_pred_prod(eq, p)
    pred_prod_hold[i,1] = K1_vals[i]
    pred_prod_hold[i,2] = pred_prod_val
    println(pred_prod_hold[i,:])
end

pred_prod_plot = Plots.plot(pred_prod_hold[:,1], pred_prod_hold[:,2], legend = false, grid = false, lw = 5.0, color = "black", xlabel = "", ylabel = "", 
xflip = true, xtickfont = font(14), ytickfont = font(14))

### Basal (primary) production 

basal_prod_hold= zeros(length(K1_vals), 2)

for i=1:length(K1_vals)
    p = ModelPar(K1 = K1_vals[i])
    u0 = [0.8, 0.4, 0.8, 0.4, 0.3]
    prob = ODEProblem(model!, u0, t_span, p)
    sol = solve(prob, reltol = 1e-8, abstol = 1e-8)
    grid = sol(ts)
    eq = mean(grid, dims=2)
    basal_prod_val = basal_prod(eq, p)
    basal_prod_hold[i,1] = K1_vals[i]
    basal_prod_hold[i,2] = basal_prod_val
    println(basal_prod_hold[i,:])
end

basal_prod_plot = Plots.plot(basal_prod_hold[:,1], basal_prod_hold[:,2], legend = false, grid = false, lw = 5.0, color = "black", xlabel = "", ylabel = "", xflip = true, xtickfont = font(14), ytickfont = font(14))

### Intermediate Production

int_prod_hold = zeros(length(K1_vals), 2)

for i in 1:length(K1_vals)
    p = ModelPar(K1 = K1_vals[i])
    u0 = [0.8, 0.4, 0.8, 0.4, 0.3]
    prob = ODEProblem(model!, u0, t_span, p)
    sol = solve(prob, reltol = 1e-8, abstol = 1e-8)
    grid = sol(ts)
    eq = mean(grid, dims=2)
    int_prod_val = intermediate_prod(eq, p)
    int_prod_hold[i, 1] = K1_vals[i]
    int_prod_hold[i, 2] = int_prod_val
    println(int_prod_hold[i, :])
end

int_prod_plot = Plots.plot(int_prod_hold[:,1], int_prod_hold[:,2], legend = false, grid = false, lw = 5.0, color = "black", xlabel = "", ylabel = "", xflip = true, xtickfont = font(14), ytickfont = font(14))
###### Stochastic model for CV

## Adding stochasticity to model using gaussian white noise (SDEproblem)
function stochmodel!(du, u, p2, t)
    @unpack  noise = p2

    du[1] = noise * u[1]
    du[2] = noise * u[2]
    du[3] = noise * u[3]
    du[4] = noise * u[4]
    du[5] = noise * u[5]

    return du 
end

## time series for stochastic food web model with ModelPar
let
    u0 = [0.8, 0.4, 0.8, 0.4, 0.3]
    t_span = (0, 2000.0)
    p = ModelPar(K1 = 0.76)
    prob = SDEProblem(model!, stochmodel!, u0, t_span, p)
    sol = solve(prob)
    model_ts = figure()
    for i in 1:5
        PyPlot.plot(sol.t, getindex.(sol.u, i), label = ["R1", "C1", "R2", "C2", "P"][i])
    end
    xlabel("Time")
    ylabel("Density")
    legend()
    return model_ts
end

##cv data frame
cv_hold = zeros(length(K1_vals), 2)
stdhold = fill(0.0, length(K1_vals), 2)
meanhold = fill(0.0, length(K1_vals), 2)

#### CV of predator population
for i=1:length(K1_vals)
    p = ModelPar(K1 = K1_vals[i])
    u0 = [0.8, 0.4, 0.8, 0.4, 0.3]
    prob_stoch = SDEProblem(model!, stochmodel!, u0, t_span, p)
    sol_stoch = solve(prob_stoch, reltol = 1e-15)
    grid_sol = sol_stoch(ts)

    # Recording the statistics
    cv_hold[i, 1] = K1_vals[i]
    stdhold[i, 1] = std(grid_sol[5, :])
    meanhold[i, 1] = mean(grid_sol[5, :])
    cv_hold[i, 2] = stdhold[i, 1] / meanhold[i, 1]

    println(cv_hold[i,:])
end

cv_plot = Plots.plot(cv_hold[:,1], cv_hold[:,2], legend = false, grid = false, lw = 5.0, color = "black", xlabel = "", ylabel = "", xflip = true, xtickfont = font(14), ytickfont = font(14))


