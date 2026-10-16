using LinearAlgebra, Statistics
using BenchmarkTools, Interpolations, LaTeXStrings, Random, Roots
using Plots

function coleman_egm(g, k_grid, β, u_prime, u_prime_inv, f, f_prime, shocks)

    #allocate memory for consumption on grid points
    c = similar(kgrid)

    #update consumption value
    for (i,k) in enumerate(k_grid)
        vals = u_prime.(g.(f(k)*shocks)) .* f_prime(k) .* shocks
        c[i] = u_prime_inv(β * mean(vals))
    end

    #determine endogenous grid
    y = k_grid + c #y_i = k_i + c_i 
    # update policy function, return
    Kg = LinearInterpolation(y, c, extrapolation_bc = Line())
    Kg_f(x) = Kg(x)
    return Kg_f
end

function Model (; α = 0.65, β = 0.95, γ = 1.0, μ = 0.0, s = 0.1, grid_min = 1e-6, grid_max = 4.0, grid_size = 200, 
                u = γ == 1? log : c -> (c^(1-γ) - 1)/(1-γ), u_prime = c -> c^(-γ), 
                f = k -> k^α, f_prime = k -> α * k^(α-1), 
                grid = range(grid_min, grid_max, length = grid_size))
                return (; α, β, γ, μ, s, grid_min, grid_max, grid_size, u, u_prime, f, f_prime, grid)
end
