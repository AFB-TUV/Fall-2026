using LinearAlgebra, Statistics
using BenchmarkTools, Interpolations, LaTeXStrings, Plots, Roots
using Optim, Random

#c_star_new = K(c_star, grid, β, dudc, f, f_prime, shocks)

function K!(Kg, g, grid, β, dudc, f, f_prime, shocks)
    g_func = LinearInterpolation(grid, g, extrapolation_bc = Line())

    for (i,y) in enumerate(grid)
        function h(c)
            vals = dudc.(g_func.(f(y-c)*shocks)) .* f_prime(y-c) .* shocks
            return dudc(c) - β * mean(vals)
        end
        Kg[i] = find_zero(h, (1e-10, y - 1e-10))
    end
    return Kg
end

function K(g, grid, β, dudc, f, f_prime, shocks)
    K!(similar(g), g, grid, β, dudc, f, f_prime, shocks)
end

function T(w, grid, β, u, f, shocks, Tw = similar(w); compute_policy = false)

    w_func = LinearInterpolation(grid, w, extrapolation_bc = Line())

    if compute_policy
        sigma = similar(w)
    end

    for (i,y) in enumerate(grid)
        objective(c) = u(c) + β * mean(w_func.(f(y-c) .* shocks))
        res = maximize(objective, 1e-10, y)

        if compute_policy
            sigma[i] = Optim.maximizer(res)
        end
        Tw[i] = Optim.maximum(res)
    end

    if compute_policy
        return Tw, sigma
    else
        return Tw
    end
end

isoelastic(c,γ) = isone(γ) ? log(c) : (c^(1-γ) - 1)/(1-γ)
function Model(; α = 0.65, β = 0.95, γ = 1.0, μ = 0.0, s = 0.1, grid = range(1e-6, 4, length = 200), grid_min = 1e-6, grid_max = 4.0, grid_size = 200, u = (c, γ = γ) -> isoelastic(c, γ), dudc = c -> c^(-γ), f = k -> k^α, f_prime = k -> α * k^(α-1))
    return (; α, β, γ, μ, s, grid, grid_min, grid_max, grid_size, u, dudc, f, f_prime)
end

m = Model();

Random.seed!(42)

shock_size = 250
shocks = collect(exp.(m.μ .+ m.s * randn(shock_size)));

function verify_true_policy(m, shocks, c_star)
    (; grid, β, dudc, f, f_prime) = m
    c_star_new = K(c_star, grid, β, dudc, f, f_prime, shocks)

    plot(grid, c_star, label = L"optimal policy $c^*$")
    plot!(grid, c_star_new, label = L"Kc^*")
    plot!(legend = :topleft)
end

verify_true_policy

c_star = (1 - m.α * m.β) * m.grid
verify_true_policy(m, shocks, c_star)

function check_convergence(m, shocks, c_star, g_init; n_iter = 15)
    (; grid, β, dudc, f, f_prime) = m
    g = g_init
    plot(m.grid, g, lw = 2, α = 0.6, label = L"initial condition $c(y) = y$")
    for i in 1:n_iter
        new_g = K(g, grid, β, dudc, f, f_prime, shocks)
        g = new_g
        plot!(grid, g, lw = 2, α = 0.6, label = " ")
    end
    plot!(grid, c_star, color = :black, lw = 2, α = 0.8, label = L"true policy function $c^*$")
    plot!(legend = :topleft)
end

check_convergence(m, shocks, c_star, m.grid, n_iter = 15)