using LinearAlgebra, Statistics, Interpolations, NLsolve
using BenchmarkTools, LaTeXStrings, Optim, Random
using Distributions
using Optim: converged, maximum, maximizer, minimizer, iterations
using Plots

u(x) = log(x)
du(x) = 1/x

function consumer_problem(; r = 0.01, β = 0.96, Π = [0.6 0.4; 0.05 0.95], zgrid = [0.5, 1.0], b = 0.0, grid_max = 16, grid_size = 50)
    R = 1+r 
    agrid = range(-b, grid_max, length = grid_size)
    return (; r, R, β, b, Π, z_vals, agrid) 
end

function T!(cp, V, out; ret_policy = false)
    (; R, Π, β, b, agrid, zgrid) = cp
    z_idx = 1:length(zgrid)

    #value function when shock index is z_i
    vf = LinearInterpolation((agrid, z_idx), V; extrapolation_bc = Interpolations.Flat())
    opt_lb = 1e-8

    #RHS of Bellman
    for (i_z, z) in enumerate(zgrid)
        for (i_a, a) in enumerate(agrid)
            function obj(c)
                EV = dot(vf.(R * a + z - c, z_idx), Π[i_z, :])
                return u(c)+ β*EV 
            end
            res = maximize(obj, opt_lb, R .* a .+ z .+ b)
            converged(res) || error("Didn't converge")

            if ret_policy
                out[i_a, i_z] = maximizer(res)
            else
                out[i_a, i_z] = maximum(res)
            end
        end
    end
    return out
end

T(cp, V; ret_policy = false) = T!(cp, V, similar(V); ret_policy = ret_policy)

get_greedy!(cp, V, out) = T!(cp, V, out, ret_policy = true)
get_greedy(cp, V) = T(cp, V, ret_policy = true)

function K!(cp, c, out)
    (; R, Π, β, b, agrid, zgrid) = cp
    z_idx = 1:length(zgrid)
    γ = R*β

    #policy function | z = z_i
    cf = LinearInterpolation((agrid, z_idx), c; extrapolation_bc = Interpolations.Flat())

    #lower bound for optimization
    opt_lb = 1e-8

    for (i_z, z) in enumerate(zgrid)
        for (i_a, a) in enumerate(agrid)
            function h(t)
                cps = cf.(R * a + z - t, z_idx)
                expectation = dot(du.(cps), Π[i_z, :])
                return abs(du(t) - max(γ * expectation, du(R * a + z + b)))
            end
            opt_ub = R * a + z + b
            res = optimize(h, min(opt_lb, opt_ub - 1e-2), opt_ub, method = Optim.Brent())
            out[i_a, i_z] = minimizer(res)
        end
    end
    return out
end

K(cp, c) = K!(cp, c, similar(c))

function initialize(cp)
    (; R, β, b, agrid, zgrid) = cp
    shape = length(agrid), length(zgrid)
    V, c = zeros(shape...), zeros(shape...)

    for (i_z, z) in enumerate(zgrid)
        for (i_a, a) in enumerate(agrid)
            c_max = R * a + z + b
            c[i_a, i_z] = c_max
            V[i_a, i_z] = u(c_max)/(1-β)
        end
    end

    return V, c
end

cp = consumer_problem()
N = 80

V, c = initialize(cp)
println("Starting value function iteration")
for i in 1:N
    V = T(cp, V)
end
c1 = T(cp, V, ret_policy = true)

V2, c2 = initialize(cp)
println("Starting policy function iteration")
for i in 1:N
    c2 = K(cp, c2)
end

plot(cp.asset_grid, c1[:, 1], label = "value function iteration")
plot!(cp.asset_grid, c2[:, 1], label = "policy function iteration")
plot!(xlabel = "asset level", ylabel = "Consumption (low income)")
plot!(legend = :topleft)