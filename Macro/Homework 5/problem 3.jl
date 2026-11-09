using Plots, LinearAlgebra, Random, Distributions

Random.seed!(2191997)

function sim_ϕ(μϕ, ρ, ϕgrid, Nϕ, N_sim, ϕ0, μϵ, σϵ)
    Dϵ = Normal(μϵ,σϵ)
    
    #eliminate influence of starting guess
    for i in 1:5000
        ϕ0 = ϕ_evo(μϕ, ρ, ϕ0, rand(Dϵ))
    end

    ϵvec = rand(Dϵ, N_sim)
    ϕvec = zeros(N_sim)
    ϕindices = zeros(Int, N_sim)

    ϕvec[1] = ϕ_evo(μϕ, ρ, ϕ0, ϵvec[1])
    ϕindices[1] = findIndex(ϕvec[1], ϕgrid)
    ϕmat = zeros(Nϕ, Nϕ)

    for i in 2:N_sim
        ϕvec[i] = ϕ_evo(μϕ, ρ, ϕvec[i-1], ϵvec[i])
        ϕindices[i] = findIndex(ϕvec[i], ϕgrid)
        ϕmat[ϕindices[i-1], ϕindices[i]] += 1.0
    end

    for i in 1:Nϕ
        if sum(ϕmat[:,i]) > 0.0
            ϕmat[:,i] = ϕmat[:,i] ./ sum(ϕmat[:,i])
        end
    end

    ϕπ = stationary_markov(ϕmat)
    #return ϕmat, ϕπ
    return ϕmat, ϕπ, ϕvec
end

function findIndex(ϕ, vec)
    ϕret = searchsortedlast(vec, ϕ)
    if ϕret == 0
        ϕret = length(vec)
    end
    return ϕret
end

function D_entrants(μ, σ, ϕgrid)
    D = cdf.(Normal(μ, σ), ϕgrid)
    for i in 2:length(D)
        D[i] -= D[i-1]
    end
    return D
end

#finds a stationary distribution of a Markov chain
function stationary_markov(Π,tol=1E-6)
    n = size(Π)[1]
    π = fill(1/n,n)

    for i in 1:10_000
        π_new = Π' * π
        if maximum(abs.(π_new - π)) < tol
            return π_new
        end
        π = π_new
    end
end

function ϕ_evo(μϕ, ρ, ϕ, ϵ)
    return (1 - ρ)*μϕ + ρ*ϕ + ϵ
end

Πtest, πtest = sim_ϕ(log(1.5), 0.9, range(-2.0, 2.0, 401), 501, 1_000_000, mean(LogNormal(0.0, 0.15)), 0.0, 0.15)

#commented out income process plotting because it is computer resource-intensive
#simplot = plot(ϕvec)
#png(simplot, "sim_plot")

function model(; σ = 4.0, β = 0.95, w = 1.0, cf = 0.25, ce = 0.25, Lbar = 1.0,
                ρ = 0.9, μϕ = log(1.5), μϵ = 0.0, σϵ = 0.15, N_sim = 10_000_000, 
                ϕmin = -2.0, ϕmax = 2.0, Nϕ = 401,
                μe = 0.0, σe = 0.1 )
    ϕgrid = range(ϕmin, ϕmax, Nϕ)
    πe = D_entrants(0.0, 0.1^0.5, ϕgrid)
    Π, π = sim_ϕ(μϕ, ρ, ϕgrid, Nϕ, N_sim, mean(LogNormal(μe, σϵ)), μϵ, σϵ)
    ϕgrid = exp.(ϕgrid)
    return (; σ, β, w, cf, ce, Lbar, ρ, μϕ, μϵ, N_sim, ϕmin, ϕmax, Nϕ, μe, σϵ, πe, Π, π, ϕgrid)
end

#series of functions that define profits through building up the CES theory
function price(ϕ, σ, w)
    return (σ/(σ-1))*(w/ϕ)
end

function quantity(A, ϕ, Lbar, w, σ)
    return w*Lbar*A*price(ϕ, σ, w)^(-σ)
end

function revenue(A, ϕ, Lbar, w, σ)
    return quantity(A, ϕ, Lbar, w, σ)*price(ϕ, σ, w)
end

function profit(A, ϕ, Lbar, w, σ, cf)
    return revenue(A, ϕ, Lbar, w, σ)/σ - w*cf
end


#iterate the value function to convergence
function vf_ss(A, model; tol = 1E-5, maxits = 10_000)
    (; Lbar, w, σ, ϕgrid, cf) = model
    V_old = profit.(A, ϕgrid, Lbar, w, σ, cf)
    W = similar(V_old)
    for i in 1:maxits
        V, W = backward_iteration(A, V_old, model)
        if i % 10 == 0
            D = maximum(abs.(V - V_old))
            #println("Iteration ", i, "Difference: ", D)
            if D < tol
                χ = zeros(length(W))
                χ[W .> 0.0] .= 1.0
                return V, W, χ
            end
        end
        V_old = V
    end    
end

#Bellman operator
function backward_iteration(A, V_old, model)
    (; β, Lbar, w, σ, ϕgrid, cf, Π) = model

    W = β .* Π * V_old
    V_new = profit.(A, ϕgrid, Lbar, w, σ, cf) + max.(0,W)
    
    return V_new, W
end

#find the A that pins down the stationary distribution
function solve_A(model)
    (; ϕgrid, w, ce, π) = model

    function solveBell(Aguess)
        V, W, χ = vf_ss(Aguess, model)
        return dot(π,W) - (w * ce) 
    end

    A = bisection(solveBell)
    V, W, χ = vf_ss(A, model)
    ϕ = ϕgrid[searchsortedlast(χ, 0.0) + 1]

    return V, W, χ, A, ϕ
end

#bisection algorithm
function bisection(f; A_a = 1e-6, A_b = 5.0, tol = 1e-8, maxits = 10_000)
    fa, fb = f(A_a), f(A_b)
    i = 0
    local fc
    for i in 1:maxits
        A_c = (A_a +A_b)/2
        fc = f(A_c)
        if abs(A_b - A_a) < tol
            return A_c
        elseif fa*fc > 0
            A_a = A_c
            fa = fc
        else
            A_b = A_c
            fb = fc
        end
    end
    println("Maximum iterations exceeded!")
end

#solve law of motion for fixed size of M 
function entryTest(Π, μ, χ, πe, M; maxits = 10_000, tol = 1e-6)
    μ_prime = zeros(length(μ))
    for i = 1:maxits
        if i % 10 == 0
            println("Iteration ", i, " : ", μ_prime)
        end

        μ_prime = (χ .* Π)*μ + M * πe
        if maximum(abs.(μ_prime - μ)) < tol
            return μ_prime
        end
        μ = μ_prime
    end
    println("Firm distribution failed to converge :(")
end

function eq_labor(A, χ, model)
    (; π, πe, ce, cf, ϕgrid, Lbar, w, σ, Π) = model
    n = quantity.(A, ϕgrid, Lbar, w, σ) ./ ϕgrid
    μ = π

    function entry(M; maxits = 10_000, tol = 1e-6)
        μ_prime = zeros(length(μ))
        for i = 1:maxits
            μ_prime = (χ .* Π)*μ + M * πe
            if i % 10 == 0
                if maximum(abs.(μ_prime - μ)) < tol
                    return μ_prime
                end
            end
            μ = μ_prime
        end
        println("Firm distribution failed to converge :(")
    end
    
    function L_clear(M)
        return M*ce + dot(n .+ cf, entry(M)) - Lbar
    end
    M = bisection(L_clear)
    μ = entry(M)
    return M, μ, n
end

function Hopenhayn(model)
    V, W, χ, A, ϕ = solve_A(model)
    M, μ, n = eq_labor(A, χ, model)
    return (; V, W, χ, A, ϕ, M, μ, n)
end

#=
m = model()
V1, W1, χ1, A1, ϕ1 = solve_A(m)
μ1 = entryTest(m.Π, m.π, χ1, m.π, 0.2)
M = eq_labor(A1, χ1, m)

sol1 = Hopenhayn(m)

πvec = profit.(1.0, m.ϕgrid, m.Lbar, m.w, m.σ, m.cf)


πplot = plot(m.ϕgrid, πvec)
vplot = plot(m.ϕgrid, sol1.V)
wplot = plot(m.ϕgrid, sol1.W)

png(πplot, "prof_plot")
png(vplot, "v_plot")
png(wplot, "w_plot")
=#