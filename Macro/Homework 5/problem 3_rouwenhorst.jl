using Plots, LinearAlgebra, Random, Distributions

Random.seed!(2191997)

function sim_ϕ(μϕ, ρ, ϕgrid, Nϕ, N_sim, ϕ0, σϵ)
    Dϵ = Normal(0,σϵ)
    
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

    return ϕmat
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

function ϕ_evo(μϕ, ρ, ϕ, ϵ)
    return (1 - ρ)*μϕ + ρ*ϕ + ϵ
end

#=
function model(; σ = 4.0, β = 0.95, w = 1.0, cf = 0.25, ce = 0.25, Lbar = 1.0,
                ρ = 0.9, μϕ = log(1.5), σϵ = 0.15, N_sim = 100_000, ϕmin = -2.5, ϕmax = 2.5, Nϕ = 501,
                μe = 0.0, σe = 0.1 )
    ϕgrid = range(ϕmin, ϕmax, Nϕ)
    πe = D_entrants(0.0, 0.1^0.5, ϕgrid)
    Π = sim_ϕ(μϕ, ρ, ϕgrid, Nϕ, N_sim, mean(LogNormal(μe, σe)), σϵ)
    ϕgrid = exp.(ϕgrid)
    return (; σ, β, w, cf, ce, Lbar, ρ, μϕ, μϵ, N_sim, ϕmin, ϕmax, Nϕ, μe, σe, πe, Π, ϕgrid)
end
=#
#code to set up a model using an AR(1) discretized by a Rouwenhorst


function model(; σ = 4.0, 
                β = 0.95, 
                w = 1.0, 
                cf = 0.25, 
                ce = 0.25,
                ρ = 0.9,
                σϵ = 0.15,
                μϕ = log(1.5),
                Lbar = 1.0,
                n = 51)

    zgrid, π, Π = discretize_productivity(ρ, σϵ, n)
    return (; σ, β, w, cf, ce, ρ, σϵ, μϕ, Lbar, n, zgrid, π, Π)
end

#Implement the Rouwenhorst method
function rouwenhorst_Π(N, p)
    Π = [p 1-p; 1-p p]

    for i in 3:N
        Π_old = Π
        Π = zeros(i,i)

        Π[1:i-1, 1:i-1] = p*Π_old
        Π[1:i-1, 2:i] .+= (1-p)*Π_old
        Π[2:i, 1:i-1] .+= (1-p)*Π_old
        Π[2:i,2:i] .+= p*Π_old
        Π[2:i-1,:] ./= 2
    end
    
    return Π
end

#finds a stationary distribution of a Markov chain
function stationary_markov(Π,tol=1E-14)
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


#function to discretize an AR(1) 
function discretize_productivity(ρ,σ,n_e)
    p = (1+ρ)/2
    nshift = (n_e - 1)/2
    e = range(-nshift,nshift,n_e)
    α = 2*σ/sqrt(n_e-1)
    e = α .* e

    Π = rouwenhorst_Π(n_e, p)
    π = stationary_markov(Π)

    z = exp.(e)
    #z = e .+ 0.000001
    return z,π,Π
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
    (; Lbar, w, σ, zgrid, cf) = model
    V_old = profit.(A, zgrid, Lbar, w, σ, cf)
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
    (; β, Lbar, w, σ, zgrid, cf, Π) = model

    W = β .* Π * V_old
    V_new = profit.(A, zgrid, Lbar, w, σ, cf) + max.(0,W)
    
    return V_new, W
end

#find the A that pins down the stationary distribution
function solve_A(model)
    (; zgrid, w, ce, π) = model

    function solveBell(Aguess)
        V, W, χ = vf_ss(Aguess, model)
        return dot(π,W) - (w * ce) 
    end

    A = bisection(solveBell)
    V, W, χ = vf_ss(A, model)
    ϕ = zgrid[searchsortedlast(χ, 0.0) + 1]

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
        #if i % 10 == 0
        #    println("Iteration ", i, " : ", μ_prime)
        #end

        μ_prime = (χ .* Π)*μ + M * πe
        if maximum(abs.(μ_prime - μ)) < tol
            return μ_prime
        end
        μ = μ_prime
    end
    println("Firm distribution failed to converge :(")
end

function eq_labor(A, χ, model)
    (; π, ce, cf, zgrid, Lbar, w, σ, Π) = model
    n = quantity.(A, zgrid, Lbar, w, σ) ./ zgrid
    μ = π

    function entry(M; maxits = 10_000, tol = 1e-6)
        μ_prime = zeros(length(μ))
        for i = 1:maxits
            μ_prime = (χ .* Π)*μ + M * π
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

m = model()
μ = entryTest(m.Π, m.π, χ1, m.π, 0.2)
#M = eq_labor(A1, χ1, m)

#part 1
πvec = profit.(1.0, m.zgrid, m.Lbar, m.w, m.σ, m.cf)
πplot = plot(m.zgrid, πvec, xlabel = "Productivity", title = "Profits", legend = false)
png(πplot, "prof_plot")

#part 2
V, W, χ = vf_ss(1.0, m)
vplot = plot(m.zgrid, V, xlabel = "Productivity", title = "Value Function (A = 1)", legend = false)
png(vplot, "v_plot")

#part 3
wplot = plot(m.zgrid, W)
png(wplot, "w_plot")
ϕ = m.zgrid[searchsortedlast(χ, 0.0) + 1]

#part 4
V1, W1, χ1, A, ϕ1 = solve_A(m)

#part 5
μ = entryTest(m.Π, m.π, χ1, m.π, 0.2)
μplot = plot(m.zgrid, μ, xlabel = "Productivity", title = "Firm Size Distribution (M = 0.2)", legend = false)
png(μplot, "mu_plot")

#part 6 and 7

sol1 = Hopenhayn(m)
mass_1 = sum(sol1.μ)
ϕbar_1 = dot(m.zgrid, sol1.μ)
distplot = plot(m.zgrid, sol1.μ, xlabel = "Productivity", title = "Stationary Firm Size Distribution", legend = false)
png(distplot, "dist_plot")

#part 8

m2 = model(cf = 0.3)
sol2 = Hopenhayn(m2)
mass_2 = sum(sol2.μ)
ϕbar_2 = dot(m.zgrid, sol2.μ)
ΔA = (sol2.A - sol1.A)/sol1.A 
Δϕ = (sol2.ϕ - sol1.ϕ)/sol1.ϕ
ΔM = (sol2.M - sol1.M)/sol1.M 
Δmass = (mass_2 - mass_1)/mass_1
Δϕbar = (ϕbar_2 - ϕbar_1)/ϕbar_1