using Plots, LinearAlgebra

#set up parameters

function model(; σ = 4.0, 
                β = 0.95, 
                w = 1.0, 
                cf = 0.25, 
                ce = 0.25,
                ρ = 0.9,
                σϵ = 0.15,
                μ = log(1.5),
                Lbar = 1.0,
                n = 9)

    zgrid, π, Π = discretize_productivity(ρ, σϵ^.5, n)
    return (; σ, β, w, cf, ce, ρ, σϵ, μ, Lbar, n, zgrid, π, Π)
end

#Implement the Rouwenhorst method; ends up not being used
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
    (; β, Lbar, w, σ, zgrid, cf, Π) = model
    V_old = profit.(A, zgrid, Lbar, w, σ, cf)
    W = similar(V_old)
    for i in 1:maxits
        V, W = backward_iteration(A, V_old, model)
        if i % 10 == 0
            D = maximum(abs.(V - V_old))
            println("Iteration ", i, "Difference: ", D)
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
function bisection(f; A_a = 0.0000001, A_b = 5.0, tol = 1e-8, maxits = 10_000)
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
    return 0
end

#solve the full model!
function Hopenhayn(model)
    
end

m = model()
V1, W1, χ1, A1, ϕ1 = solve_A(m)

πvec = profit.(1.0, m.zgrid, m.Lbar, m.w, m.σ, m.cf)
πplot = plot(m.zgrid, πvec)
vplot = plot(m.zgrid, V1)
wplot = plot(m.zgrid, W1)

png(πplot, "prof_plot")
png(vplot, "v_plot")
png(wplot, "w_plot")