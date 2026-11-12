using Plots, LinearAlgebra

function ϵ_D(s::Float64, σ::Float64, ρ::Float64)
    return ((1 - s)/σ + s/ρ)^(-1)
end

function markup(ϵ::Float64)
    return ϵ/(ϵ - 1)
end

function price(σ::Float64, ϕ::Float64, τ::Float64)
    return τ*(σ/(σ - 1))/ϕ
end

function quantity(p::Float64, P::Float64, σ::Float64, ρ::Float64, E::Float64)
    return (p^(-σ))*(P^(σ - ρ))*E
end

function agg_p(p::Matrix{Float64}, σ::Float64)
    return (sum(p .^ (1 - σ)))^(1/(1 - σ))
end

function shares(p, q, D)
    return (p*q)/D
end

function model(; E_H = 100.0, E_F = 100.0, ϕ_H = [1.5 1.35 1.25 1.2], ϕ_F = [1.55 1.5 1.45 1.25],
     ρ = 4.0, σ = 10.0, τ_f = 1.2)
     ϕ = [ ϕ_H; ϕ_F ]
     E = [E_H; E_F]
     τ = [1.0; τ_f]
    return (; E_H, E_F, E, ϕ_H, ϕ_F, ϕ, ρ, σ, τ)
end

function find_shares(model; maxits = 10_000, tol = 1e-20)
    (; E_H, E_F, E, ϕ_H, ϕ_F, ϕ, ρ, σ, τ) = model
    s = 0.2*ones(Float64, size(ϕ))  
    s_new = zeros(Float64, size(ϕ))
    ϵ_demand = zeros(Float64, size(ϕ))
    markups = zeros(Float64, size(ϕ))
    q = zeros(Float64, size(ϕ))
    p = zeros(Float64, size(ϕ))
    P = 0.0
    D = 0.0

    for i in 1:maxits
        #println("Iteration ", i)
        ϵ_demand = ϵ_D.(s, σ, ρ)
        markups = markup.(ϵ_demand)
        p = price.(σ, ϕ)
        P = agg_p(p, σ)
        q = quantity.(p, P, σ, ρ, E)
        D = dot(p,q)
        s_new = shares.(p, q, D)

        if i % 10 == 0
            if maximum(abs.(s - s_new)) < tol
                return (; s = s_new, p = p, q = q, P = P, μ = markups, ϵ = ϵ_demand)
            end
        end
        s = s_new
    end
end

function HHI(shares::Matrix{Float64})
    return sum(shares .^ 2)
end

function cstatics_trade(τ_i::Float64, m::Function)
    model = m(τ = τ_i)
    sol = find_shares(model)
    h = HHI(sol.s)
    return h
end

#part 1
m = model()
sol1 = find_shares(m)

#part 2

HHI1 = HHI(sol1.s)

#part 3

τvec = range(1.0, 2.0, 101)
HHIvec = zeros(length(τvec))

for i in 1:length(τvec)
    #println("Trade cost: ", τvec[i])
    HHIvec[i] = cstatics_trade(τvec[i], model)
end

HHIplot = plot(HHIvec)
