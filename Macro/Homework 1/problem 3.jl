cd("C://Users/Alejandro/Documents/Actuarial Science/Fall 2026/Macro/Homework 1")
include("Problem 1.jl")
β_ss = 0.996
z_ss = 1.0
b_ss = 0.4
η_ss = 0.5
γ_ss = η_ss
σ_ss = s_bar
χ_ss = f_bar
λ_ss = χ_ss

function κ_cf(β, λ, γ, b, z, σ)
    return (β * σ * (1 - γ) * (z - b)/(1 - β*((1-σ) - γ * λ)))
end



function w(b, γ, z, κ, θ)
    return (1 - γ)*b + γ*(z + κ * θ)
end

function J(γ, S)
    return γ * S
end

function S(z, b, β, σ, γ, λ)
    return (z - b) / (1 - β * (1 - σ) + β * γ * λ)
end

function u(σ, λ)
    return σ/(σ + λ)
end

function ss(β, λ, γ, b, z, σ, θ)
    κ_ss = κ_cf(β_ss, λ_ss, γ_ss, b_ss, z_ss, σ_ss)
    w_ss = w(b_ss, γ_ss, z_ss, κ_ss, 1)
    u_ss = u(σ_ss, λ_ss)
    S_ss = S(z_ss, b_ss, β_ss, σ_ss, γ_ss, λ_ss)
    J_ss = J(γ_ss, S_ss)
    return κ_ss, w_ss, u_ss, S_ss, J_ss
end

κ_star, w_star, u_star, S_star, J_star = ss(β_ss, λ_ss, γ_ss, b_ss, z_ss, σ_ss, 1)
star_md = 1/λ_ss
