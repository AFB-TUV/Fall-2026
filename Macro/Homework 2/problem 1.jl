#Part A

β_p1 = 0.996
r_p1 = 1/β_p1 - 1
z_p1 = 1.0
b_p1 = 0.4
η_p1 = 0.5
γ_p1 = η_p1
σ_p1 = 0.02
χ_p1 = 0.53
κ_p1 = 0.75

function wc(b, γ, z, κ, θ)
    return (1 - γ)*b + γ*(z + κ * θ)
end

function jc(z, κ, r, η, σ, χ, θ)
    return z - (r + σ)*κ*(θ^η)/χ
end

function eq(b, γ, z, κ, r, η, σ, χ, θ)
    return z - (r+σ) * κ * (θ^η)/χ - (1 - γ) * b - γ * (z + κ*θ)
end

function λ(η, χ, θ)
    return χ*(θ^(-η))
end

function β(r)
    return 1/(1+r)
end

function S(b, γ, z, β, σ, λ)
    return (z - b) / (1 - β * (1 - σ) + β * γ * λ)
end

function J(γ, S)
    return γ * S
end

function u(σ, λ)
    return σ/(σ + λ)
end

function ss(b, γ, z, κ, r, η, σ, χ, θa = 0.0, θb = 20.0)
    θ_star = equilibrium(θa, θb, b, γ, z, κ, r, η, σ, χ)
    w_star = wc(b, γ, z, κ, θ_star)
    λf_star = λ(η, χ, θ_star)
    λw_star = θ_star * λf_star
    u_star = u(σ, θ_star * λf_star)
    S_star = S(b, γ, z, β(r), σ, λw_star)
    J_star = J(γ, S_star)
    return θ_star, w_star, λf_star, λw_star, u_star, S_star, J_star
end

function equilibrium(θa, θb, b, γ, z, κ, r, η, σ, χ, tol = 1e-5, maxits = 1000)
    fa, fb = eq(b, γ, z, κ, r, η, σ, χ, θa), eq(b, γ, z, κ, r, η, σ, χ, θb)
    i = 0
    local θc
    for i in 1:maxits
        θc = (θa+θb)/2
        fc = eq(b, γ, z, κ, r, η, σ, χ, θc)
        if θb - θa < tol
            return θc
        elseif fa*fc > 0
            θa = θc
            fa = fc
        else
            θb = θc
            fb = fc
        end
    end
    println("Maximum iterations exceeded!")
    return 0
end

function fe(κ, λf, β, γ, S)
    return κ/λf - β*(1-γ)*S
end

θ_star, w_star, λf_star, λw_star, u_star, S_star, J_star = ss(b_p1, γ_p1, z_p1, κ_p1, r_p1, η_p1, σ_p1, χ_p1)
fec = fe(κ_p1, λf_star, β(r_p1), γ_p1, S_star)

#Part B

bvec = range(0.0,0.96,97)
ustar_vec = zeros(length(bvec))
for i in 1:length(bvec)
    ustar_vec[i] = ss(bvec[i], γ_p1, z_p1, κ_p1, r_p1, η_p1, σ_p1, χ_p1)[5]
end
dmp_uplot = plot(bvec,ustar_vec, xlabel = "b", ylabel = "unemployment", legend = false)
png(dmp_uplot, "u_vs_b")