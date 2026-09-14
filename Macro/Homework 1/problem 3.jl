r_ss = 1/0.996 - 1
z_ss = 1.0
b_ss = 0.4
η_ss = 0.5
γ_ss = η_ss
σ_ss = s_bar
χ_ss = f_bar
κ_ss = (1-γ_ss)*(z_ss-b_ss)/(γ_ss+(r_ss+σ_ss)/χ_ss)

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
    S_star = S(b, γ, z, β(r), σ, λf_star)
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

function ϵθz(b, γ, z, κ, r, η, σ, χ, θ)
    return (z*(1 - γ)/θ)/((η*(r + σ)*κ/χ) *(θ^(η-1)) + (γ * κ))
end

θ_star, w_star, λf_star, λw_star, u_star, S_star, J_star = ss(b_ss, γ_ss, z_ss, κ_ss, r_ss, η_ss, σ_ss, χ_ss)
star_umd = 1/λw_star
star_vmd = 1/λf_star

#θvec = range(0.000001, 2.0, 500)
#jcvec = jc.(z_ss, κ_ss, r_ss, η_ss, σ_ss, χ_ss, θvec)
#wcvec = wc.(b_ss, γ_ss, z_ss, κ_ss, θvec)
#eqplot = plot(θvec, [jcvec wcvec])

#C2 

zl = 0.98
zh = 1.02
bh = 0.7
σh = maximum(s_vec)

θ_zl, w_zl, λf_zl, λw_zl, u_zl, S_zl, J_zl = ss(b_ss, γ_ss, zl, κ_ss, r_ss, η_ss, σ_ss, χ_ss)
θ_zh, w_zh, λf_zh, λw_zh, u_zh, S_zh, J_zh = ss(b_ss, γ_ss, zh, κ_ss, r_ss, η_ss, σ_ss, χ_ss)
θ_bh, w_bh, λf_bh, λw_bh, u_bh, S_bh, J_bh = ss(bh, γ_ss, z_ss, κ_ss, r_ss, η_ss, σ_ss, χ_ss)
θ_σh, w_σh, λf_σh, λw_σh, u_σh, S_σh, J_σh = ss(b_ss, γ_ss, z_ss, κ_ss, r_ss, η_ss, σh, χ_ss)

ϵ_zl = ϵθz(b_ss, γ_ss, zl, κ_ss, r_ss, η_ss, σ_ss, χ_ss, θ_zl)
ϵ_zh = ϵθz(b_ss, γ_ss, zh, κ_ss, r_ss, η_ss, σ_ss, χ_ss, θ_zh)


bvec = range(0.0,0.95,96)
ustar_vec = zeros(length(bvec))
θstar_vec = zeros(length(bvec))

for i in 1:length(bvec)
    ustar_vec[i] = ss(bvec[i], γ_ss, z_ss, κ_ss, r_ss, η_ss, σ_ss, χ_ss)[5]
    θstar_vec[i] = ss(bvec[i], γ_ss, z_ss, κ_ss, r_ss, η_ss, σ_ss, χ_ss)[1]
end
dmp_uplot = plot(bvec,ustar_vec)
ϵvec = ϵθz.(bvec, γ_ss, z_ss, κ_ss, r_ss, η_ss, σ_ss, χ_ss, θstar_vec)
dmp_ϵplot = plot(bvec,ϵvec)

zvec = zeros(38)
for i in 1:length(zvec)
    if 1 < i < 13
        zvec[i] = 1.02
    else
        zvec[i] = 1.0
    end
end

function u_evo(u, λw, σ)
    return (1 - λw)*u + σ*(1 - u)
end

function eq_evo(b, γ, zs, κ, r, η, σ, χ, u0)
    us = zeros(length(zs))
    us[1] = u0
    vs = zeros(length(zs))
    λws = zeros(length(zs))
    θs = zeros(length(zs))
    θs[1] = equilibrium(0.0, 20.0, b, γ, zs[1], κ, r, η, σ, χ)[1]
    vs[1] = θs[1]*us[1]
    λws[1] = λ(η, χ, vs[1]/us[1])
    
    for i in 2:length(zs)
        θs[i] = equilibrium(0.0, 20.0, b, γ, zs[i], κ, r, η, σ, χ)
        λws[i] = λ(η, χ, θs[i]) * θs[i]
        us[i] = u_evo(us[i-1], λws[i], σ)
        vs[i] = us[i] * θs[i]
    end
    bevs = zeros(length(zs))
    bevs = σ ./ (σ .+ λws)
    return us, vs, bevs, θs, λws
end

utrans, vtrans, bevcurve, θtrans, λtrans = eq_evo(b_ss, γ_ss, zvec, κ_ss, r_ss, η_ss, σ_ss, χ_ss, u_star)
transition_plot = plot([utrans[2:38] bevcurve[2:38]], vtrans[2:38], xlabel = "Unemployment", ylabel = "Vacancies", label = ["Unemployment Transition" "Beveridge Curve"])
transition_plot_alt = plot([utrans bevcurve], vtrans, xlabel = "Unemployment", ylabel = "Vacancies", label = ["Unemployment Transition" "Beveridge Curve"])

vtrans_plot = plot(vtrans[2:38], label = "Vacancies over time")
vtrans_plot_alt = plot(vtrans, label = "Vacancies over time")

firsttrans_diffs = utrans[1:12] .- minimum(utrans[1:12])
halftrans = maximum(firsttrans_diffs)/2
transline = halftrans .* ones(12)
diffsplot = plot([firsttrans_diffs transline], label = ["Gap from steady-state unemployment" "Halfway point"])

png(transition_plot, "mit_uv")
png(transition_plot_alt, "mit_uv_full")
png(vtrans_plot, "mit_vacancies")
png(vtrans_plot_alt, "mit_vacancies_full")
png(diffsplot, "mit_diffs")