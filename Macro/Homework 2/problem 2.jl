#Problem 2

function calib(b, γ, z, κ, r, η, σ, χ, θ, θa = 0.0, θb = 2.0)
    return equilibrium(0.0, 2.0, b, γ, z, κ, r, η, σ, χ) - θ
end

function calibrate_b(γ, z, κ, r, η, σ, χ, θ, b_a = 0.90, b_b = 0.98, tol = 1e-8, maxits = 1000)
    fa, fb = calib(b_a, γ, z, κ, r, η, σ, χ, θ), calib(b_b, γ, z, κ, r, η, σ, χ, θ)
    i = 0
    local b_c
    for i in 1:maxits
        b_c = (b_a+b_b)/2
        fc = calib(b_c, γ, z, κ, r, η, σ, χ, θ)
        if abs(b_b - b_a) < tol
            return b_c
        elseif fa*fc > 0
            b_a = b_c
            fa = fc
        else
            b_b = b_c
            fb = fc
        end
    end
    println("Maximum iterations exceeded!")
    return 0
end


function ϵθz(b, γ, z, κ, r, η, σ, χ, θ)
    return (z*(1 - γ)/θ)/((η*(r + σ)*κ/χ) *(θ^(η-1)) + (γ * κ))
end

γ_new = 0.05
b_new = calibrate_b(γ_new, z_p1, κ_p1, r_p1, η_p1, σ_p1, χ_p1, θ_star)
θ_new, w_new, λf_new, λw_new, u_new, S_new, J_new = ss(b_new, γ_new, z_p1, κ_p1, r_p1, η_p1, σ_p1, χ_p1)
ϵ_new = ϵθz(b_new, γ_new, z_p1, κ_p1, r_p1, η_p1, σ_p1, χ_p1, θ_new)
ϵ_old = ϵθz(b_p1, γ_p1, z_p1, κ_p1, r_p1, η_p1, σ_p1, χ_p1, θ_star)

ϵvec_old = zeros(length(bvec))
ϵvec_new = zeros(length(bvec))
θ_vec_old = zeros(length(bvec))
θ_vec_new = zeros(length(bvec))
u_vec_old = zeros(length(bvec))
u_vec_new = zeros(length(bvec))


for i in 1:length(bvec)
    θ_vec_old[i] = equilibrium(0.0, 5.0, bvec[i], γ_p1, z_p1, κ_p1, r_p1, η_p1, σ_p1, χ_p1)
    θ_vec_new[i] = equilibrium(0.0, 100.0, bvec[i], γ_new, z_p1, κ_p1, r_p1, η_p1, σ_p1, χ_p1)
    ϵvec_old[i] = ϵθz(bvec[i], γ_p1, z_p1, κ_p1, r_p1, η_p1, σ_p1, χ_p1, θ_vec_old[i])
    ϵvec_new[i] = ϵθz(bvec[i], γ_new, z_p1, κ_p1, r_p1, η_p1, σ_p1, χ_p1, θ_vec_new[i])
    u_vec_old[i] = u(σ_p1, θ_vec_old[i]* λ(η_p1, χ_p1, θ_vec_old[i]))
    u_vec_new[i] = u(σ_p1, θ_vec_old[i]* λ(η_p1, χ_p1, θ_vec_new[i]))
end

ϵplot = plot(bvec, [ϵvec_old ϵvec_new 19*ones(length(bvec)) ϵ_new*ones(length(bvec))], xlabel = "b", ylabel = "Elasticity", label = ["γ = 0.5" "γ = 0.05" "ϵ ≈ 19" "steady-state, γ = 0.5"])
ϵplot_zoom = plot(bvec[70:97], [ϵvec_old[70:97] ϵvec_new[70:97] 19*ones(length(bvec[70:97]))], xlabel = "b", ylabel = "Elasticity", label = ["γ = 0.5" "γ = 0.05" "ϵ ≈ 19"])
png(ϵplot, "elasticities")


function b_deriv(bh, bl, uh, ul)
    return (uh - ul)/(log(bh) - log(bl))
end

dubvec_old = zeros(96)
dubvec_new = zeros(96)

for i in 1:length(dubvec_old)
    dubvec_old[i] = b_deriv(bvec[i+1], bvec[i], u_vec_old[i+1], u_vec_old[i])
    dubvec_new[i] = b_deriv(bvec[i+1], bvec[i], u_vec_new[i+1], u_vec_new[i])
end

ϵbplot = plot(bvec[1:96], [dubvec_old dubvec_new 0.3*ones(96)], xlabel = "b", ylabel = "Elasticity", label = ["γ = 0.5" "γ = 0.05" "ϵ = 0.3"])
png(ϵbplot, "belasticities")
