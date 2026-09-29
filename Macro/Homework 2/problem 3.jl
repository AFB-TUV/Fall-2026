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

Π_z = rouwenhorst_Π(5,0.99)

z_grid = [0.93771, 0.96835, 1.0, 1.03268, 1.06643]

function draw_One(cdf)
    n = rand()
    for i in 1:length(cdf)
        if n <= cdf[i]
            return i
        end
    end
    return length(cdf)
end

zvec = zeros(Int, 55000)

function sim_z(zs, Π)
    zs[1] = 3
    for i in 2:length(zs)
        cdf = cumsum(Π[zs[i-1],:])
        zs[i] += draw_One(cdf)
    end
    return zs
end

function surplus(b, γ, zs, κ, β, η, σ, χ, Π, θs, Ss)
    Ss_new = zeros(length(Ss))
    for i in 1:length(zs)
        Ss_new[i] = zs[i] - b + β*(1 - σ - γ*λws(η, χ, θs[i]) * dot(Π[i,:],Ss))
    end
    return Ss_new
end

function entry(γ, zs, κ, β, η, χ, Π, Ss)
    θs_new = zeros(length(Ss))
    for i in 1:length(θs_new)
        θs_new[i] = ((β * χ * (1-γ) * dot(Π[i,:],Ss))/κ)^(1/η)
    end
    return θs_new
end

function λws(η, χ, θ)
    return χ*θ^(1-η)
end

function u_evo(σ, η, χ, θ, u)
    return u + σ * (1 - u) - λws(η, χ, θ) * u
end

z_temp = sim_z(zvec, Π_z)
z_sim = zeros(Int, length(z_temp) - 5000)
for i in 1:(length(z_temp) - 5000)
    z_sim[i] = z_temp[i + 5000]
end

zs_plot = plot(z_sim, legend = false)
png(zs_plot, "zs_plot")


function stss(b, γ, zs, κ, β, η, σ, χ, Π, θs = 5*ones(5), Ss = 5*ones(5), tol = 1e-6, maxits = 1000)
    local Ss_new
    local θs_new

    for i in 1:maxits
        θs_new = entry(γ, zs, κ, β, η, χ, Π, Ss)
        Ss_new = surplus(b, γ, zs, κ, β, η, σ, χ, Π, θs, Ss)
        #println("iteration ", i, " ", θs_new, " ", Ss_new)
        if maximum(abs.(Ss_new - Ss)) < tol
            return Ss_new, θs_new
        else
            θs = θs_new
            Ss = Ss_new
        end
   8 end
    println("Failed to converge!")    
    return 0
end

function simulate(simz, zs, θs, u, b, γ, κ, η, σ, χ)
    simθ = zeros(length(simz))
    simu = zeros(length(simz))
    simw = zeros(length(simz))

    #solve for period 1 specially since we have to manually designate u0
    simθ[1] = θs[simz[1]]
    simw[1] = wc(b, γ, zs[simz[1]], κ, simθ[1])
    simu[1] = u_evo(σ, η, χ, θs[simz[1]], u)

    for i in 2:length(simz)
        #println("one")
        simθ[i] = θs[simz[i]]
        #println("two")
        simw[i] = wc(b, γ, zs[simz[i]], κ, simθ[i])
        #println("three")
        simu[i] = u_evo(σ, η, χ, θs[simz[i]], simu[i-1])
        #println("four")
    end
    return simθ, simw, simu
end

Ss_star, θs_star = stss(b_p1, γ_p1, z_grid, κ_p1, β_p1, η_p1, σ_p1, χ_p1, Π_z)
ws_star = wc.(b_p1, γ_p1, z_grid, κ_p1, θs_star)
us_star = u.(σ_p1, θs_star .* λ.(η_p1, χ_p1, θs_star))

θ_sim, w_sim, u_sim = simulate(z_sim, z_grid, θs_star, mean(us_star), b_p1, γ_p1, κ_p1, η_p1, σ_p1, χ_p1)
v_sim = θ_sim .* u_sim

sdz = std(log.(z_sim))
sdθ = std(log.(θ_sim))
sdu = std(log.(u_sim))
sdw = std(log.(w_sim))
r_θz = sdθ/sdz
r_wz = sdw/sdz
ρ_uv = cor(log.(u_sim), log.(v_sim))

#part c

function J_max(zs, β, σ, w, Π, Js)
    Js_new = zeros(length(Js))
    for i in 1:length(zs)
        Js_new[i] = zs[i] - w + β*(1-σ)*dot(Π[i,:], Js)
        Js_new[i] = max(0.0, Js_new[i])
    end
    return Js_new
end

function solve_J(zs, β, σ, w, Π, Js, tol = 1e-6, maxits = 1000)
    Js_new = zeros(length(Js))
    for i in 1:maxits
        #println(Js_new)
        Js_new = J_max(zs, β, σ, w, Π, Js)
        if maximum(abs.(Js_new - Js)) < tol
            return Js_new
        end
        Js = Js_new
    end
    println("Failed to converge!")    
    return 0
end


w_rigid = γ_p1 * (z_p1 + κ_p1 * θ_star) + (1 - γ_p1) * b_p1
ϵ_pred = (1/(1-w_rigid))/η_p1
Js_star = solve_J(z_grid, β_p1, σ_p1, w_star, Π_z, ones(length(z_grid)))
θs_J = entry(γ_p1, z_grid, κ_p1, β_p1, η_p1, χ_p1, Π_z, Js_star)
us_J = u.(σ_p1, θs_J .* λ.(η_p1, χ_p1, θs_J))

θ_Jsim, w_Jsim, u_Jsim = simulate(z_sim, z_grid, θs_J, mean(us_star), b_p1, γ_p1, κ_p1, η_p1, σ_p1, χ_p1)
v_Jsim = θ_Jsim .* u_Jsim

sdz_J = std(log.(z_sim))
sdθ_J = std(log.(θ_Jsim))
sdu_J = std(log.(u_Jsim))
r_θz = sdθ_J/sdz_J
ρ_uv = cor(log.(u_Jsim), log.(v_Jsim))

zls = zeros(length(z_sim))
for i in 1:length(z_sim)
    if z_sim[i] == 1
        zls[i] += 1.0
    end
end
zl_frac = mean(zls)