using Plots, LinearAlgebra

function model(; τ0 = 0.2, τ = 0.0, β = 0.96, θ = 0.36, T = 60)
    s = β * θ * (1 - τ)
    s0 = β * θ * (1 - τ0)
    k0 = k_ss(s0, θ)
    kT = k_ss(s, θ)
    c0 = c(s0, k0, θ)
    cT = c(s, kT, θ)
    return (; τ0, τ, β, θ, T, s, s0, k0, kT, c0, cT)
end

function k_evo(k, s, θ)
    return s*(k^θ)
end

function k_ss(s, θ)
    return s^(1/(1-θ))
end

function c(s,k,θ)
    return (1 - s)*k^θ
end

function r(θ, k)
    return θ * k^(θ - 1)
end

function euler(c, β, τ, r)
    return c/(β * (1 - τ) * r)
end

function constraint(k,c,θ)
    return k^θ - c
end

p2 = model()

#part 1

#part 2

kvals_2 = zeros(p2.T + 1)
cvals_2 = zeros(p2.T + 1)
kvals_2[1] = p2.k0

#set up consumption to show impact response
cvals_2[1] = c(p2.s0, kvals_2[1], p2.θ)

for i in 1:p2.T
    kvals_2[i+1] = k_evo(kvals_2[i], p2.s, p2.θ)
    cvals_2[i+1] = kvals_2[i]^p2.θ - kvals_2[i+1]
end

imrs_c = (cvals_2[2] - cvals_2[1])/cvals_2[1]

kplot_2 = plot(kvals_2, title = "Capital", legend = false)
cplot_2 = plot(cvals_2, title = "Consumption", legend = false)
png(kplot_2, "k_plot_cf")
png(cplot_2, "c_plot_cf")


#part 3


function transition_path(model; ξ = 0.1, mink = 1e-6, maxits = 100_000, tol = 1e-10, noise = false)
    (; τ0, τ, β, θ, T, s, s0, k0, kT, c0, cT) = model
    
    ks = Vector(range(k0, kT, T+1))
    ks_up = zeros(T+1)    
    ks_up[1], ks_up[T+1] = k0, kT
    cs = zeros(T+1)
    cs[1], cs[T+1] = c0, cT
    rs = r.(θ, ks)
    its = maxits
    
    for i in 1:maxits
        if noise == true
           println("Iteration ", i)    
        end
        #iterate holding c_0 fixed to show impact response
        for j in 1:T-1
            cs[T+1-j] = euler(cs[T+2-j], β, τ, rs[T+2-j])
        end

        for j in 1:T-1
            ks_up[j+1] = max(constraint(ks[j], cs[j], θ), mink)
        end
        
        if maximum(abs.(ks_up - ks)) < tol
            its = i
            break
        end
        
        ks = ξ*ks_up + (1.0 - ξ)*ks
        rs = r.(θ, ks)
    end

    #set up array to return k_{t+1} and not k_t 
    ks_ret = zeros(length(ks_up))
    ks_ret[length(ks_ret)] = kT
    for j in 1:T
        ks_ret[j] = ks_up[j+1]
    end

    return ks_ret, cs, its
end

k_iter, c_iter, numits = transition_path(p2)
kplot_3 = plot([k_iter kvals_2], title = "Capital", label = ["Closed-Form" "Computed"])
cplot_3 = plot([c_iter cvals_2], title = "Consumption", label = ["Closed-Form" "Computed"])
png(kplot_3, "k_plot_iter")
png(cplot_3, "c_plot_iter")

cdiff = maximum(abs.(c_iter - cvals_2))
kdiff = maximum(abs.(k_iter - kvals_2))

k_iter_2, c_iter_2, numits_2 = transition_path(p2; ξ = 0.1, mink = -Inf)

k_iter_3, c_iter_3, numits_3 = transition_path(p2; ξ = 1.0, mink = 1e-6)
kplot_4 = plot([k_iter_3 kvals_2], title = "Capital", label = ["Closed-Form" "Computed"])
cplot_4 = plot([c_iter_3 cvals_2], title = "Consumption", label = ["Closed-Form" "Computed"])
png(kplot_4, "k_plot_bad")
png(cplot_4, "c_plot_bad")

#this line fails; I include it for the sake of reference but comment it out accordingly
#k_iter_4, c_iter_4, numits_4 = transition_path(p2; ξ = 1.0, mink = -Inf, noise = true)
