using LinearAlgebra, BasicInterpolators, Plots, Statistics

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

function discretize_productivity(ρ,σ,n_e)
    p = (1+ρ)/2
    nshift = (n_e - 1)/2
    e = range(-nshift,nshift,n_e)
    α = 2*σ/sqrt(n_e-1)
    e = α .* e

    Π = rouwenhorst_Π(n_e, p)
    π = stationary_markov(Π)

    y = exp.(e)
    return y,π,Π
end


function backward_iteration(Va, model, B)

    (; β, Π, w, ν, r, a_grid, z_grid) = model
    Wa = (β .* Π) * Va

    c_endog = Wa .^ -1.0
    h_endog = ((w .* z_grid ./ B) .* Wa).^ν
    coh = w .* z_grid .* h_endog .+ (1.0 + r) .* a_grid'
        #println(size(coh))
    a = zeros(size(coh))
    for e in 1:length(z_grid)
        itp_c = LinearInterpolator(c_endog[e,:] + a_grid, a_grid,NoBoundaries())
        a[e,:] = itp_c.(coh[e,:])
    end

    a = max.(a,a_grid[1])
    c = coh - a
    h = ((w .* z_grid ./ B) .* (c .^ -1)) .^ ν
    Va = (1+r) .* (c .^ -1)

    return Va, a, c, h
end

function policy_ss(model, B, tol = 1E-9)

    (; w, r, a_grid, z_grid) = model
    coh = w .* z_grid .+ (1.0 + r) .* a_grid'
    c = 0.05 .* coh
    Va = (1.0+r) ./ c
    a_old = zeros(size(coh))
    for i in 1:10_000
        Va, a, c, h = backward_iteration(Va, model, B)
        if i > 1 && maximum(abs.(a - a_old)) < tol
            return Va, a, c, h
        end
        a_old = a
    end
end

function get_lottery(a, a_grid)
    a_i = searchsortedlast(a_grid, a)
    a_π = (a_grid[a_i + 1] - a)/(a_grid[a_i + 1] - a_grid[a_i])
    if a_π == 1.0
        a_π = 0.0
    end
    return a_i, a_π
end

function forward_policy(D, a_i, a_π)
    Dend = zeros(size(D))
    for i in 1:size(a_i)[1]
        for j in 1:size(a_i)[2]
            Dend[i,a_i[i,j]] += a_π[i,j]*D[i,j]
            Dend[i,a_i[i,j]+1] += (1-a_π[i,j])*D[i,j]
        end
    end
    return Dend
end

function forward_iteration(D,Π,a_i,a_π)
    Dend = forward_policy(D, a_i, a_π)
    return Π' * Dend
end

function distribution_ss(model, a, tol = 1E-6)
    (; a_grid, Π) = model
    a_i = Matrix{Int64}(undef, size(a))
    a_π = zeros(size(a))
    for i in 1:size(a)[1]
        for j in 1:size(a)[2]
            a_i[i,j], a_π[i,j] = get_lottery(a[i,j], a_grid)
        end
    end

    π = stationary_markov(Π)
    D = π*ones(size(a_grid))' ./ length(a_grid)

    for i in 1:10_000
        #println("iteration ", i)
        D_new = forward_iteration(D, Π, a_i, a_π)
        #println(maximum(abs.(D_new - D)))
        if maximum(abs.(D_new - D)) < tol
            return D_new
        end
        D = D_new
    end
end

function household(model, B)
    Va, a, c, h = policy_ss(model, B)
    D = distribution_ss(model, a)
    return Va, a, c, h, D
end

function calib(model, B)
    Va, a, c, h, D = household(model, B)
    return dot(h,D) - (1/3)
end
    

function calibrate_B(model; B_a = 0.4, B_b = 30.0, tol = 1e-8, maxits = 10_000)
    fa, fb = calib(model, B_a), calib(model, B_b)
    i = 0
    local B_c
    for i in 1:maxits
        B_c = (B_a + B_b)/2
        fc = calib(model, B_c)
        if abs(B_b - B_a) < tol
            return B_c
        elseif fa*fc > 0
            B_a = B_c
            fa = fc
        else
            B_b = B_c
            fb = fc
        end
    end
    println("Maximum iterations exceeded!")
    return 0
end

function P2(; β = 0.94, r = 0.04, w = 1.0, ν = 0.5, ρ = 0.95, σ = 0.4, n_e = 7, amin = 0.0, amax = 50.0, n_a = 1001)
    z_grid, π, Π = discretize_productivity(ρ, σ, n_e)
    a_grid = range(amin, amax, n_a)
    return (; β, r, w, ν, ρ, σ, n_e, z_grid, π, Π, a_grid, amin, amax, n_a)
end

p2 = P2()
B_p2 = calibrate_B(p2)

Va_p2, a_p2, c_p2, h_p2, D_p2 = household(p2, B_p2)

hbar_p2 = dot(h_p2, D_p2)
L_p2 = dot(p2.z_grid .* h_p2,D_p2)

#part 1

#[1] = 0.0, [21] = 1.0, [41] = 2.0, [101] = 5.0, [201] = 10.0
h_out_2 = [h_p2[4,1], h_p2[4,21], h_p2[4,41], h_p2[4,101], h_p2[4,201]]
corrs = zeros(length(p2.z_grid))
for (i,z) in enumerate(p2.z_grid)
    corrs[i] = cor(h_p2[i,:], a_p2[i,:])
end

#part 2

#part 3

p2c = P2(w = 1.02)
Va_p2c, a_p2c, c_p2c, h_p2c, D_p2c = household(p2c, B_p2)

L_impact = dot(p2.z_grid .* h_p2c,D_p2)
ΔL = (L_impact - L_p2)/L_p2
ϵL = ΔL/(p2c.w - p2.w)

println("Parts I through III done!")

#part 4

function backward_iteration_d(Va, model)

    (; β, Π, w, ν, r, a_grid, z_grid, hbar, B) = model
    Wa = (β .* Π) * Va

    c_endog = Wa .^ -1.0
    coh_e = w .* z_grid .* hbar .+ (1.0 + r) .* a_grid'
    coh_n = ones(length(z_grid)) .* (1.0 + r) .* a_grid'
    a_e = zeros(size(coh_e))
    a_n = zeros(size(coh_e))
    for i in 1:size(coh_e)[1]
        itp_e = LinearInterpolator(c_endog[i,:] + a_grid, a_grid,NoBoundaries())
        for j in 1:size(coh_e)[2]
            a_e[i,j] = itp_e(coh_e[i,j])
            a_n[i,j] = itp_e(coh_n[i,j])
        end
    end

    a_e = max.(a_e,a_grid[1])
    a_n = max.(a_n,a_grid[1])
    c_e = coh_e - a_e .+ B
    c_n = coh_n - a_n
    h_e = (1/3) .* ones(size(coh_e))
    h_n = 0.0 .* ones(size(coh_n)) 
    Va_e = (1+r) .* (c_e .^ -1)
    Va_n = (1+r) .* (c_n .^ -1)

    Va, a, c, h = zeros(size(coh_e)), zeros(size(coh_e)), zeros(size(coh_e)), zeros(size(coh_e))
    #we want to pick the higher value function => should pick the lower marginal utility
    for i in 1:length(z_grid)
        for j in 1:length(a_grid)
            if Va_e[i,j] < Va_n[i,j]
                Va[i,j] = Va_e[i,j]
                a[i,j]  = a_e[i,j]
                c[i,j] = c_e[i,j]
                h[i,j] = h_e[i,j]
            else
                Va[i,j] = Va_n[i,j]
                a[i,j] = a_n[i,j]
                c[i,j] = c_n[i,j]
                h[i,j] = h_n[i,j]
            end
        end
    end
    return Va, a, c, h
end

function policy_ss_d(model, tol = 1e-9)
    (; w, r, a_grid, z_grid) = model
    coh = w .* z_grid .+ (1.0 + r) .* a_grid'
    c = 0.05 .* coh
    Va = (1.0+r) ./ c
    a_old = zeros(size(coh))
    for i in 1:10_000
        Va, a, c, h = backward_iteration_d(Va, model)
        if i > 1 && maximum(abs.(a - a_old)) < tol
            return Va, a, c, h
        end
        a_old = a
    end
end    

function P2d(; β = 0.94, r = 0.04, w = 1.0, ν = 0.5, B = 1.0, hbar = 1/3, ρ = 0.95, σ = 0.4, n_e = 7, amin = 0.0, amax = 50.0, n_a = 1001)
    z_grid, π, Π = discretize_productivity(ρ, σ, n_e)
    a_grid = range(amin, amax, n_a)
    return (; β, r, w, ν, ρ, σ, n_e, B, hbar, z_grid, π, Π, a_grid, amin, amax, n_a)
end

function household_d(model)
    Va, a, c, h = policy_ss_d(model)
    D = distribution_ss(model, a)
    return Va, a, c, h, D
end

p2d = P2d()
Va_p2d, a_p2d, c_p2d, h_p2d, D_p2d = household_d(p2d)
