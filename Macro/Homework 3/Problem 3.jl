using LinearAlgebra, BasicInterpolators, Plots, Statistics

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

function backward_iteration(Va, model)

    (; β, Π, μ, r, y_grid, a_grid) = model
    Wa = (β .* Π) * Va

    c_endog = Wa .^ -μ
    coh = y_grid .+ (1.0 + r).*a_grid'
    a = zeros(size(coh))
    for e in 1:size(y_grid)[1]
        itp_c = LinearInterpolator(c_endog[e,:] + a_grid, a_grid,NoBoundaries())
        a[e,:] = itp_c.(coh[e,:])
    end

    a = max.(a,a_grid[1])
    c = coh - a

    Va = (1+r) .* c .^ (-1/μ)
    return Va, a, c
end

function policy_ss(model, tol = 1E-9)

    (; a_grid, r, μ, y_grid) = model
    coh = y_grid .+ (1.0+r) .* a_grid'
    c = 0.05 .* coh
    Va = (1+r) .* c .^ (-1/μ)
    #println(Va)
    a_old = zeros(size(coh))
    for i in 1:10_000
        Va, a, c = backward_iteration(Va, model)
        if i > 1 && maximum(abs.(a - a_old)) < tol
            return Va, a, c
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
        D_new = forward_iteration(D, Π, a_i, a_π)
        if maximum(abs.(D_new - D)) < tol
            return D_new
        end
        D = D_new
    end
end

function household(model)
    Va, a, c = policy_ss(model)
    D = distribution_ss(model, a)
    return (; Va, a, c, D)
end

function P3(; β = 0.996, r = 0.003, w = 1.0, μ = 2.0, λ = 0.45, σ = 0.02, b = 0.4, amin = 0.0, amax = 50.0, n_a = 1001)
    Π = [1.0 - σ σ; λ 1.0-λ]
    y_grid = [w; b]
    a_grid = range(amin, amax, n_a)
    return (; β, r, w, μ, σ, λ, Π, y_grid, a_grid, amin, amax, n_a)
end

p3 = P3()
sp3 = household(p3)

#part 1

#[1] = 0.0, [21] = 1.0, [41] = 2.0, [101] = 5.0, [201] = 10.0
s_p3 = sp3.a - ones(2)*p3.a_grid'

c_out_3 = [sp3.c[:,1], sp3.c[:,21], sp3.c[:,41], sp3.c[:,101], sp3.c[201]]
s_out_3 = [s_p3[:,1], s_p3[:,21], s_p3[:,41], s_p3[:,101], s_p3[:,201]]
abar = dot(sp3.a, sp3.D)
#meda = median(sp3.a .* sp3.D)

#part 2

#part 3