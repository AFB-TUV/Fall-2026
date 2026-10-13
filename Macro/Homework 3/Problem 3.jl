using LinearAlgebra, BasicInterpolators, Plots

w_p3 = 1.0
b_p3 = 0.4
μ_p3 = 2.0
β_p3 = 0.996
r_p3 = 0.003
λ_p3 = 0.45
σ_p3 = 0.02

Π_p3 = [1 - σ_p3 σ_p3; λ_p3 1 - λ_p3]
y_p3 = [1.0; 0.4]

function discretize_assets(amin,amax,n_a)
    ubar = log(1+log(1+amax-amin))
    ugrid = range(0,ubar,n_a)
    return amin .+ exp.(exp.(ugrid) .- 1) .- 1 
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

function backward_iteration(Va, Π, a_grid, y, r, β, σ)
    Wa = (β .* Π) * Va

    c_endog = Wa .^ -σ
    coh = broadcast(+,y,(1+r).*a_grid')

    a = zeros(size(coh))
    for e in 1:size(y)[1]
        #println(c_endog[e,:])
        itp_c = LinearInterpolator(c_endog[e,:] + a_grid, a_grid,NoBoundaries())
        a[e,:] = itp_c.(coh[e,:])
    end

    a = max.(a,a_grid[1])
    c = coh - a

    Va = (1+r) .* c .^ (-1/σ)
    return Va, a, c
end

function policy_ss(Π, a_grid, y, r, β, σ, tol = 1E-9)
    coh = y_p3 .+ (1+r_p3).*a_grid'
    c = 0.05 .* coh
    Va = (1+r) .* c .^ (-1/σ)
    #println(Va)
    a_old = zeros(size(coh))
    for i in 1:10_000
        Va, a, c = backward_iteration(Va,Π, a_grid, y, r, β, σ)
        if i > 1 && maximum(abs.(a - a_old)) < tol
            return Va, a, c
        end
        a_old = a
    end
end


a_grid = discretize_assets(0,10_000,500)
Va,a,c = policy_ss(Π_p3,a_grid,y_p3,r_p3,β_p3,σ_p3)

function get_lottery(a, a_grid)
    a_i = searchsortedlast(a_grid, a)
    a_π = (a_grid[a_i + 1] - a)/(a_grid[a_i + 1] - a_grid[a_i])
    if a_π == 1.0
        a_π = 0.0
    end
    return a_i, a_π
end

#get_lottery(a[6,1], a_grid)
#print(get_lottery(a[1,1], a_grid))

a_i, a_π = zeros(size(a)), zeros(size(a))
for i in 1:size(a)[1]
    for j in 1:size(a)[2]
        a_i[i,j], a_π[i,j] = get_lottery(a[i,j], a_grid)
    end
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

function distribution_ss(Π, a, a_grid, tol = 1E-6)
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
        println("iteration ", i)
        D_new = forward_iteration(D, Π, a_i, a_π)
        println(maximum(abs.(D_new - D)))
        if maximum(abs.(D_new - D)) < tol
            return D_new
        end
        D = D_new
    end
end

D = distribution_ss(Π_p3, a, a_grid)
cplot = plot((ones(2)*a_grid')',c')
cplot_zoom1 = plot((ones(2)*a_grid')'[1:121,:],c'[1:121,:])
splot_zoom1 = plot((ones(2)*a_grid')'[1:301,:], a'[1:301,:] - (ones(2)*a_grid')'[1:301,:])

function plot_cwd(a_grid,D,cap)
    Dpv = zeros(cap)
    
    for i in 1:cap
        Dpv[i] += sum(D[:,i])
    end
    
    return plot(a_grid[1:cap],cumsum(Dpv))
end

cwdplot = plot_cwd(a_grid,D,length(a_grid))
cwdplot_zoom1 = plot_cwd(a_grid,D,argmax(a_grid .> 20))
cwdplot_zoom2 = plot_cwd(a_grid,D,argmax(a_grid .> 2))