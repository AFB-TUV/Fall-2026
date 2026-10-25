using LinearAlgebra, BasicInterpolators, Plots

function discretize_assets(amin,amax,n_a)
    ubar = log(1+log(1+amax-amin))
    ugrid = range(0,ubar,n_a)
    return amin .+ exp.(exp.(ugrid) .- 1) .- 1 
end

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
        if i % 10 == 0
            if maximum(abs.(π_new - π)) < tol
                return π_new
            end
        end
        π = π_new
    end
end

function discretize_income(ρ,σ,n_e)
    p = (1+ρ)/2

    e = range(0,n_e-1,n_e)
    α = 2*σ/sqrt(n_e-1)
    e = α .* e

    Π = rouwenhorst_Π(n_e, p)
    π = stationary_markov(Π)

    y = exp.(e)
    y = y./dot(π,y)
    return y,π,Π
end

#y,π,Π = discretize_income(0.975,0.7,7)

#m = rouwenhorst_Π(7,1.975/2)
#d = dot(π,y)
#meanlogy = dot(π,log.(y))
#sdlogy = sqrt(dot(π,(log.(y) .- meanlogy).^2))
#sdlogy

function backward_iteration(Va, Π, a_grid, y, r, β, σ)
    Wa = (β .* Π) * Va

    c_endog = Wa .^ -σ
    coh = broadcast(+,y,(1+r).*a_grid')

    a = zeros(size(coh))
    for e in 1:size(y)[1]
        #println(c_endog[e,:])
        itp_c = LinearInterpolator(c_endog[e,:] + a_grid, a_grid,NoBoundaries())
        a[e,:] = itp_c.(coh[e,:])
        #a[e, :] = np.interp(coh[e, :], c_endog[e, :] + a_grid, a_grid)
    end

    a = max.(a,a_grid[1])
    c = coh - a

    Va = (1+r) .* c .^ (-1/σ)
    return Va, a, c
end

function policy_ss(Π, a_grid, y, r, β, σ, tol = 1E-9)
    coh = y .+ (1.0 + r).*a_grid'
    c = 0.05 .* coh
    Va = (1+r) .* c .^ (-1/σ)
    #println(Va)
    a_old = zeros(size(coh))
    for i in 1:10_000
        Va, a, c = backward_iteration(Va,Π, a_grid, y, r, β, σ)
        if i % 10 == 0
            if maximum(abs.(a - a_old)) < tol
                return Va, a, c
            end
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

#=
#get_lottery(a[6,1], a_grid)
#print(get_lottery(a[1,1], a_grid))

#a_i, a_π = zeros(size(a)), zeros(size(a))

for i in 1:size(a)[1]
    for j in 1:size(a)[2]
        a_i[i,j], a_π[i,j] = get_lottery(a[i,j], a_grid)
    end
end
#println(a_π[1,:])
=#

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

function distribution_ss(Π, a, a_grid, tol = 1E-10)
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
        if i % 10 == 0
            if maximum(abs.(D_new - D)) < tol
                return D_new
            end
        end
        D = D_new
    end
end


function steady_state(model)
    (; Π, a_grid, y, r, β, σ) = model
    Va, a, c = policy_ss(Π, a_grid, y, r, β, σ)
    D = distribution_ss(Π, a, a_grid)
    A = dot(a,D)
    C = dot(c,D)
    return (; Va = Va, a = a, c = c, D = D, a_grid = a_grid, y = y, Π = Π, A = A, C = C)
end

function Model(; a_min = 0.0, a_max = 10_000.0, n_a = 500, r = 0.01/4, β = 1-0.08/4, σ = 1.0, ρ = 0.975, σy = 0.7, N_y = 7)
    a_grid = discretize_assets(a_min,a_max,n_a)
    y,π,Π = discretize_income(ρ, σy, N_y)
    return (; a_min, a_max, n_a, a_grid, r, β, σ, ρ, σy, N_y, y, π, Π)
end

model = Model()
Va,a,c = policy_ss(model.Π,model.a_grid,model.y,model.r,model.β,model.σ)
D = distribution_ss(model.Π, a, model.a_grid)
ss = steady_state(model)

cplot = plot((ones(7)*model.a_grid')',ss.c')
cplot_zoom1 = plot((ones(7)*model.a_grid')'[1:121,:],ss.c'[1:121,:])
cplot_zoom2 = plot((ones(7)*model.a_grid')'[1:121,1:4],ss.c'[1:121,1:4])
splot_zoom1 = plot((ones(7)*model.a_grid')'[1:301,:], ss.a'[1:301,:] - (ones(7)*model.a_grid')'[1:301,:])


function plot_cwd(a_grid,D,cap)
    Dpv = zeros(cap)
    
    for i in 1:cap
        Dpv[i] += sum(D[:,i])
    end
    
    return plot(a_grid[1:cap],cumsum(Dpv))
end

cwdplot = plot_cwd(model.a_grid,ss.D,length(model.a_grid))
cwdplot_zoom1 = plot_cwd(model.a_grid,ss.D,argmax(model.a_grid .> 20))
cwdplot_zoom2 = plot_cwd(model.a_grid,ss.D,argmax(model.a_grid .> 2))

println("Conducting comparative statics on interest rates...")

rs = range(-0.02,0.015,15) .+ model.r
As = zeros(size(rs))

for i in 1:length(rs)
    atemp = policy_ss(model.Π, model.a_grid, model.y, rs[i], model.β, model.σ)[2]
    Dtemp = distribution_ss(model.Π, atemp, model.a_grid)
    As[i] = dot(atemp,Dtemp)
end

rcsplot = plot(rs,As)

println("Conducting comparative statics on income volatility...")
    

sds = range(0.3,1.2,8)
As2 = zeros(size(sds))
for i in 1:length(sds)
    y_new,π_new,Π_new = discretize_income(0.975,sds[i],7)
    atemp = policy_ss(Π_new,model.a_grid,y_new,model.r,model.β,model.σ)[2]
    Dtemp = distribution_ss(Π_new,atemp,model.a_grid)
    As2[i] = dot(atemp,Dtemp)
end

sdcsplot = plot(sds,As2)

println("Conducting comparative statics on EIS...")
    

σs = range(0.4,2,10)
As3 = zeros(size(σs))
for i in 1:length(σs)
    atemp = policy_ss(model.Π, model.a_grid, model.y, model.r, model.β, σs[i])[2]
    Dtemp = distribution_ss(model.Π, atemp, model.a_grid)
    As3[i] = dot(atemp, Dtemp)
end

σcsplot = plot(σs,As3)

function calibrate_β(M; βa = 0.98, βb = 0.99, tol = 1e-12, maxits = 1000)
    Ma = steady_state(M(β = βa))
    Mb = steady_state(M(β = βb))
    fa, fb = dot(Ma.a, Ma.D) - 5.6, dot(Mb.a, Mb.D) - 5.6
    local βc
    for i in 1:maxits
        βc = (βa+βb)/2
        Mc = steady_state(Model(β = βc))
        fc = dot(Mc.a, Mc.D) - 5.6
        if abs(βb - βa) < tol
            return M(β = βc), steady_state(M(β = βc))
        elseif fa*fc > 0
            βa = βc
            Ma = Mc
            fa = fc
        else
            βb = βc
            Mb = Mc
            fb = fc
        end
    end
    println("Maximum iterations exceeded!")
    return 0
end

#cM, cS = calibrate_β(Model)

function GEModel(; a_min = 0.0, a_max = 10_000.0, n_a = 500, r = 0.01/4, β = 1-0.08/4, σ = 1.0, ρ = 0.975, σy = 0.7, N_y = 7, B = 5.6)
    a_grid = discretize_assets(a_min,a_max,n_a)
    y,π,Π = discretize_income(ρ, σy, N_y)
    τ = r*B
    y = (1-τ)*y    
    return (; a_min, a_max, n_a, a_grid, r, β, σ, ρ, σy, N_y, y, π, Π, τ, B)
end

function calibrate_β_GE(M; βa = 0.98, βb = 0.99, tol = 1e-12, maxits = 1000)
    B = M().B
    Ma = steady_state(M(β = βa))
    Mb = steady_state(M(β = βb))
    fa, fb = dot(Ma.a, Ma.D) - B, dot(Mb.a, Mb.D) - B
    local βc
    for i in 1:maxits
        βc = (βa+βb)/2
        Mc = steady_state(Model(β = βc))
        fc = dot(Mc.a, Mc.D) - B
        if abs(βb - βa) < tol
            return M(β = βc), steady_state(M(β = βc))
        elseif fa*fc > 0
            βa = βc
            Ma = Mc
            fa = fc
        else
            βb = βc
            Mb = Mc
            fb = fc
        end
    end
    println("Maximum iterations exceeded!")
    return 0
end

println("Calibrating general equilibrium...")

gemodel = GEModel()
gM, gS = calibrate_β_GE(GEModel)
clear = (gS.A - gM.B, gS.C - 1)

MPCs = zeros(size(gS.c))
MPCs[:, 2:end-1] .= (gS.c[:, 3:end] .- gS.c[:, 1:end-2]) ./ (gM.a_grid[3:end] .- gM.a_grid[1:end-2]) ./ (1 + gM.r)