using LinearAlgebra, BasicInterpolators

function discretize_assets(amin,amax,n_a)
    ubar = log(1+log(1+amax-amin))
    ugrid = range(0,ubar,n_a)
    return amin .+ exp.(exp.(ugrid) .- 1) .- 1 
end

#agrid = discretize_assets(0,10_000,50)

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
    println(sum(Va))
    return Va, a, c
end

function policy_ss(Π, a_grid, y, r, β, σ, tol = 1E-9)
    coh = broadcast(+,y,(1+r).*a_grid')
    c = 0.05 .* coh
    Va = (1+r) .* c .^ (-1/σ)
    #println(Va)
    a_old = zeros(size(coh))
    for i in 1:10_000
        Va, a, c = backward_iteration(Va,Π, a_grid, y, r, β, σ)
        #if i < 200
        #    println("iteration ", i, ":   ", sum(a))
        #end
        if i > 1 && maximum(abs.(a - a_old)) < tol
            return Va, a, c
        end
        a_old = a
    end
end


a_grid = discretize_assets(0,10_000,500)
y,π,Π = discretize_income(0.975,0.7,7)
r = 0.01/4
β = 1-0.08/4
σ = 1
Va,a,c = policy_ss(Π,a_grid,y,r,β,σ)
println(sum(Va))
