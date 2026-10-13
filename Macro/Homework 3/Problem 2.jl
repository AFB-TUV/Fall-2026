β_p2 = 0.94
r_p2 = 0.04
w_p2 = 1.0
ν_p2 = 0.5
ρ_p2 = 0.95
σ_p2 = 0.4

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

function VFI(V,  ; maxits = 10_000, tol=1E-6)
    Vnew = ones(size(V))

    for i in 1:maxits
        Vnew, cpol, hpol, apol = backit(V)
        if i > 1 && maximum(abs.(Vnew - V)) < tol
            return Vnew, cpol, hpol, apol
        end
        V = Vnew
    end
    println("Value function failed to converge!")
end

function backit(Vprime, )
    Wprime = (β .* Π) * Vprime
    V = 


end