#compute solution to McCall model using given parameters
#could be analytical but easier to do by computer

β = 0.99
σvec = [0 0.02 0.05]
b = 0.3
function coeff(β, σ)
    return 0.5*β/(1 - β*(1 - σ))
end

function sqrtarg(c1,c2,c3)
    return c2^2 - 4*c1*c3 
end

function w_star_plus(C, b)
    c1 = C
    c2 = -1*(2*C + 1)
    c3 = C + b
    return (-c2 + sqrt(sqrtarg(c1,c2,c3)))/(2*c1)
end

function w_star_minus(C,b)
    c1 = C
    c2 = -1*(2*C + 1)
    c3 = C + b
    return (-c2 - sqrt(sqrtarg(c1,c2,c3)))/(2*c1)
end

function λ(w)
    return 1 - w
end

function u(σ, λ)
    return σ/(σ + λ)
end

Cvec = zeros(size(σvec))
coeffvec = zeros(size(Cvec))
plus_wstar = zeros(size(Cvec))
minus_wstar = zeros(size(Cvec))
λvec = zeros(size(Cvec))
uvec = zeros(size(Cvec))
for i in 1:length(σvec)
    Cvec[i] += coeff(β, σvec[i])
    plus_wstar[i] += w_star_plus(Cvec[i], b)
    minus_wstar[i] += w_star_minus(Cvec[i], b)
    λvec[i] += 1 - minus_wstar[i]
    uvec[i] += σvec[i]/(σvec[i] + λvec[i])
end

function calc_b(λ, C)
    w = 1 - λ
    return -1*(C*((w - 1)^2) - w)
end

b2 = calc_b(f_bar, coeff(β, b))