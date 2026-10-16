T_p1 = 40
β_p1 = 0.96
α_p1 = 0.7
δ_p1 = 0.04
h0_p1 = 1.0
w_p1 = 1.0

function h_evo(δ, a, α, s, h)
    return (1 - δ)*h + a * (s*h)^α
end

function x(δ, a, α, β, T, t, w)
    return (a*α*β*Q(δ, β, T, t, w))
end

function Q(δ, β, T, t, w)
    return w*((1 - (β*(1 - δ))^(T - 1 - t) )/(1 - β*(1 - δ)))
end

function s(x, h)
    return x/h
end

function earnings(s,h)
    return h*(1-s)
end

#part a

tvec = range(0.0, 41.0, 42)
svec = zeros(length(tvec))
hvec = zeros(length(tvec))
hvec[1] = h0_p1
xvec = zeros(length(tvec))
evec = zeros(length(tvec))


for i in 1:length(tvec) - 1
    xvec[length(tvec) - i] = x(δ_p1, 0.15, α_p1, β_p1, T_p1, length(tvec) - i, w_p1)
end

svec[1] = s(xvec[1],hvec[1])
evec[1] = earnings(s,h)

for i in 2:length(xvec)
    hvec[i] = h_evo(δ_p1, 0.15, α_p1, svec[i-1], hvec[i-1])
    svec[i] = s(xvec[i],hvec[i])
    evec[i] = earnings(svec[i], hvec[i])
end

bins = range(0.12,0.18,7)