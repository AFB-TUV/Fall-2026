#problem 1

b = 1.0
r_1 = 0.04
g_1 = 0.02
r_2 = 0.015
g_2 = 0.02

function surplus(b, r, g)
    return ((1+r)/(1+g))*b - b
end

s1 = surplus(b, r_1, g_1)
s2 = surplus(b, r_2, g_2)

#problem 2

α = 1.54
θ = 0.3224
cy = 0.75
τ1 = 0.4
τ2 = 0.6

function h(α, θ, cy, τ)
    return ((1 - θ)*(1 - τ))/((1-θ)*(1-τ)+(α*cy))
end

h1 = h(α, θ, cy, τ1)
h2 = h(α, θ, cy, τ2)