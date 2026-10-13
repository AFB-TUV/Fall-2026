using Plots
#fundamental parameters
const β = 0.96
const θ = 1/1.7
const α = 0.29
const γ = 2
const ω = 2
const A = 10
const σ = 0.051
const λ = 0.44
const p = 0.5*(1+λ)

#constants for algorithm
const tol = 0.001
#const h = 0.00000001

#define the capital grid
cap_grid = collect(range(250, length=100, stop=1000))

#define functions
#Labor FOC
l(k,ϵ) = (1-α)^((ω + α)/((θ+α) * (ω + α) + α - 1)) * (α)^((ω + α)/((θ+α) * (ω+α) + α - 1)) * A^((ω + 2*α)/((θ+α) * (ω+α) + α - 1)) * k^( α * (ω + 2*α)/((θ+α) * (ω+α) + α - 1)) * exp((α*ϵ)/((θ+α) * (ω+α) + α - 1))
#Capital utilization FOC
h(k,ϵ) = (α*A)^(1/(ω+α)) * (k/l(k,ϵ))^((α-1)/(ω+α)) * exp(ϵ/(ω+α))
#Production Function
f(k,ϵ) = A*((k*h(k,ϵ))^α)*(l(k,ϵ)^(1-α))
#Consumption
c(k0,k1,ϵ0,ϵ1) = f(k0,ϵ0) - exp(ϵ1)*k1 + k0*(1-(h(k0,ϵ0)^ω)/ω)*exp(ϵ0)
#Bellman
V(k0,k1,ϵ0,ϵ1,V1) = (1/(1-γ))*(c(k0,k1,ϵ0,ϵ1) - (l(k0,ϵ0)/(1+θ))^(1+θ))^(1-γ) + β * V1
ExV(k0,k1,ϵ0,V1) = p*V(k0,k1,ϵ0,ϵ0,V1) + (1-p)*V(k0,k1,ϵ0,-ϵ0,V1)

#Maximize the Bellman equation for a given V' and k_t.
#Return V and k_{t+1}
function maxBell(v0::Number, shock::Number, k::Number, kgrid::Array{Float64,1})
    local n = length(kgrid)
    local vnew = zeros(n)

    for i=1:n
        if shock==1
            if c(k,kgrid[i],σ,σ) < 0 || c(k,kgrid[i],σ,-σ) < 0
                vnew[i] = -1000000
            else
                #vnew[i] = p*V(k,kgrid[i],σ,σ,v0) + (1-p)*V(k,kgrid[i],σ,-σ,v0)
                vnew[i] = ExV(k,kgrid[i],σ,v0)
            end
        else
            if c(k,kgrid[i],-σ,σ) < 0 || c(k,kgrid[i],-σ,-σ) < 0
                vnew[i] = -1000000
            else
                #vnew[i] = p*V(k,kgrid[i],-σ,-σ,v0) + (1-p)*V(k,kgrid[i],-σ,σ,v0)
                vnew[i] = ExV(k,kgrid[i],-σ,v0)
            end
        end
    end
    local sol = findmax(vnew)
    return sol[1], kgrid[sol[2]]
end

#Policy function iteration; no improvements
function solveBell(kgrid::Array{Float64,1})
    local n = length(kgrid)
    local vold = zeros(n,2)
    local vnew
    local savings
    local z=0

    while true
        z=z+1
        vnew = zeros(n,2)
        savings = zeros(n,2)
        for i=1:n
            for j=1:2
                local sol = maxBell(vold[i,j],j,kgrid[i],kgrid)
                vnew[i,j] = sol[1]
                savings[i,j] = sol[2]
            end
        end

        local vcheck = vnew - vold
        if abs(maximum(vcheck)) < tol
            break
        else
            vold = zeros(n,2)
            vold = vnew
        end
    end
    return vnew, savings, z
end

sol = solveBell(cap_grid)
v = sol[1]
k = sol[2]
v1 = v[:,1]
v2 = v[:,2]
k1 = k[:,1]
k2 = k[:,2]
iterations = sol[3]
#@time solveBell(cap_grid)

plot_v = plot(cap_grid,[v1 v2], title="Value Function")
plot_k = plot(cap_grid,[k1 k2], title="Savings")
