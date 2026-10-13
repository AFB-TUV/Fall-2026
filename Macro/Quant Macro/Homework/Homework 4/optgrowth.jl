using Plots
#fundamental parameters
const β = 0.95
const δ = 0.1
const α = 0.33
const σ = 2

#constants for algorithm
const tol = 0.0001
const h = 0.00000001
#const h = 0.00000001

#define the capital grid
cap_grid = collect(range(-20, length=100, stop=20))

#define functions
#Production Function
f(k) = (k^α)
#Consumption
c(k0,k1) = f(k0) + (1-δ)*k0 - k1
#Bellman
V(k0,k1,V1) = (1-σ)*(c(k0,k1) - 1)^(1-σ) + β * V1

function Newton(f::Function, k_0::Number, k_1::Number, V0::Number, h::Number, tol::AbstractFloat)
    #initialize iteration tracker and next capital
    i=0
    local k_1
    #run while loop until broken inside
    while true
        #Do Newton's method
        V1 = V(k_0,k_1,V0)
        #Increase iterations
        i = i + 1
        # break if convergence is found
        if abs(V1 - V0 < tol)
            break
        end
        #if no convergence if found, replace the initial capital and start over
        V0 = V1
    end
    #return the new capital and number of iterations
    return V1
end

#Maximize the Bellman equation for a given V' and k_t.
#Return V and k_{t+1}
function maxBell(v0::Number, k::Number, kgrid::Array{Float64,1})
    local n = length(kgrid)
    local vnew = zeros(n)
    local eat = zeros(n)
    for i=1:n
        eat[i] = c(k,kgrid[i])
        if eat[i] > 0
            vnew[i] = Newton(k,kgrid[i],v0)
        else
            vnew[i] = -Inf
        end
    end

    sol_help =

    local sol = findmax(vnew)
    return sol[1], kgrid[sol[2]], eat[sol[2]]
end

#Policy function iteration; no improvements
function solveBell(kgrid::Array{Float64,1})
    local n = length(kgrid)
    local vold = zeros(n)
    local vnew
    local savings
    local eat
    local z=0

    while true
        z=z+1
        vnew = zeros(n)
        savings = zeros(n)
        eat = zeros(n)
        for i=1:n
            local sol = maxBell(vold[i],kgrid[i],kgrid)
            vnew[i] = sol[1]
            savings[i] = sol[2]
            eat[i] = sol[3]
        end

        local vcheck = vnew - vold
        if abs(maximum(vcheck)) < tol
            break
        else
            vold = zeros(n)
            vold = vnew
        end
    end

    return vnew, savings, eat, z
end

@time sol = solveBell(cap_grid)
v = sol[1]
k = sol[2]
eat = sol[3]
iterations = sol[4]
#@time solveBell(cap_grid)

plot_v = plot(cap_grid,v, title="Value Function")
plot_k = plot(cap_grid,k, title="Savings")
plot_c = plot(cap_grid,eat, title="Consumption")
