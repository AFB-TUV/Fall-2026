using Plots, LinearAlgebra
#fundamental parameters
const β = 0.93
const ρ = 0.975
const σ = 2
const σa = 10
const R = 1.025
#Income states
ygrid = [16000.0; 80.0; 15.0]

#Transition Matrix
Γ = [ 0.02 0.98 0; 0.001 0.78 0.219; 0 0.036 0.964 ]

#constants for algorithm
const tol = 0.00001
#const h = 0.00000001

#functions

#consumption
c(a0,a1,y) = R*a0 + y - a1
#utility
u(a0,a1,y) = (c(a0,a1,y)^(1-σ))/(1-σ)
#Bellman
V(a0,a1,y,V1) = u(a0,a1,y) + ρ*β*V1

#solution methods

function elemax(a::Array{Float64,1}, b::Array{Float64,1})
    n = length(a)
    c = zeros(n)
    for i=1:n
        c[i] = max(a[i],b[i])
    end
    return c
end

function elesum(a::Array{Float64,2}, j::Number)
    n = size(a)[1]
    vec = zeros(n)
    for i=1:n
        vec[i] = sum(a[:,j])
    end
    return vec
end


function solvefluc(aold::Array{Float64,2}, y::Array{Float64,1}, trans::Array{Float64,2}, agrid::Array{Float64,1})

    local n = length(agrid)
    local m = length(y)
    local anew
    local its = 0

    while true

        anew = zeros(n,m)
        its = its + 1

        for j=1:m
            for i=1:n
                c0 = R*agrid[i] + y[j] - agrid[1]
                c1 = R*agrid[1]*ones(m) + y - aold[1,:]
                Euler = c0^(-σ) - ρ*β*R* dot(elesum(trans,j), c1.^(-σ))

                if Euler > 0
                    vnew[i,j] = agrid[1]
                else
                    local δ = 1
                    local Δ = zeros(n)
                    local k=0

                    while (δ > 0 && k < n)
                        k = k+1
                        c0 = 0
                        c1 = zeros(m)
                        c0 = R*agrid[i] + y[j] - agrid[k]
                        c1 = R*agrid[k]*ones(m) + y - aold[k,:]
                        δ = 0
                        δ = max(0,c0^(-σ)) - ρ*β*R*dot(elesum(trans,j), elemax(zeros(m),c1.^(-σ)))
                        Δ[k] = δ
                    end

                    if k < 2
                        anew[i,j] = agrid[k]
                    elseif abs(Δ[k]) > abs(Δ[k-1])
                        anew[i,j] = agrid[k-1]
                    else
                        anew[i,j] = agrid[k]
                    end

                end
                #ends of for loops:
            end
        end
        #check break loop condition:
        if (maximum(anew - aold) < tol || minimum(anew-aold) < tol)
            break
        else
            aold = zeros(n,m)
            aold = anew
        end
    end

    return anew, its
end

#=
function maxBell(v0::Number, a::Number, agrid::Array{Float64,1})
    local n = length(agrid)
    local vnew = zeros(n)
    local eat = zeros(n)
    for i=1:n
        eat[i] = c(a,agrid[i],y2)
        if eat[i] > 0
            vnew[i] = V(a,agrid[i],y2,v0)
        else
            vnew[i] = -100000 - 1000*abs(eat[i])
        end
    end
    local sol = findmax(vnew)
    return sol[1], agrid[sol[2]], eat[sol[2]]
end

function solveBell(agrid::Array{Float64,1})
    local n = length(agrid)
    local vold = zeros(n)
    local vnew
    local eat
    local savings
    local z=0

    while true
        z=z+1
        vnew = zeros(n)
        savings = zeros(n)
        eat = zeros(n)
        for i=1:n
            local sol = maxBell(vold[i],agrid[i],agrid)
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
=#

#main code

#asset grid
const a_l = 0
const a_la = -30
const a_max = 3*ygrid[1]

assets = collect(range(a_l, length=30, stop=a_max))
ex_assets = collect(range(a_l, length=500, stop=a_max))
assets_a = collect(range(a_la, length=30, stop=a_max))
ex_assets_a = collect(range(a_la, length=500, stop=a_max))

a_guess = zeros(30,3)
a_guess_a = zeros(500,3)

sol1 = solvefluc(a_guess, ygrid, Γ, assets)
asol1 = sol1[1]
its1 = sol1[2]
@time solvefluc(a_guess,ygrid,Γ,assets)

plot_a = plot(assets,asol1, title="Savings Policy Rule")
