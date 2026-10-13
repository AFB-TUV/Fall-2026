using Distributions, LinearAlgebra, Plots

#fundamental parameters
const ρ1 = 0.95
const σ1 = 0.20
const μ1 = 0

#grid parameters
λ1 = 2
m1 = 10
m2 = 50

function Tauchen(λ::Number, m::Number, ρ::Number, σ::Number, μ::Number)
    #discretize AR(1) process
    local z1 = μ - λ * σ * (1-ρ^2)^(-0.5)
    local step = -2*z1*(m-1)^(-1)
    local zm = z1 + (m-1)*step
    local AR1 = collect(range(z1, length=m, stop=zm))

    #Compute transition matrix
    local TransMat = zeros(m,m)
    local quant

    for j=1:m
        for k=1:m

            quant=0
            if k==1
                quant = cdf(Normal(), ((AR1[k] + 0.5*step - ρ*AR1[j])/σ) )
            elseif k==m
                quant = 1 - cdf(Normal(), ((AR1[k] - 0.5*step - ρ*AR1[j])/σ) )
            else
                quant = cdf(Normal(), ((AR1[k] + 0.5*step - ρ*AR1[j])/σ) ) - cdf(Normal(), ((AR1[k] - 0.5*step - ρ*AR1[j])/σ) )
            end
            TransMat[j,k] = quant
        end
    end

    return AR1, TransMat
end

function staDist(M::Array{Float64,2})
    local n = size(M)[1]
    local Vold = transpose(collect(range(1/n, length=n, stop=1/n)))
    local Vnew
    for i=1:10000
        Vnew = zeros(n)
        Vnew = Vold * M
        Vold = zeros(n)
        Vold = Vnew
    end
    #create return value of type Array and not type Transpose
    local Vret = zeros(n)
    for i=1:n
        Vret[i] = Vnew[i]
    end
    return Vret
end

function trueIncome(v::Array{Float64,1})
    local temp = zeros(length(v))
    for i=1:length(v)
        temp[i] = exp(v[i])
    end
    return temp
end

function tmat(A::Array{Float64,2})
    local B = zeros(size(A'))
    for i=1:size(B)[1]
        for j=1:size(B)[2]
            B[i,j] = A[j,i]
        end
    end
    return B
end

function compStats(z::Array{Float64,1}, p::Array{Float64,1})
    local y = trueIncome(z)
    local z2 = z.^2
    local y2 = y.^2

    #Variance
    local σz2 = dot(p,z2) - dot(p,z)^2
    #Coefficient of Variation
    local σy2 = dot(p,y2) - dot(p,y)^2
    local coeff = sqrt(σy2)/dot(p,y)

    #percentiles
    local medfinder = findfirst(x -> x > 0.5, cumsum(p))
    local tenfinder = findfirst(x -> x > 0.1, cumsum(p))
    local ninfinder = findfirst(x -> x > 0.9, cumsum(p))
    local med = y[medfinder]
    local ten = y[tenfinder]
    local nin = y[ninfinder]

    #mean/median
    local mm = dot(p,y)/med
    #90/50
    local nm = nin/med
    #50/10
    local mt = med/ten

    #Gini coefficient
    local g = zeros(length(y),length(y))
    for i=1:length(y)
        for j=1:length(y)
            g[i,j] = abs(y[i] - y[j]) *p[i] * p[j]
        end
    end
    local Gini = (2*dot(y,p))^(-1)*sum(g)

    #Prepare simulated Lorenz curve data
    local F = [0; cumsum(p)]
    local Ltemp = [0; cumsum(y.*p)]
    local L = Ltemp/Ltemp[length(Ltemp)]
    local lp = collect(range(0,length=length(L),stop=1))

    return σz2, coeff, Gini, mm, nm, mt, F, L, lp

end

function compDists(n::Number, λ::Number, m::Number, ρ::Number, σ::Number, μ::Number)
    σσ = collect(range(σ, length=n, step=0.005))
    #local zs = Array{Array{Float64,1}, n}
    #local ps = Array{Array{Float64,2}, n}

    #local zs = fill(zeros(n), 1, n)
    local z = Tauchen(λ, m, ρ, σσ[1], μ)[1]
    local ps = fill(zeros(n,n), 1, n)
    for i=1:n
        #zs[i] = Tauchen(λ, m, ρ, σσ[i], μ)[1]
        ps[i] = Tauchen(λ, m, ρ, σσ[i], μ)[2]
    end

    #local dists = Array{Array{Float64,1}, n}
    local dist = staDist(ps[1])
    #for i=1:n
    #   dists[i] = staDist(ps[i])
    #end

    return z, ps, dist
end

function multMats(ps = Array{Array{Float64,2}}, z = Array{Float64,1})
    n = length(z)
    m = length(ps)
    local zold = z
    local znew
    #local zs = fill(zeros(n), 1, n)
    local zlist = fill(zeros(n), 1, m)
    zlist[1] = z

    #transpose matrices
    local plist = fill(zeros(n,n), 1, n)
    for i=1:m
        for j=1:n
            for k=1:n
                plist[i][j,k] = ps[i][k,j]
            end
        end
    end

    #multiply once by every transition matrix after the first (1989)
    for i=2:m
        znew = zeros(n)
        znew = plist[i] * zold
        zlist[i] = znew
        zold = zeros(n)
        zold = znew
    end

    #create return values of type Array and not type Transpose
    #local Vret = zeros(n)
    #for i=1:n
    #    Vret[i] = Vnew[i]
    #end
    return zlist
end

#Problem 1
IP1 = Tauchen(λ1, m1, ρ1, σ1, μ1)
z1 = IP1[1]
y1 = trueIncome(z1)
M1 = IP1[2]

V1 = collect(range(0.1, length=10, stop=0.1))
V2 = V1'*IP1[2]

#Part C

#compute matrix exponential
matinf1 = M1^1000
#get first row (any will do; all are equivalent)
stad1 = matinf1[1,:]
#check against another way of finding π
Vf1 = staDist(M1)

#Part D
stats1 = compStats(z1,stad1)

#Part E
IP2 = Tauchen(λ1, m2, ρ1, σ1, μ1)
z2 = IP2[1]
y2 = trueIncome(z2)
M2 = IP2[2]
Vf2 = staDist(IP2[2])
matinf2 = M2^1000
#get first row (any will do; all are equivalent)
stad2 = matinf2[1,:]
stats2 = compStats(z2,stad2)

lorenzes = plot(stats1[7], stats1[8])
plot!(lorenzes, stats1[9], stats1[9])
savefig("p1_1")
plot!(lorenzes, stats2[7], stats2[8])
savefig("p1_2")

#Problem 3

const ncomps = 2016 - 1989 + 1
seqDists = compDists(ncomps, λ1, m2, ρ1, σ1, μ1)

grid = trueIncome(seqDists[1])
mats = seqDists[2]
dist = seqDists[3]

evp = multMats(mats, dist)

#=
    #Gini coefficient
    local g = zeros(length(y),length(y))
    for i=1:length(y)
        for j=1:length(y)
            g[i,j] = abs(y[i] - y[j]) *p[i] * p[j]
        end
    end
    local Gini = (2*dot(y,p))^(-1)*sum(g)
=#

meanIncs = zeros(ncomps)
medIncs = zeros(ncomps)
cumIncs = fill(zeros(length(grid)), 1, ncomps)
Ginis = zeros(ncomps)

for i=1:ncomps
    meanIncs[i] = dot(evp[i], grid)
    cumIncs[i] = cumsum(evp[i])
    medIncs[i] = grid[findfirst(x -> x > 0.5, cumsum(evp[i]))]
    local g = zeros(length(grid), length(grid))
    for j=1:length(grid)
        for k=1:length(grid)
            g[j,k] = abs(grid[j] - grid[k]) * evp[i][j] * evp[i][k]
        end
    end
    Ginis[i] = (2*dot(grid,evp[i]))^(-1) * sum(g)
end

times = collect(range(1989, step=1, stop=2016))
plot([times meanIncs])
savefig("p2_1")
plot([times Ginis])
savefig("p2_2")
