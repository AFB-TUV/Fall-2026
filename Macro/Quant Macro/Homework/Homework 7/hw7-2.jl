using Distributions, LinearAlgebra, Plots

const α = 2/3
const r = 0.2
const δ = 0.1
const δ1 = 0.0
const δ2 = 0.2
const cf = 20
const ce = 40
const ce1 = 30
const ce2 = 50
const θ = 100
#const ρ = 0.9
#const ρ = 0.75
const ρ = 0.5
const p_in = 0.9
const lz = 1.39
const σ = 0.2
const λ = 1
const z_length = 5
const tol = 0.00001
const h = 0.0001

D(p) = θ/p

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

zp = Tauchen(λ,z_length,ρ,σ,0)
z = zeros(z_length)
Γ = zp[2]

for i=1:z_length
  z[i] = exp(zp[1][i]+1.39)
end
G = (1/z_length)*ones(z_length)

#m_in = (1/z_length)*ones(z_length)
#m_plus = transpose(m_in)*Γ
#m_star = transpose(m_in)*(Γ^1000)


nmin = 0.0
nmax = α*0.9*z[z_length]^(1/(1-α))
n_length = 50
n = collect(range(nmin, length=n_length, stop=nmax))

V_0 = zeros(n_length, z_length)

V(i,j,k,Vold,p,δ) = p*z[i]*n[k]^α - n[k] - δ*max(n[j]-n[k],0) - cf + (1/(1+r))*dot(Γ[i,:],max.(Vold[:,k],-δ*n[k]))
Ve(ce,V) = -ce + (1/(1+r))*dot(G,max.(V[:,1],0))

function solveBell(p::Number,δ::Number)
  local Vold = zeros(z_length, n_length)
  local Vnew = zeros(z_length, n_length)
  local Vind = zeros(z_length, n_length)
  local Vtemp = zeros(n_length)
  local diff = 0

  while true
    for i=1:z_length
      for j=1:n_length
        for k=1:n_length
          Vtemp[k] = V(i,j,k,Vold,p,δ)
        end
        Vnew[i,j] = maximum(Vtemp)
        Vind[i,j] = argmax(Vtemp)
        Vtemp = zeros(n_length)
      end
    end

    diff = maximum(Vnew-Vold)
    if diff < tol
      Vind = convert(Array{Int,2},Vind)
      break
    end

    Vind = zeros(z_length, n_length)
    Vold = zeros(z_length, n_length)
    Vold = Vnew
    Vnew = zeros(z_length, n_length)
  end

  n_pol = zeros(z_length, n_length)
  x_pol = ones(z_length, n_length)
  for i=1:z_length
    for j=1:n_length
      V_aux = zeros(z_length, n_length)
      V_aux[i,j] = max(Vnew[i,Vind[i,j]],-δ*n[Vind[i,j]])
      if abs(V_aux[i,j] - Vnew[i,j]) < tol
        x_pol[i,j] = 0
      end
      n_pol[i,j] = n[Vind[i,j]]
    end
  end

  return Vnew, Vind, n_pol, x_pol
end

#Part A

sol1 = solveBell(p_in,δ)
Val1 = sol1[1]
dom1 = sol1[2]
n1 = sol1[3]
x1 = sol1[4]


#plot_V1 = plot([Val1[1,:] Val1[3,:] Val1[5,:]], title="Value Functions")
#savefig("val")
#plot_n1 = plot([n1[1,:], n1[3,:], n1[5,:]], title="Hiring Policy Functions")
#savefig("npol")
#plot_x1 = plot([x1[1,:], x1[3,:], x1[5,:]], title="Exit Policy Functions")
#savefig("xpol")

#Part B

#ve1 = Ve(40, Val1)
#ve2 = Ve(30, Val1)
#ve3 = Ve(50, Val1)

function Newton_ve(p,δ,ce)
  local p_it = 0
  local ve_f = 0
  while true
    val = solveBell(p,δ)[1]
    valh = solveBell(p+h,δ)[1]
    p_it = p - h*Ve(ce,val)/(Ve(ce,valh) - Ve(ce,val))
    if abs(p_it - p) < tol
      break
    end
    p = 0
    p = p_it
    p_it = 0
  end
  return p_it
end

#p1 = Newton_ve(p_in,δ,ce)
#p2 = Newton_ve(p_in,δ,ce1)
#p3 = Newton_ve(p_in,δ,ce2)

#Part c

function Λ(m0::Number, δ::Number, ce::Number)
  Ω = ones(z_length,n_length)
  p_star = Newton_ve(1.2,δ,ce)
  sol = solveBell(p_star, δ)
  val = sol[1]
  domv = sol[2]
  npol = sol[3]
  xpol = sol[4]

  In = zeros(z_length, n_length)
  In_n = zeros(z_length,n_length,n_length)

  for i=1:n_length
    if n[i]==0
      In[i] = 1
    end
  end

  for i=1:z_length
    for j=1:n_length
      for k=1:n_length
        if npol[i,j] == n[k]
          In_n[i,j,k] = 1
        end
      end
    end
  end

  Λ0 = ones(z_length, n_length, 1, 1)
  Λ1 = zeros(z_length, n_length, z_length, n_length)
  Λret = zeros(z_length, n_length, 1, 1)
  while true
    for i=1:z_length
      for j=1:n_length
        for k=1:z_length
          for l=1:n_length
            Λ1[i,j,k,l] = Λ0[k,l,1,1]*(1-xpol[k,l]) * Γ[k,i]*In_n[k,l,j] + m0*G[i]*In[j]
          end
        end
      end
    end

    Λt1 = sum(Λ0, dims=4)
    Λret = sum(Λt1, dims=3)
    if abs(maximum(Λret - Λ0))<tol
      break
    else
      Λ0 = zeros(z_length, n_length, 1, 1)
      Λ0 = Λret
      Λ1 = zeros(z_length, n_length, z_length, n_length)
      Λret = zeros(z_length, n_length, 1, 1)
    end
  end
  return Λret, Λ1, npol, xpol
end

solc = Λ(0.10, δ, ce)
cca = solc[1]
ccb = solc[2]
ccc = solc[3]
ccd = solc[4]

Λplot = zeros(z_length, n_length)
for i=1:z_length
  for j=1:n_length
    Λplot[i,j] = ccb[i,j,1,1]
  end
end

#Part d
function compstats(trans::Array{Float64,4}, otrans::Array{Float64,4}, npol::Array{Float64,2}, xpol::Array{Float64,2})
Q1 = zeros(z_length, n_length)
Q2 = zeros(z_length, n_length)
for i=1:z_length
  for j=1:n_length
    Q1[i,j] =  otrans[i,j,1,1] * z[i] *npol[i,j]^α * (1 - xpol[i,j])
    Q2[i,j] = otrans[i,j,1,1] * z[i] * npol[i,j]^(α-1) * (1 - xpol[i,j])
  end
end

sQ1 = sum(Q1)
sQ2 = sum(Q2)
total = sum(otrans)
exit = sum(otrans.*xpol)/total
L = sum(otrans.*npol.*(ones(z_length, n_length)-xpol))

adjc = zeros(z_length,n_length)
for i=1:z_length
  for j=1:n_length
    adjc = otrans[i,j,1,1] * δ * max(n[j]-npol[i,j],0)*(1-xpol[i,j])
  end
end
tadjc = sum(adjc)

return sQ1, sQ2, total, exit, L, tadjc
end


#use this again for Problem 2 Part b with appropriate ρ

stats14 = compstats(cca, ccb, ccc, ccd)
#Question 2

#part a

sol2_1 = Λ(0.02, δ1, ce)
sol2_1s = compstats(sol2_1[1], sol2_1[2], sol2_1[3], sol2_1[4])

sol2_2 = Λ(0.0017, δ2, ce)
sol2_2s = compstats(sol2_2[1], sol2_2[2], sol2_2[3], sol2_2[4])
