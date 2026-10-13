using Distributions, LinearAlgebra, Plots

const α = 2/3
const r = 0.20
const δ = 0.1
const c_f = 20
const c_e = 40
const θ = 100
const ρ = 0.9
const lz = 1.39
const σ = 0.2

D(p) = θ/p

function V(z,n_t,n_plus,V_vec,Γ_vec)
  local del_vec = ones(5)*(-δ*n_plus)
  local exit_vec = max.(del_vec,V_vec)
  if n_t > n_plus
    ret = p*z*n_plus^α - n_plus - c_f + ((1+r)^(-1)) * (dot(Γ_vec,exit_vec))
  else
    ret = p*z*n_plus^α - n_plus - δ - c_f + ((1+r)^(-1)) * (dot(Γ_vec,exit_vec))
  end
  return ret
end



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




#Part A
zp = Tauchen(λ,20,ρ,σ,zl)
