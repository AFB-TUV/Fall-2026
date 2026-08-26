using StatsPlots, LaTeXStrings, Distributions, LinearAlgebra
using Roots, NLsolve, Statistics

n = 50
dist = BetaBinomial(n, 200, 100)
@show support(dist)
p = pdf.(dist,support(dist))
w = range(10.0,60.0,length=n+1)
plt1 = plot(w, p, xlabel = "wages", ylabel = "probabilities", legend = false)

wage(i) = w[i+1]
w = wage.(support(dist))
E_w = dot(p,w)
E_w_2 = dot(p,w .^2) - E_w^2
@show E_w, E_w_2

@show sum(p .* w)

c = 25
β = 0.99
num_plots = 6

T(v) = max.(w/(1-β), c + β * dot(v,p))
vs = zeros(n+1, 6)
vs[:,1] .= w/(1-β)

for col in 2:num_plots
    v_last = vs[:,col-1]
    vs[:,col] .= T(v_last)
end
plt2 = plot(vs) 

function compute_reservation_wage_direct(params;
    v_iv = collect(w ./ (1-β)),
    max_iter = 500, tol = 1e-6)
    ( ; c, β, w) = params
    T(v) = max.(w / (1-β), c + β * dot(v,p))

    v = copy(v_iv)
    v_next = similar(v)
    i = 0
    error = Inf
    while i < max_iter && error > tol
        v_next .= T(v)
        error = norm(v_next - v)
        i += 1
        v .= v_next
    end
    return (1 - β) * (c + β * dot(v,p))
end

function mccall_model(; c = 25.0, β = 0.99, w = range(10.0,60.0,length = n+1))
    (; c, β, w)
end

compute_reservation_wage_direct(mccall_model())

grid_size = 25
R = zeros(grid_size, grid_size)

c_vals = range(10.0, 30.0, length = grid_size)
β_vals = range(0.9, 0.99, length = grid_size)

for (i,c) in enumerate(c_vals)
    for (j,β) in enumerate(β_vals)
        R[i,j] = compute_reservation_wage_direct(mccall_model(; c, β))
    end
end

plt3 = contour(c_vals,β_vals,R',
    title = "Reservation Wage",
    xlabel = L"c",
    ylabel = L"\beta",
    fill = true)