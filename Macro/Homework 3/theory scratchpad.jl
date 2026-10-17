#problem 2

function TP2(; B = 0.3, hbar = 1/3, w = 1.0, r = 0.04)
    return (; B, hbar, w, r)
end

function e_star(model, a)
    (; B, r, w, hbar) = model
    return ((1 - exp(-B)* (1+r))/(exp(-B)*w*hbar))*a
end

function int_e_star(model, a)
    (; B, r, w, hbar) = model
    return 0.5 * a * e_star(model, a)
end

function ϵN(model, a1, a2)
    N = 0.5* max(1 - e_star(model, a1), 0.0) + 0.5*max(1 - e_star(model, a2), 0.0)
    ϵ = (0.5*int_e_star(model, a1) + 0.5*int_e_star(model, a2))/N
    return ϵ, N
end

tp2 = TP2()

#part a

κ = e_star(tp2, 1.0)

#part c

ϵi, Ni = ϵN(tp2, 0.5, 0.5)
ϵii, Nii = ϵN(tp2, 0.05, 0.95)
ϵiii, Niii = ϵN(tp2, 0.0, 1.0)


#part d

ϵiv, Niv = ϵN(tp2, 0.1, 0.9)

#problem 5

function bargaining(model, a, w)
    (; b, γ, z) = model
    return (1-γ)*(log(a+w) - log(a+b)) - γ * (z-w)/(a+w)
end

function bisection(model, a; w_a = 1e-2, w_b = 10.0, tol = 1e-8, maxits = 10_000)
    fa, fb = bargaining(model, a, w_a), bargaining(model, a, w_b)
    i = 0
    local w_c
    for i in 1:maxits
        w_c = (w_a + w_b)/2
        fc = bargaining(model, a, w_c)
        if abs(w_b - w_a) < tol
            return w_c
        elseif fa*fc > 0
            w_a = w_c
            fa = fc
        else
            w_b = w_c
            fb = fc
        end
    end
    println("Maximum iterations exceeded!")
    return 0
end

function RN(model)
    (; b, γ, z) = model
    return γ * z + (1-γ) * b
end

function TP5(; z = 1, b = 0.4, γ = 0.5)
    return (; z, b, γ)
end

tp5 = TP5()
tp5d = TP5(b = 0.5)
w_RN = RN(tp5)
w_1 = bisection(tp5, 0.0)
w_2 = bisection(tp5, 10.0)

w_RNd = RN(tp5d)
w_1d = bisection(tp5d, 0.0)
w_2d = bisection(tp5d, 10.0)