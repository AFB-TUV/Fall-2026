capcost_0 = 75000
subs_0 = 200
price_0 = 350
fixcost_0 = 25000
varcost_0 = 80
interest_0 = 0.15

cf_0 = (subs_0 * price_0) - (fixcost_0 + varcost_0 * subs_0)

function npv_rates(i)
    return (1+i)^(-1) + (1+i)^(-2) + (1+i)^(-3) + (1+i)^(-4) + (1+i)^(-5)
end

function npv(i, cf, cc)
    return npv_rates(i) * cf - cc
end

function breakeven_capcost(a, b, cf, cc, tol = 1e-5, maxits = 1000)
    fa, fb = npv(a, cf, cc), npv(b, cf, cc)
    i = 0
    local fc
    for i in 1:maxits
        c = (a+b)/2
        fc = npv(c, cf, cc)
        if b - a < tol
            return c
        elseif fa*fc > 0
            a = c
            fa = fc
        else
            b = c
            fb = fc
        end
    end
    println("Maximum iterations exceeded!")
    return 0
end

i_be = breakeven_capcost(0.0, 1.0, cf_0, capcost_0)