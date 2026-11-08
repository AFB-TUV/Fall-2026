#begin
	
using Parameters 
using Printf
using Distributions
using Roots
using LinearAlgebra
using Plots # plots: comment out if you wont use it, since it is heavy.
#using PlutoUI    # printing stuff in the Pluto notebook

function tauchen(N, ρ, σ; μ = 0.0, m = 3.0)
	s1    = μ/(1-ρ) - m*sqrt(σ^2/(1-ρ^2))
   	sN    = μ/(1-ρ) + m*sqrt(σ^2/(1-ρ^2))
    s = collect(range(s1, sN, length = N))
    step    = (s[N]-s[1])/(N-1)  #evenly spaced grid
    P      = fill(0.0, N, N)

    for i = 1:ceil(Int, N/2)
    	P[i, 1] = cdf.(Normal(), (s[1] - μ - ρ*s[i] + step/2.0)/σ)
        P[i, N]  = 1 - cdf.(Normal(), (s[N] - μ - ρ*s[i]  - step/2.0)/σ)
        for j = 2:N-1
        	P[i,j]  = cdf.(Normal(), (s[j] - μ - ρ*s[i]  + step/2.0)/σ) -
                            cdf.(Normal(), (s[j] - μ - ρ*s[i] - step/2.0)/σ)
        end
        P[floor(Int, (N-1)/2+2):end, :]=P[ceil(Int ,(N-1)/2):-1:1, end:-1:1]
	end

    ps = sum(P, dims = 2)
    P = P./ps

    return s, P
end
	
	
	
	
#end


function set_par(;
    beta = 0.8,
    rho = 0.9,
    sigma = 0.2,
    phi_mean = 1.0,
    alpha = 2/3,
    c_e = 40.0,
    c_f = 20.0,
    D_bar = 100.0,
    mu_g  = 0.0,
    sigma_g = 0.2,  
    nPhi = 101,
    w = 1.0)
    

    # === SHOCK DISCRETIZATION
    phi_mean = phi_mean*(1-rho) 
    gPhi, F_trans = tauchen(nPhi, rho, sigma; μ = phi_mean, m = 4.0)
	gPhi = @. exp(gPhi)	
	invPhi = F_trans^1000 # invariant distribution
	invPhi = invPhi[1,:]	
	
	# === ENTRANTS DISTRIBUTION: Assume they draw from the invariant distribution.
	G_prob = invPhi
	
	# Alternative: assume they draw from a log-normal
    #p_cdf = cdf.(LogNormal(mu_g, sigma_g), gPhi)
  
    #G_prob = fill(0.0, nPhi)
    #G_prob[1] = (p_cdf[2]-p_cdf[1])/2.0
    #for i = 2:nPhi-1
    #    G_prob[i] = (p_cdf[i+1]-p_cdf[i])/2.0 + (p_cdf[i]-p_cdf[i-1])/2.0
    #end
    #G_prob[nPhi] = (1.0 - p_cdf[nPhi])/2.0 + (p_cdf[nPhi]-p_cdf[nPhi-1])/2.0
    #G_prob[1] = G_prob[1] +  (1.0 - p_cdf[nPhi])/2.0 #  assign the extra mass back in the first grid point

    return (beta = beta, F_trans = F_trans, gPhi = gPhi, alpha = alpha, c_e = c_e, 
            c_f = c_f, D_bar = D_bar, nPhi = nPhi, w=w, G_prob = G_prob)
end

param = set_par()

function solve_bellman(p_guess, param)
    @unpack beta, F_trans, gPhi, alpha, c_f, nPhi, w = param
    
    # static decision:
    gN = @. (p_guess*alpha*gPhi/w)^(1/(1.0-alpha))
    gPi = @. p_guess*gPhi*gN^alpha - w*gN - c_f*w
    
    # solve bellman
    tol = 10.0^-9
    max_iter = 500
    iter_count = 10
    print_it = false
    V = gPi # initial guess
  
    function vfi(V)
        v_guess = fill(0.0, nPhi)
        for iter = 1:max_iter 
            v_guess .= V
            
            V = gPi + beta*max.(0,F_trans*v_guess) 

            sup = maximum(abs.(V - v_guess))  # check tolerance
            if sup<tol*maximum(abs.(V))
                if print_it; @printf("Iter: %1.0f. Tol. achieved: %1.2E \n", iter, sup); end;
            break; end
            if (iter == max_iter) & print_it; @printf("Max iterations achieved. VF did not converge: %1.2E\n", sup); end
            if (iter%iter_count==0) & print_it; @printf("Iter: %1.0f. Tol: %1.2E\n", iter, sup); end  
        end
        
        chi = fill(0.0, nPhi)
        chi[F_trans*V .< 0.0] .= 1.0 # recover exit policy function      
        return V, chi
    end
    
    V, chi = vfi(V)
    return V, chi, gN, gPi
end

function solve_price(param)
        
    function entry(p_guess)
        @unpack beta, c_e, G_prob, w = param
        V, chi, gN, gPi = solve_bellman(p_guess, param)
        excess_entry = beta*sum(V.*G_prob)-c_e*w
        return excess_entry
    end
    
    p = find_zero(entry, 5.0)
    V, chi, gN, gPi = solve_bellman(p, param)

    return (p = p, V = V, chi = chi, gN = gN, gPi = gPi)
end
    
function solve_m(param, solution)
    @unpack F_trans, nPhi, gPhi, alpha, G_prob, D_bar = param
    @unpack chi, gN, p = solution
    
    # construct transition probability:
    Phat = ((1 .- chi).*F_trans)'
    
    # invariant distribution is just a homogeneous function of M
    inv_dist(M) = M*inv(I - Phat)*G_prob # this is a function
    
    # supply: integrate the total production
    y = @. gPhi*gN^alpha
    supply = sum(inv_dist(1).*y) # just use the function with an arbitrary M

    # demand
    demand = D_bar/p
    
    # find mass of entrants (exploit linearity of the invariant distribution)
    M = demand/supply 
	mu = inv_dist(M)
	
    return M, mu
end

function ModelStats(param, sol_price, M, mu)
    @unpack F_trans, nPhi, gPhi, alpha, G_prob, D_bar = param
	@unpack gN, p, chi, gPi = sol_price
	
	# productivity distribution
	pdf_dist = mu./sum(mu)
	cdf_dist = cumsum(pdf_dist)
	
	# employment distribution
	emp_dist = mu.*gN
	pdf_emp = emp_dist./sum(emp_dist)
	cdf_emp = cumsum(pdf_emp)
	
	# exit productivity
	cut_index = findfirst(chi .== 0)
	phicut=	param.gPhi[cut_index]

	# stats
	avg_firm_size = sum(emp_dist)/sum(mu)
	exit_rate = M/sum(mu)
	Y = sum((gPhi.*gN.^alpha).*mu) ## agg production
	emp_prod = sum(emp_dist) # employment used in production
	Pi =  sum(gPi.*mu)   # aggregate profits
	
	# employment share
	#size_array = [10, 20, 50, 100, 500]
	
	return (pdf_dist, cdf_dist, pdf_emp, cdf_emp, avg_firm_size, exit_rate, Y, emp_prod, phicut, Pi)
	
end

### PS. The function is not considering the corner solution where there is no entry. Change the function later.

function SolveModel(param)
	
	# Solve For Prices
	sol_price =solve_price(param)  
	M, mu = solve_m(param, sol_price)
	
	if M<=0;
		println("Warning: No entry, eq. not found.")
	end
	
	stats = ModelStats(param, sol_price, M, mu)
	return (sol_price, M, mu, stats)
end

#begin
	solution_info = SolveModel(set_par(; alpha = 2/3))
	
	function printStats(solution_info)
		sol_price = solution_info[1]
		M = solution_info[2]
		mu = solution_info[3]
		stats = solution_info[4]
	
		p = sol_price.p
		avg_firm_size = stats[5]
		exit_rate = stats[6]
		Y = stats[7]
		cutoff = stats[9]
		Pi = stats[10]
	
		#with_terminal() do # this is just to print stuff in the notebook
		#	println("Price: $p")
		#	println("Avg. Firm Size: $avg_firm_size")
		#	println("Exit/entry Rate: $exit_rate")
		#	println("Productivity Cutoff: $cutoff")	
		#	println("Aggregate Output: $Y")	
		#	println("Aggregate Profits: $Pi")	
		#end
	end
	
	printStats(solution_info)
	
#end


	#begin
	@unpack gN = solution_info[1]
	
	plot(gN, [solution_info[4][1], solution_info[4][3]], label=["Firm Share" "Emp. Share" ], lw = 3)
	title!("Stationary Distribution")
	xlims!(0,2000) # little mass over 50
	#ylims!(0,0.05) # Ps the mass of constrained HH is about 20% of individuals in S)1
	xlabel!("Size")
	ylabel!("density")
	
	
#end

# increase Entry cost: 
#begin
	solution_info2 = SolveModel(set_par(; c_e = 60))
	printStats(solution_info2)
	#solution_info2[4]
#end

# increase fixed cost: 
#begin
#	solution_info3 = SolveModel(set_par(; c_f = 30))
#	printStats(solution_info3)
#end