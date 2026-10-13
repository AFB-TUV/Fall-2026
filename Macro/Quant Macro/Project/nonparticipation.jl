using LinearAlgebra
#parameters
#impatience
const β = 0.95
#bargaining power
const α = 0.72
#firing rate
const λ = 0.145
#hiring rates
const q1 = 0.45
#not used because difference is not large
const q2 = 0.05
const q3 = 0.95
#transfer to bad state
const ξ1 = 0.30
#transfer to good state during employment
const ξ2e = 0.20
#transfer to good state during unemployment
const ξ2u = 0.05
#benefits
const b1 = 1.0
const b2 = 2.0
const b3 = 4.0
const b4 = 1.9
#human capital
const hu = 5.0
const hd = 2.0
#simulation parameters
const tol = 0.000001
const hn = 0.00001
const dampen = 0.2
#production function
f(h) = h^0.7
#value functions
Jf_u(w,Jf_up,Jf_dp) = f(hu) - w + β*(1-λ) * (1-ξ1) * Jf_up + β*(1-λ)*ξ1*Jf_dp
Jf_d(w,Jf_up,Jf_dp) = f(hd) - w + β*(1-λ) * (1-ξ2e) * Jf_dp + β*(1-λ)*ξ2e*Jf_up
Je_u(w,Je_up, Je_dp, Ju_up, Ju_dp) = w + β*(1-λ)*( (1-ξ1)*Je_up + ξ1*Je_dp ) + β*λ*( (1-ξ1)*Ju_up + ξ1*Ju_dp )
Je_d(w,Je_up, Je_dp, Ju_up, Ju_dp) = w + β*(1-λ)*( (1-ξ2e)*Je_dp + ξ2e*Je_up ) + β*λ*( (1-ξ2e)*Ju_dp + ξ2e*Ju_up )
Ju_u(w,Je_up, Je_dp, Ju_up, Ju_dp,χ,b,q) = b + β*(1-q)*( (1-ξ1)*Ju_up + ξ1*Ju_dp ) + β*χ*q*( (1-ξ1)*Je_up + ξ1*Je_dp )
Ju_d(w,Je_up, Je_dp, Ju_up, Ju_dp,χ,b,q) = b + β*(1-q)*( (1-ξ2u)*Ju_dp + ξ2u*Ju_up ) + β*χ*q*( (1-ξ2u)*Je_dp + ξ2u*Je_up )
#Nash Bargaining FOC
nb(Jf,Je,Ju) = (1-α)/α*(Je-Ju) - Jf

#helper function for calculating wages from bargaining
function y(sol::NTuple{6,Float64})
  local y = [ nb(sol[1],sol[3],sol[5]) nb(sol[2],sol[4],sol[6]) ]
  return y
end

#iterating the value function
function itVal(wu::Number, wd::Number,χ1::Number,χ2::Number,b::Number,q::Number)
  local Jfu = 0
  local Jfd = 0
  local Jeu = 0
  local Jed = 0
  local Juu = 0
  local Jud = 0

  local Jfup = 0
  local Jfdp = 0
  local Jeup = 0
  local Jedp = 0
  local Juup = 0
  local Judp = 0
  #loop until fixed-point convergence
  while true
    #update value functions
    Jfu = Jf_u(wu, Jfup, Jfdp)
    Jfd = Jf_d(wd, Jfup, Jfdp)
    Jeu = Je_u(wu, Jeup, Jedp, Juup, Judp)
    Jed = Je_d(wd, Jeup, Jedp, Juup, Judp)
    Juu = Ju_u(wu, Jeup, Jedp, Juup, Judp, χ1,b,q)
    Jud = Ju_d(wd, Jeup, Jedp, Juup, Judp, χ2,b,q)

    newvec = [Jfu; Jfd; Jeu; Jed; Juu; Jud]
    oldvec = [Jfup; Jfdp; Jeup; Jedp; Juup; Judp]
    #check for convergence
    if abs(maximum(newvec - oldvec)) < tol
      break
    else
      #get ready for another iteration
      Jfup = 0
      Jfdp = 0
      Jeup = 0
      Jedp = 0
      Juup = 0
      Judp = 0

      Jfup = Jfu
      Jfdp = Jfd
      Jeup = Jeu
      Jedp = Jed
      Juup = Juu
      Judp = Jud

      Jfu = 0
      Jfd = 0
      Jeu = 0
      Jed = 0
      Juu = 0
      Jud = 0

    end
  end

  return Jfu, Jfd, Jeu, Jed, Juu, Jud
end
#Newton's method in 2D
function Newton_nb(wu::Number, wd::Number,χ1::Number,χ2::Number,b::Number,q::Number)
  local wnew = [0.0; 0.0]
  local wold = [wu; wd]
  #initialize direction
  local dev = diagm(0 => hn*[wu; wd])

  while true
    #get value functions at current point, steps in each direction
    sol = itVal(wold[1], wold[2], χ1, χ2,b,q)
    solhup = itVal(wold[1]+hn,wold[2], χ1, χ2,b,q)
    solhdp = itVal(wold[1],wold[2]+hn, χ1, χ2,b,q)
    solhum = itVal(wold[1]-hn,wold[2], χ1, χ2,b,q)
    solhdm = itVal(wold[1],wold[2]-hn, χ1, χ2,b,q)
    #get value functions
    yold = y(sol)
    yup = y(solhup)
    ydp = y(solhdp)
    yum = y(solhum)
    ydm = y(solhdm)
    #approximate Jacobian
    dy0 = zeros(2,2)
    dy0[:,1] = (yup - yum)/(2*dev[1,1])
    dy0[:,2] = (ydp - ydm)/(2*dev[2,2])
    #update wages
    d0 = (-yold)/dy0
    wnew = (1-dampen)*wold + dampen*(wold + transpose(d0))
    ynew = y(itVal(wnew[1],wnew[2],χ1,χ2,b,q))

    #debugging!
    #println("New: ", wnew)
    #println(" ")
    #println("Old: ", wold)
    #println(" ")

    if abs(maximum(ynew - yold)) < tol
      break
    else
      wold = [0;0]
      wold = wnew
      wnew = [0;0]
    end

  end

  return wnew[1],wnew[2]
end
#solve the Bellman equation
function solveBell(wu::Number,wd::Number,b::Number,q::Number)
  #start with all hired
  local χ1o = 1
  local χ2o = 1

  local wunew = 0
  local wdnew = 0
  #solve
  sol1 = Newton_nb(wu,wd,χ1o,χ2o,b,q)
  wunew = sol1[1]
  wdnew = sol1[2]
  val1 = itVal(wunew,wdnew,χ1o,χ2o,b,q)

  local χ1 = 0
  local χ2 = 0
  #change indicators if not correct
  if val1[3]>val1[5]
    χ1 = χ1+1
  end

  if val1[4]>val1[6]
    χ2 = χ2+1
  end

  #println("Old: ", χ1o, " ", χ2o)
  #println(" ")
  #println("New: ", χ1, " ", χ2)
  #println(" ")
  #println("Wages: ", wunew, " ", wdnew)

  if χ1 == χ1o & χ2 == χ2o
    return wunew, wdnew, χ1o, χ2o
  else
    #if indicators are wrong, change and redo.
    #from trying, the second time always works
    χ1o = χ1
    χ2o = χ2

    sol2 = Newton_nb(wu,wd,χ1o,χ2o,b,q)
    wunew = sol2[1]
    wdnew = sol2[2]
    #val2 = itVal(wunew,wdnew,χ1o,χ2o,b)
    return wunew, wdnew, χ1o, χ2o
  end

end
#calculate the distribution of employment and h
function calDist(χ1::Number, χ2::Number,q)
  #calculate transition rates
  uuTud = ξ1*(1-q*χ1)
  uuTeu = (1-ξ1)*q*χ1
  uuTed = ξ1*q*χ1
  uuTuu = 1 - uuTud - uuTeu - uuTed

  udTuu = ξ2u*(1-q*χ2)
  udTed = (1-ξ2u)*q*χ2
  udTeu = ξ2u*q*χ2
  udTud = 1-udTuu-udTed-udTeu

  euTuu = (1-ξ1)*λ
  euTud = ξ1*λ
  euTed = ξ1*(1-λ)
  euTeu = 1 - euTuu - euTud - euTed

  edTuu = ξ2e*λ
  edTud = (1-ξ2e)*λ
  edTeu = ξ2e*(1-λ)
  edTed = 1 - edTuu - edTud - edTeu
  #construct transition matrix
  M = [ uuTuu uuTud uuTeu uuTed;
        udTuu udTud udTeu udTed;
        euTuu euTud euTeu euTed;
        edTuu edTud edTeu edTed  ]
  #stationary distribution
  MM = M^1000

  return M, MM[1,:]
end
#equilibria
eq1 = solveBell(f(5.0),f(2.0),b1,q1)
deq1 = calDist(eq1[3],eq1[4],q1)
eq2 = solveBell(f(5.0),f(2.0),b2,q1)
deq2 = calDist(eq2[3],eq2[4],q1)
eq3 = solveBell(f(5.0),f(2.0),b3,q1)
deq3 = calDist(eq3[3],eq3[4],q1)
eq4 = solveBell(f(5.0),f(2.0),b4,q1)
deq4 = calDist(eq4[3],eq4[4],q1)
