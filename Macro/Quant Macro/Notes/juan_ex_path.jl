k_ss = ((1/beta - (1 - delta)) / (A*theta) )^(1/(theta-1))
#(Fix) k = (0.25:0.75/(size-1):1).^(1) [start; step size;stop]
Av = A*ones(1,100)
Av = Av'
k = k'
k_ini = k
k_old = k

k_all = []

dist=1
iter=1

#options=optimset('Display', 'iter');

while (dist > 0.001 && iter < 1000)
  k_all = [k_all k_old]

for i=1:length(k)-2

  f(y) = LHS(euler)-RHS(euler)
  [out,fval] = fsolve(f,k(i,1),options)
  k(i+1,1) = out
  if(abs(fval) > 0.001)
    disp('problem');
    fval
    #pause
  end
  #stuff?
end
