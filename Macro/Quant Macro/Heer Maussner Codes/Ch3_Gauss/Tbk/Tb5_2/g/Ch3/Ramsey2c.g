@ --------------------------------------- Ramsey2c.g ---------------------


   Alfred Maussner
   05 January 2008

   Purpose: Solve the deterministic infinite-horizon Ramsey model
            of Heer and Maussner, 2nd Edition, by using
            forward iteration as described in Section 3.1.

            u(c)=(c^(1-eta)-1)/(1-eta)
            F(N,K)=K^alpha
            C=F(N,K)+(1-delta)K - K'

-------------------------------------------------------------------------  @

new;

/* Parameters */
alpha=0.27;
beta=0.994;
delta=0.011;
eta=2.0;
Tmax=250;

/* Stationary Solution */
kstar=((1-beta*(1-delta))/(alpha*beta))^(1/(alpha-1));
cstar=kstar^alpha - delta*kstar;

/* Capital stock at t=0 */
rel=0.10;
k0=rel*kstar;
x0=ones(tmax,1).*kstar;
bounds=ones(tmax,1).*0.01~ones(tmax,1).*(10*kstar);
cls;
@{x1,crit}=eqSolve(&Sys,x0);@
@goto ende;@
cls;output file=Ramsey2.txt reset;
{x1,crit}=FixvMN2(x0,bounds,&Sys);
if crit[1]==0;
    do while rel>0.10;
        rel=rel-0.01;
        k0=rel*kstar;
        x0=x1;
        {x1,crit}=FixvMN2(x0,bounds,&Sys);
        if crit[1]/=0; break; endif;
    endo;
endif;
if crit[1]==0;
    GraphSettings;
    _plwidth=7;
    xlabel("Period");
    ylabel("Capital Stock");
    t=seqa(1,1,tmax);
    xy(t,x1);
else;
    ?"Not able to find solution";
endif;
ende:
end;

end;

/* ---------------------------------- Procedures ------------------------- */

/* Sys(x): returns the rhs of the system of Tmax equations
**         that constitute the Eulerequations of the model.
**         See equtions of the book.
*/
proc(1)=Sys(k);

    local fx, x, c0, c1, t;

    fx=zeros(rows(k),1);
    x=k0|k|kstar;

    c0=(x[1]^alpha)+(1-delta)*x[1] - x[2];
    if c0<0; fx[1]=miss(1,1); retp(fx); endif;
  
    for t (1,tmax,1);
        c1=(x[t+1]^alpha)+(1-delta)*x[t+1] - x[t+2];
        if c1<0; fx[t]=miss(1,1); retp(fx); endif;
        fx[t]=(c0^(-eta)) - beta*(c1^(-eta))*(1-delta+alpha*(x[t+1]^(alpha-1)));
        c0=c1;
    endfor;

    retp(fx);

endp;
    
 