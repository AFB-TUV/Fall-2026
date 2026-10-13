@ --------------------------------------- Ramsey2c.g ---------------------


   Alfred Maussner
   05 January 2008
   28 September 2011 (last bugfix)

   Purpose: Solve the stochastic infinite-horizon Ramsey model
            of Heer and Maussner, 2nd Edition, by using the method of
            deterministic extended path in Section 3.2.

            u(c)=(c^(1-eta)-1)/(1-eta)
            F(N,K)=K^alpha
            C=F(N,K)+(1-delta)K - K'

            if eta=1 and delta=1 the program plots the
            time path from the analytic solution together
            with the time path obtained from the deterministic
            extended path method.

  external procedures:

    - FixvMN2 (in NLEQ.src from the Gauss programs to Chapter 11)

-------------------------------------------------------------------------  @

new;

/* Parameters */
alpha=0.27;
beta=0.994;
delta=0.011; delta=1;
eta=2.0;     eta=1;
rho=0.90;
sigma=0.0072;
Tmax=150;

nobs=100;   @ the length of the comuted time path in the simulation of the model @


zt=zeros(Tmax+1,1);      @ a vector of Tmax elements, external to GetKC and Sys @
x0=zeros(Tmax,1);      @ this vector stores the inital values for the non-linear equations solver @
x1=zeros(Tmax,1);      @ this vector stores the solution @
bounds=zeros(Tmax,2);  @ matrix of lower and upper bounds for non-linear equations solver @


/* Stationary Solution */
kstar=((1-beta*(1-delta))/(alpha*beta))^(1/(alpha-1));
cstar=kstar^alpha - delta*kstar;

/* Fix lower and upper bounds */
bounds[.,1]=ones(Tmax,1).*(0.0001);
bounds[.,2]=ones(Tmax,1).*(10*kstar);

/* Initial values */
x0=ones(Tmax,1).*kstar;
k0=kstar; @ k0 is also an external matrix to GetKC and Sys @

/* Compute a long path, given a path of productivity levels */
eps=sigma*rndn(nobs+1,1);
zvec=zeros(nobs+1,1);

zvec[1]=eps[1];
for t (1,nobs,1);
    zvec[t+1]=rho*zvec[t]+eps[t+1];
endfor;

kvec=zeros(nobs+1,1);
cvec=zeros(nobs,1);

kvec[1]=kstar;
i=DosWinOpen("Ramsey3c",0|0|15|1|1);
_MNR_PRINT=0;

for t (1,nobs,1);
    locate 2,2;?"Period t= " ftos(t,"*.*lf",5,0);      
    {k1,c0}=GetKC(zvec[t],kvec[t]);
    if ismiss(k1); ?"Not able to solve for k1 "; goto ende; endif;
    kvec[t+1]=k1;
    cvec[t]=c0;
    x0=x1;  @ use the previous solution as starting value for the next computation @

endfor;

if eta==1 and delta==1;
    kvec2=zeros(nobs+1,1);
    kvec2[1]=kstar;    
    for t (1,nobs,1);
        kvec2[t+1]=alpha*beta*exp(zvec[t])*(kvec2[t]^alpha);        
    endfor;
    GraphSettings;
    _plwidth=7;
    _paxht=0.20;
    xlabel("Period");
    ylabel("Capital Stock");
    _plegstr="Forward Iteration\000Analytic Solution";
    _plegctl=(1~5~10~maxc(kvec));
    _pltype=6;
    _pcolor={0,0};
    _pltype={6,4};
    xy(seqa(0,1,rows(kvec)),kvec~kvec2);
    output file=Ramsey3c.txt reset;
    MyDate;
    test=(kvec-kvec2)./kvec2;
    ?"Maximum relative distance between the two solutions";
    ?maxc(abs(test));
    output off;
endif;
ende:

end;

/* ---------------------------------- Procedures ------------------------- */

/* GetKC: determines the time path of the model under the assumption
**        that no further shock occurs.
**        The program returns the next-period capital stock k1
**        and consumption obtained from that path.
*/

proc(2)=GetKC(z,k);

    local t, crit, k1, c0;
    external matrix zt, x0, x1, bounds, k0;

    @ compute the exptected path of the log of productivity @
    zt[1]=z;
    for t (1,Tmax,1);
        zt[t+1]=rho*zt[t];
    endfor;
    k0=k;
    {x1,crit}=FixvMN2(x0,bounds,&Sys);    
    if crit/=0;
        retp(miss(1,1),miss(1,1));
    else;
        k1=x1[1];
        c0=exp(z)*(k^alpha)+(1-delta)*k - k1;
        retp(k1,c0);
    endif;

endp;

/* Sys(x): returns the rhs of the system of Tmax equations
**         that constitute the Eulerequations of the model.
**         See equtions of the book.
**         The vector zt must be given in the main program.
*/

proc(1)=Sys(k);

    local fx, x, c0, c1, t;

    fx=zeros(rows(k),1);
    x=k0|k|kstar;

    c0=(exp(zt[1])*x[1]^alpha)+(1-delta)*exp(zt[1])*x[1] - x[2];
    if c0<0; fx[1]=miss(1,1); retp(fx); endif;
  
    for t (1,tmax,1);
        c1=(exp(zt[t+1])*x[t+1]^alpha)+(1-delta)*exp(zt[t+1])*x[t+1] - x[t+2];
        if c1<0; fx[t]=miss(1,1); retp(fx); endif;
        fx[t]=(c0^(-eta)) - beta*(c1^(-eta))*(1-delta+alpha*exp(zt[t+1])*(x[t+1]^(alpha-1)));
        c0=c1;
    endfor;

    retp(fx);

endp;
    
 