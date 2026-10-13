@ ------------------------------------ SGNNI_a.g ----------------------------------

 
  12 February 2008
  Alfred Maussner

  Last changes: 21 November 2008, Alfred Maussner


  Purpose: Solve the stochastic infinite horizon Ramsey model
           of Heer and Maussner, 2nd ed. with a non-negativity
           constraint on investment via iterations over
           the policy function on a discrete grid.

                      1. one-period utility function
                      
                         u(c)=(c^(1-eta)/(1-eta))/(1-eta)
                         
                         or

                         u(C)=ln(c)
                         
                         for eta=1

                      2. technology Y = Z*K^alpha

                      3. transition equation

                         K' = Y + (1-delta)*K - C

                      4. Productivity Shock ln(Z[t]) = rho*ln(Z[t-1])+ eps[t];

  In this version of the program, the previous value function is used
  to initialize the value function for the next, finer grid.

  Different from a previous version of this program, the procedure SolveVI
  is able to implement the constraint for all settings of _VI_IP (0, 1, 2).


----------------------------------------------------------------------------------- @

new; @ clear memory @

/* Parameters of the model */
alpha=0.27;                   @ elasticity of production with respect to capital    @
delta=0.011;                  @ rate of capital depreciation                        @
beta=0.994;                   @ discount factor: 6.5% annual rate of return on capital@
eta=2;                        @ elasticity of marginal utility                    @
rho=0.90;
sigma=0.0072;
sigma=0.05;

/* Parameters of the algorithm */
nz=15;                        @ number of grid points for the productivity shock @
size=4.5;                     @ size of the grid for the productivity shock @
nk=250;                       @ number of grid points for the capital stock @
kmin_g=0.10;                  @ lower bound of the grid for the capital stock @
kmax_g=2.50;                  @ upper bound of the grid for the capital stock @
_VI_IP=1;                     @ =0 without interpolation, =1 linear interpolation, =2 cubic interpolation @
_VI_MPI=0;                    @ =0 without modified policy iteration, =1 with modified policy iteration @
_VI_MPI_K=30;                 @ number of iterations over the policy function @
_VI_nc=300;                    @ number of iterations with unchanged policy function before iterations are terminated @
_cumulative=1;                @ in successive computations, compute total run time, =0 do not @

/* Parameters for the computation of Euler equation residuals */
kmin_e=0.8;               @ kmin_e*kstar is the lower bound @
kmax_e=1.2;               @ kmax_e*kstar is the upper bound @
nobs_e=200;               @ the number of residuals to be computed @


/* Open file and write initial information */
output file=SGNNI_a.txt reset;
MyDate;
?"Parameters of the model:";
?"alpha = " ftos(alpha,"*.*lf",8,4);
?"beta  = " ftos(beta, "*.*lf",8,4);
?"delta = " ftos(delta,"*.*lf",8,4);
?"eta   = " ftos(eta,  "*.*lf",8,4);
?"rho   = " ftos(rho,  "*.*lf",8,4);
?"sigma = " ftos(sigma,"*.*lf",8,4);
?"";
?"Parameters of the algorithm:";
?"kmin_g= " ftos(kmin_g, "*.*lf",8,3);
?"kmax_g= " ftos(kmax_g, "*.*lf",8,3);
?"kmin_e= " ftos(kmin_e, "*.*lf",8,3);
?"kmax_e= " ftos(kmax_e, "*.*lf",8,3);
?"nobs_e= " ftos(nobs_e,"*.*lf",8,0);
?"nz    = " ftos(nz    ,"*.*lf",8,0);
?"";
?"eps   = " ftos(_VI_eps ,"*.*lf",7,4);
?"stop  = " ftos(_VI_nc  ,"*.*lf",7,0);
?"IP    = " ftos(_VI_IP ,"*.*lf",7,0);
?"MPI   = " ftos(_VI_MPI,"*.*lf",7,0);
if _VI_MPI;
?"K     = " ftos(_VI_MPI_K,"*.*lf",7,0);
endif;

output off;

/* Compute Markov chain approximation */
{zgrid,pmat}=MarkovAR(size,nz,rho,sigma);

zgrid=exp(zgrid);

zmin=zgrid[1];
zmax=zgrid[nz];

zmin_i=rho*ln(zmin)+sqrt(2)*sigma*(-1.650680123);
zmax_i=rho*ln(zmax)+sqrt(2)*sigma*(1.650680123);


kmin=((1-beta*(1-delta))/(alpha*beta*zmin))^(1/(alpha-1));
kmax=((1-beta*(1-delta))/(alpha*beta*zmax))^(1/(alpha-1));

/* Compute stationary solution of deterministic model and intialize the value function */
kstar=((1- beta*(1-delta))/alpha*beta)^(1/(alpha-1));  @ stationary capital stock                          @
cstar=kstar^alpha - delta*kstar;                       @ stationary level of consumption                   @


kmin_g=kmin_g*kmin;
kmax_g=kmax_g*kmax;

kmin_e=kmin_e*kstar;
kmax_e=kmax_e*kstar;

/* Vector with different values of nk */
nvec={250, 250};

/* Iterations of nvec start here */
lmax=rows(nvec);

policy=arrayinit(lmax|nobs_e|nobs_e,0);
  emat=arrayinit(lmax|nobs_e|nobs_e,0);

/* Different ways to initialize v0, uncomment to try the others */
nk=nvec[1];

v0=rf(1,kstar,kstar)/(1-beta);  @ stationary solution @
v0=ones(nk,nz).*v0;

@     v0=zeros(nk,nz); @        @ the zeros function  @
@    for i (1,nk,1); for j (1,nz,1); knext=kgrid[i]; v0[i,j]=rf(zgrid[j],kgrid[i],knext); endfor; endfor;@ @ maintain the given capital stock @

/* Iterations over different nk start here */
s2=0; @ stores time needed to obtain initial v from previous v @

tottime=zeros(lmax,1);

for l (1,lmax,1);
    nk=nvec[l];
    kgrid=seqa(kmin_g,(kmax_g-kmin_g)/(nk-1),nk);
    
    @ outcommend if necessary @
     if l==1; _VI_IP=0; else; _VI_IP=2; endif; 

    /* Solve for the policy function */
    s1=hsec;
    {v1,hmati}=SolveVI(beta,kgrid,zgrid,pmat,v0);
    s1=hsec-s1;
    
    if _VI_IP/=0;
     hmat=hmati;
    else;
        hmat=zeros(nk,nz);
        for i (1,nk,1);
            for j (1,nz,1);
                hmat[i,j]=kgrid[hmati[i,j]];
            endfor;
        endfor;
    endif;

    tottime[l]=s1+s2;
    test="nk= " $+ ftos(nk,"*.*lf",10,0) $+  " Run time= " $+ etstr(s1+s2);
    if _cumulative;
        test=test $+ "Cumulative run time= " $+ etstr(sumc(tottime[1:l]));    
    endif;
    output on;
    @?"_VI_IP= " ftos(_VI_IP,"*.*lf",2,0);@
    ?test;
    output off;
    /* computation of policy function */
   @ for i (1,nobs_e,1);@
        @for j (1,nobs_e,1);@
            @policy[l,i,j]=BLIP(kgrid,zgrid,hmat,kvec[i],zvec[j]);@
        @endfor;@
    @endfor;@
    /* New initial v0 */    
    if l<lmax;
        if nk==nvec[l+1];
            v0=v1;
        else;
            locate 15,5; ?"Compute new initial v0";
            s2=hsec;
            nk1=nvec[l+1];
            kgnew=seqa(kmin_g,(kmax_g-kmin_g)/(nk1-1),nk1);   
            kgnew[nk1]=kgrid[nk];kgnew[1]=kgrid[1]; 
            v0=zeros(nk1,nz);    
            for j (1,nz,1);
                locate 15,30;?ftos(j,"*.*lf",9,0);
                v0[.,j]=LIP(kgrid,v1[.,j],kgnew);
            endfor;
            s2=hsec-s2;
        endif;
    endif;
cp=zeros(nk,nz);
for i (1,nk,1);
    for j (1,nz,1);
        cp[i,j]=zgrid[j]*(kgrid[i]^alpha) + (1-delta)*kgrid[i]-hmat[i,j];
    endfor;
endfor;
GraphSettings;
surface(kgrid',zgrid,cp');

endfor;
cp=zeros(nk,nz);
for i (1,nk,1);
    for j (1,nz,1);
        cp[i,j]=zgrid[j]*(kgrid[i]^alpha) + (1-delta)*kgrid[i]-hmat[i,j];
    endfor;
endfor;
@GraphSettings;@
@surface(kgrid',zgrid,cp');@
@goto ende;@
kgrid_oc_20T=kgrid;
zgrid_oc_31=zgrid;
cp_oc_20T=cp;
@save kgrid_oc_20T, zgrid_oc_31, cp_oc_20T;@

@save emat, policy, kvec, zvec;@
ende:
end;
   

/* ----------------------------------- Subroutines ----------------------- */
   
proc(1)=rf(z,k1,k2);  @ defines the utility function @

  local c;

  @if (k2<(1-delta)*k1); ?"Error"; wait; endif;@
  c=z*(k1^alpha) + (1-delta)*k1 - k2;
  if c<0.0; retp(miss(1,1)); endif;
  
  if eta==1;
      c=ln(c);
  else;
     c=(c^(1-eta))/(1-eta);
  endif;
  

retp(c);

endp;

/* Policy function via linear interpolation: kgrid, zgrind, and hmat must be given as a global matrices */

proc(1)=PF(k,z);

    local knext;

    knext=BLIP(kgrid,zgrid,hmat,k,z);

    retp(knext);

endp;


proc(2)=SolveVI(beta,xvec,zvec,pmat,v0);
                   
   local i, j, l, v1, v2, v3, h1, h2, t, nx, nz, dv, eps1, di, nc, w, js, jmin, jmax, jl, ju,
         ax, cx, bx, umat;

   external matrix _VI_MPI, _VI_MPI_K, _VI_xvec, _VI_ymat, _VI_xex, _VI_zex, _VI_eps, _VI_nc, _VI_der;
   external proc rf;

/* Step 1: Initialize */
eps1=_VI_eps*(1-beta);      @ convergence criteria @


nx=rows(xvec);   @ number of grid points in xvec @
nz=rows(zvec);   @ number of grid points in zvec @

if _VI_IP==0; h1=ones(nx,nz); endif;            @ intial policy function @
h2=zeros(nx,nz); @ new policy function @ 
w =zeros(3,1);

v1=v0;              @ old policy function @
v2=zeros(nx,nz);    @ new policy function @
dv=1;
nc=0;

if _VI_IP/=0; _VI_xvec=xvec;_VI_zvec=zvec;_VI_ymat=v0; _VI_pmat=pmat; _VI_beta=beta; endif;
if _VI_IP==2; _VI_der=zeros(nx,nz); for j (1,nz,1); _VI_der[.,j]=CSpline(xvec,v1[.,j],1,0|0); endfor; endif;
if _VI_MPI==1; umat=zeros(nx,nz); v3=umat; endif;

/* Step 2: Iterate over the value function */
DosWinOpen("Value Function Iteration",0|0|15|1|1);
cls;
t=1;
do until (t>_VI_Max) or (dv<eps1) or (nc>_VI_nc); @ begin loop over value function @

    for j (1,nz,1);  @ begin loop over zvec @
        if _VI_IP/=0; _VI_zex=j; endif;
        js=1;
        for i (1,nx,1);  @ begin loop over xvec @
            if _VI_IP/=0; _VI_xex=xvec[i]; endif; 
            if _VI_BS;    
                jmin=js;
                do while xvec[jmin]<(1-delta)*xvec[i]; jmin=jmin+1; endo;                
                jmax=nx;
                do while (jmax-jmin)>2;       @ the next lines implement the binary search algorithm @
                    jl=floor((jmin+jmax)/2);                
                    ju=jl+1;                    
                    w[1]=rf(zvec[j],xvec[i],xvec[jl])+beta*(pmat[j,.]*(v1[jl,.]'));                    
                    w[2]=rf(zvec[j],xvec[i],xvec[ju])+beta*(pmat[j,.]*(v1[ju,.]'));                    
                    if w[2]>w[1]; jmin=jl; else; jmax=ju; endif;
                endo;
                w[1]=rf(zvec[j],xvec[i],xvec[jmin])+beta*(pmat[j,.]*(v1[jmin,.]'));
                if jmax>jmin;    w[2]=rf(zvec[j],xvec[i],xvec[jmin+1])+beta*(pmat[j,.]*(v1[jmin+1,.]')); else; w[2]=w[1]; endif;
                w[3]=rf(zvec[j],xvec[i],xvec[jmax])+beta*(pmat[j,.]*(v1[jmax,.]'));
                js=maxindc(w);
                if _VI_IP==0; v2[i,j]=w[js]; endif;               
                js=jmin+js-1;
               if _VI_MPI==1; umat[i,j]=rf(zvec[j],xvec[i],xvec[js]); endif;
            else;    
                jmin=js;
                w[1]=rf(zvec[j],xvec[i],xvec[jmin])+beta*(pmat[j,.]*(v1[jmin,.]'));         
                for jl (jmin+1,nx,1);                
                     w[2]=rf(zvec[j],xvec[i],xvec[jl])+beta*(pmat[j,.]*(v1[jl,.]'));
                    if w[2]<=w[1]; js=jl-1; break; else; w[1]=w[2]; endif;                
                endfor; 
            endif;
    
            /* The next lines implement linear interpolation between grid points */
            if _VI_IP/=0;
                if js==1;  /* boundary optimum, ax=bx=a[1]  */
                    ax=xvec[1];
                    bx=ax+eps1*(xvec[2]-xvec[1]);
                    cx=xvec[2];
                    if _rhs_bellman(j,xvec[i],bx)<_rhs_bellman(j,xvec[i],ax);
                        h2[i,j]=xvec[1];
                    else;                      
                        h2[i,j]=GSS(&_VI_valuefunction,xvec[1],xvec[2]);                    
                    endif;
                elseif js==nx;   /* boundary optimum, bx=cx=a[n] */
                    ax=xvec[nx-1];
                    cx=xvec[nx];
                    bx=cx-eps1*(xvec[nx]-xvec[nx-1]);
                    if bx>(1-delta)*xvec[i];
                        if _rhs_bellman(j,xvec[i],bx)<_rhs_bellman(j,xvec[i],cx);
                            h2[i,j]=xvec[nx];
                        else;
                            if xvec[nx-1]>=(1-delta)*xvec[i];
                                h2[i,j]=GSS(&_VI_valuefunction,xvec[nx-1],xvec[nx]);                    
                            else;
                                h2[i,j]=GSS(&_VI_valuefunction,(1-delta)*xvec[i],xvec[nx]);
                            endif;
                        endif;
                    else;
                        h2[i,j]=xvec[nx];
                    endif;
                else;
                    if xvec[js-1]>=(1-delta)*xvec[i];
                        h2[i,j]=GSS(&_VI_valuefunction,xvec[js-1],xvec[js+1]);         
                    else;
                        h2[i,j]=GSS(&_VI_valuefunction,(1-delta)*xvec[i],xvec[js+1]);
                    endif;
                endif;
    
                v2[i,j]=_rhs_bellman(j,xvec[i],h2[i,j]);    
            else;
                h2[i,j]=js;
            endif;
            
        endfor; @ end loop over xvec @
    endfor;    @ end loop over zvec @
    
    if _VI_IP==0;
        /* modified policy iteration */
        if _VI_MPI==1;            
            for l (1,_VI_MPI_K,1);                
                for i (1,nx,1);
                    for j (1,nz,1);
                        v3[i,j]=umat[i,j]+beta*(pmat[j,.]*(v2[h2[i,j],.]'));
                    endfor;
                endfor;
                v2=v3;
            endfor;
        endif;
        /* compute stopping criterium 2 */        
        di=sumc(sumc(h2 ./= h1));
        if di>=1; nc=0; else; nc=nc+1; endif; 
        h1=h2;
    endif;
    dv=maxc(maxc(abs(v2-v1))); 

    locate 5,5;
    ?"Iteration #= " ftos(t,"*.*lf",6,0);
    locate 6,5;
    ?"Largest element in v1-v0= " dv;
    if _VI_IP==0;
        locate 7,5; ?"# of indices that have changed: " di;
        locate 8,5; ?"# of consecutive iterations with constant policy function=" nc;
    endif;
    v1=v2;  
    if _VI_IP/=0; _VI_ymat=v1; endif;
    if _VI_IP==2; for j (1,nz,1); _VI_der[.,j]=CSpline(xvec,v1[.,j],1,0|0); endfor; endif;
    t=t+1;
endo;
if t>_VI_Max;
   locate 10,5;
   ?"Maximum number of iterations exceeded. Change _VI_Tmax!";
   locate 11,5;
   ?"The computed solution may be inaccurate.Press any key...";wait;
endif;
if _VI_IP==0;
    if minc(minc(h1))==1;  locate 12,5; ?"Policy function hits lower bound of grid"; endif;
    if maxc(maxc(h1))==nx; locate 13,5; ?"Policy function hits upper bound of grid"; endif;
else;
    if minc(minc(h2))==xvec[1];  locate 12,5; ?"Policy function hits lower bound of grid"; endif;
    if maxc(maxc(h2))==xvec[nx]; locate 12,5; ?"Policy function hits upper bound of grid"; endif;
endif;

retp(v2,h2);

endp;

@------------------------------------ MarkovAR ---------------------------------

   Alfred Maussner
   2004

   Purpose:  Approximate AR(1)-Process by Markov chain (Algorithm 9.2.1 of Heer and Maußner, 2005)
  
   Usage: {z,p}=MarkovAR(size,m,rho,sigma)
  
   Input: size:  scalar, the mulitple of the unconditional 
                 standard deviation of the AR(1) process
                 used to define the grid size
  
             m:  integer scalar, the number of grid points
   
           rho:  scalar, the autoregressive parameter of the process
  
         sigma:  scalar, the standard deviation of the innovations
  
   Output:   z:  m x 1 vector the grid approximating the process
  
             p:  m x m matrix of transition probabilities
----------------------------------------------------------------------------------- @
   
   
proc(2)=MarkovAR(size,m,rho,sigma);

  local p, zt, i, j, sigmaz, zbar;

  sigmaz=sqrt(sigma^2/(1-rho^2));
  zbar=size*sigmaz;  
  zt=seqa(-zbar,2*zbar/(m-1),m);

  p=zeros(m,m);
  i=1;
  do until i>m;

   p[i,1]   =cdfn((zt[1]-rho*zt[i]+(zt[2]-zt[1])/2)/sigma); 
   
   j=2;
   do until j>(m-1);
      p[i,j]=cdfn((zt[j]-rho*zt[i]+(zt[j]-zt[j-1])/2)/sigma)-
             cdfn((zt[j]-rho*zt[i]-(zt[j]-zt[j-1])/2)/sigma);
           j=j+1;
   endo;
   p[i,m]=1-sumc(p[i,1:m-1]');
   i=i+1;
endo;
              
retp(zt,p);

endp;

/* GraphSettings: Changes Gauss' default initialization of graphic routines
**
** Usage: GraphSettings;
*/
proc(0)=GraphSettings;

 external matrix  _pltype, _ptitlht, _pnumht,
         _paxht, _pmcolor, _plwidth, _pcsel, _pdate;
 external proc GraphSet;
 GraphSet;
   _pdate="";        @ do not plot date               @
  _pltype=6;         @ solid lines only               @
 _plwidth=4;         @ line width always 4            @
   _paxht=0.15;      @ size of axis labels            @        
  _pnumht=0.15;      @ size of axis numbering         @
 _ptitlht=0.20;      @ size of titles                 @
 _pmcolor={0,        @ color of axis: black           @
           0,        @ color of axis' numbers         @
           0,        @ color of x-axis label          @
           0,        @ color of y-axis label          @
           0,        @ color of z-axis label          @
           0,        @ color of title                 @
           0,        @ color of boxes                 @
           0,        @ color of date and author       @
           15};      @ background color:        white @
   _pcsel={0,        @ color of first line: black     @
           9,        @ color of second line: light blue @
           10,       @ color of third line:  light green @
           12,       @ color of fourth line: light red      @
           13,       @ color of fifth line:  light magenta @
           11,       @ color of sixth line: light cyan @
           6};       @ color of seventh line: brown @
   fonts("Simplex, Simgrma");
          
retp;
endp;           

