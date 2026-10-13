@ ------------------------------------ SOE.g ---------------------------

   10 January 2008
   Alfred Maussner

   Purpose: Solve the model of a small open economy of
            Heer and Maussner, 2nd edition, Section 3.3.2
            via the deterministic extended path method.
            
   Remark: the program will terminate with undefined symbol error messages
           unless you include references to the procedures FixvMN2
           and GraphSettings or use the #include command to
           read in the Gauss source files Derivatives, Nleq, and Tools.
           
           See the Readme_Ch3.txt file for more information.

-----------------------------------------------------------------------  @


new;

// The library file pgraph is required for plotting with the Gauss publication quality graphics.
// You may wish to add your user library or any other library which stores information about
// the locaation of the procedures called by this program, in particular, FixvMN2 and GraphSettings.

library pgraph;

/* Parameters of the model */
a=1.005;
alpha=0.27;
beta=0.994;
eta=2;
delta=0.011;
nstar=0.13;
nu=5;
zeta=1/30;
rho_z=0.90;
sigma_z=0.0072;
rho_r=0.90;
sigma_r=0.01;
bstar=0;

/* Parameter for the algorithm */
Tend=150;

/* Length of impulse response to compute */
nobs=30;

/* Determination of the parameters of the adjustment cost equation */
a1=(a+delta-1)^zeta;
a2=(a+delta-1)*(zeta/(zeta-1));

/* Stationary Solution */
rstar=((a^eta)/beta)-1;
yk=(a^eta-beta*(1-delta))/(alpha*beta);
kn=yk^(1/(alpha-1));
kstar=nstar*kn;
ystar=(nstar^(1-alpha))*(kstar^alpha);
theta=(1-alpha)*(ystar/nstar)*(nstar^(-nu));
istar=(a+delta-1)*kstar;
cstar=ystar-istar;
qstar=1;
lstar=(cstar-(theta/(1+nu))*(nstar^(1+nu)))^(-eta);
wstar=(1-alpha)*(ystar/nstar);
zstar=1;
hstar=ystar-cstar-istar;

test=(a1/(1-zeta))*(istar/kstar)^(1-zeta) + a2;

/* Compute paths for z and r */
zpath=zeros(Tend+1,1);
rpath=zeros(Tend+1,1);

/* Compute impulse responses for a one-time productivty shock */
zpath[1]=sigma_z;
for t (1,Tend,1);
    zpath[t+1]=rho_z*zpath[t];
endfor;
zpath=exp(zpath);
rpath=ones(Tend+1,1).*rstar;

goto next;

/* Compute the impulse response for a one-time change of the real interest rate */
zpath=ones(Tend+1,1).*zstar;
rpath=ones(Tend+1,1).*rstar;
rpath[1]=1.01*rstar;

goto next;
/* Compute the impulse response for an autocorrleated real interest rate shock */
rpath[1]=1.01*rstar;
for t (1,Tend,1);
    rpath[t+1]=(1-rho_r)*rstar + rho_r*rpath[t];
endfor;
next:

/* Set initial conditions */
k0=kstar;
b0=bstar;

x0=ones(Tend-1,1).*kstar;
x0=x0|ones(Tend-1,1).*bstar;
x0=x0|ones(Tend+1,1).*qstar;
x0=x0|lstar;

/* Set bounds */
bounds=       (ones(Tend-1,1).*0.00001)~(ones(Tend-1,1).*(10*kstar));
if bstar<0;
    bounds=bounds|(ones(Tend-1,1).*(10*bstar))~(ones(Tend-1,1).*(-10*bstar));
else;
    bounds=bounds|(ones(Tend-1,1).*(-10*bstar))~(ones(Tend-1,1).*(10*bstar));
endif;
bounds=bounds|(ones(Tend+1,1).*0.00001)~(ones(Tend+1,1).*(10*qstar));
bounds=bounds|(0.00001~10*lstar);

@zpath=ones(Tend+1,1);@
@test=Sys(x0);@

/* Find Solution */
DosWinOpen("SOE",0|0|15|1|1);
{x1,crit}=FixvMN2(x0,bounds,&Sys);
crit;

if crit[1]==0;
    /* Plot time path of solutions */
    kvec=k0|x1[1:Tend-1]|kstar;
    bvec=b0|x1[Tend:2*Tend-2]|x1[2*Tend-3];
    qvec=x1[2*Tend-1:3*Tend-1];
    lvec=zeros(rows(qvec),1);
    lvec[1]=x1[3*Tend];
    for t (1,rows(lvec)-1,1);
        lvec[t+1]=((a^eta)/beta)*(lvec[t]/(1+rpath[t+1]));
    endfor;
    Graphsettings;
    _plwidth=7;
    _ptek="Fig1.tkf";
    begwind;
    window(2,2,1);
    t=seqa(1,1,rows(kvec));
    title("Capital Stock");
    xy(t,kvec~ones(rows(kvec),1).*kstar);
    nextwind;
    title("Net Foreign Debt");
    xy(t,bvec~ones(rows(bvec),1).*bstar);
    nextwind;
    title("Relative Price of Capital");
    xy(t,qvec~ones(rows(qvec),1).*qstar);
    nextwind;
    title("Marginal Utility of Consumption");
    xy(t,lvec~ones(rows(lvec),1).*lstar);
    endwind;
    
    /* Compute time path of variables */    
    kvec=kvec[1:nobs];
    bvec=bvec[1:nobs];
    qvec=qvec[1:nobs];
    lvec=lvec[1:nobs];
    cvec=zeros(nobs,1);    
    yvec=cvec;
    nvec=cvec;
    wvec=cvec;
    ivec=cvec;
    zvec=cvec;
    hvec=cvec;

    for t (1,nobs,1);
        nvec[t]=((1-alpha)/theta)*zpath[t]*(kvec[t]^alpha);
        nvec[t]=(nvec[t]^(1/(alpha+nu)));        
        cvec[t]=(lvec[t]^(-1/eta))+(theta/(1+nu))*(nvec[t]^(1+nu));
        ivec[t]=((a1*qvec[t])^(1/zeta))*kvec[t];
        yvec[t]=zpath[t]*(nvec[t]^(1-alpha))*(kvec[t]^alpha);
        wvec[t]=(1-alpha)*(yvec[t]/nvec[t]);
        hvec[t]=yvec[t]-cvec[t]-ivec[t];
    endfor;
    for t (1,nobs,1);
        nvec[t]=(nvec[t]-nstar)/nstar;
        lvec[t]=(lvec[t]-lstar)/lstar;
        cvec[t]=(cvec[t]-cstar)/cstar;
        ivec[t]=(ivec[t]-istar)/istar;
        yvec[t]=(yvec[t]-ystar)/ystar;
        wvec[t]=(wvec[t]-wstar)/wstar;
        qvec[t]=(qvec[t]-qstar)/qstar;
        zvec[t]=(zpath[t]-zstar)/zstar;
        kvec[t]=(kvec[t]-kstar)/kstar;        
    endfor;
             

    /* Plot impulse responses */    
    t=seqa(1,1,nobs+1);
    GraphSettings;
    _ptek="Fig2.tkf";
    begwind;
    window(2,2,1);
    _pline=(1~6~0~0~nobs+1~0~0~5);
    _plegstr="Productivity";
    _plegctl=(1~7~10~0.5*maxc(zvec));
    zvec=0|zvec;
    xy(t,zvec);
    nextwind;
    _plegstr="Real Wage\000Hours";
    _plegctl=(1~7~10~0.5*maxc(maxc(wvec~nvec)));
    wvec=0|wvec;nvec=0|nvec;
    xy(t,wvec~nvec);
    nextwind;
    _plegstr="Output\000Consumption\000Investment";
    _plegctl=(1~7~10~0.5*maxc(maxc(yvec~cvec~ivec)));
    yvec=0|yvec; cvec=0|cvec; ivec=0|ivec;
    xy(t,yvec~cvec~ivec);
    nextwind;
    ylabel("Deviation from Trend");
    _plegstr="Trade Balance";
    _plegctl=(1~7~10~0.5*maxc(hvec));
    hvec=0|hvec;
    xy(t,hvec);
    endwind;

endif;

ende:
end;

/* ----------------------------- Procedures ------------------------ */

/* phi(x): defines the adjustment cost function phi(i/k) */
fn phi(x)=(a1/(1-zeta))*(x^(1-zeta))+a2;

/* Sys(x): the zero of this system determines the deterministic extended path.
**         External at zpath, rpath, tend, k0, kstar, and b0, as well as the parameters of the model. 
*/

proc(1)=Sys(x);

    local kt, bt, qt, n0, n1, c0, c1, i0, i1, l0, l1, fx, j, t;

      fx=zeros(rows(x),1);

      kt=k0|x[1:Tend-1]|kstar;
      bt=b0|x[Tend:2*(Tend-1)]|x[2*(Tend-1)];  @ implies b[T-1]=b[T] @    
      qt=x[2*Tend-1:3*Tend-1];
      l0=x[3*Tend];

      n0=((1-alpha)/theta)*zpath[1]*(kt[1]^alpha);
      n0=n0^(1/(alpha+nu));
      c0=(l0^(-1/eta))+(theta/(1+nu))*(n0^(1+nu));
      i0=((a1*qt[1])^(1/zeta))*kt[1];

    j=1;
    for t (1,Tend,1);
        n1=((1-alpha)/theta)*zpath[t+1]*(kt[t+1]^alpha);
        n1=n1^(1/(alpha+nu));
        l1=((a^eta)/beta)*(l0/(1+rpath[t+1]));
        c1=(l1^(-1/eta))+(theta/(1+nu))*(n1^(1+nu));
        i1=((a1*qt[t+1])^(1/zeta))*kt[t+1];                
        
        fx[j  ]=a*kt[t+1]+(delta-1)*kt[t] - phi(i0/kt[t])*kt[t];
        fx[j+1]=a*bt[t+1]-zpath[t]*(n0^(1-alpha))*(kt[t]^alpha)-(1+rpath[t])*bt[t]+c0+i0;
        fx[j+2]=qt[t]-(1/(1+rpath[t+1]))*(alpha*zpath[t+1]*(n1^(1-alpha))*(kt[t+1]^(alpha-1)) 
                                         + qt[t+1]*(1-delta+phi(i1/kt[t+1])) - (i1/kt[t+1]));          
             n0=n1;
             l0=l1;
             c0=c1;
             i0=i1;
              j=j+3;
    endfor;

    retp(fx);

endp;
    

 
