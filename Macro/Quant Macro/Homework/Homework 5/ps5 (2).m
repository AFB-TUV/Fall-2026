clear all
clc
n_l=25;
n_h=80;
n_ret=65;
sigma=2;
beta=0.92;
r=0.02;
R=1+r;
z=[log(15) log(80) log(500)];
Gamma=[0.6 0.35 0.05;0.1 0.8 0.1;0.05 0.35 0.6];
for n=n_l:n_h;
    f(n)=log(-1.0135+0.1086*n-0.001122*n^(2));
end
eps=[log(0.8) log(1.2)];
pr_eps=[0.5 0.5];
Eeps=pr_eps*eps';
% eps=[Eeps Eeps];
lee=length(eps);
stat_z=[1/3 1/3 1/3]*Gamma^(1000);
Ez=stat_z*z';
% z=[log(500) log(500) log(500)];
lez=length(z);
a_low=0;
a_high=5000;
lea=450;
a=linspace(a_low,a_high,lea);
V=zeros(lea,lez,lee,n_h);
Ind=zeros(lea,lez,lee,n_h);
%value function at the age of 80
for u=1:lee;
    for j=1:lea;
        for i=1:lez;
            Rz(i)=0.7*exp(z(i)+f(n_ret));
            V(j,i,u,n_h)=(1/(1-sigma))*(R*a(j)+Rz(i))^(1-sigma);
            if R*a(j)+Rz(i)<0;
               V(j,i,u,n_h)=-999999-99*abs(R*a(j)+Rz(i));
            end
        end
    end
end

%value function from 66 to 79
for u=1:lee;
    for n=1:(n_h-n_ret-1);
        for j=1:lea;
            for i=1:lez;
                for jj=1:lea;
                    V1m(j,jj,i,u,n_h-n)=(1/(1-sigma))*(R*a(j)+Rz(i)-a(jj))^(1-sigma)+beta*V(jj,i,u,n_h-(n-1));
                    if R*a(j)+Rz(i)-a(jj)<0;
                       V1m(j,jj,i,u,n_h-n)=-999999-99*abs(R*a(j)+Rz(i)-a(jj));
                    end
                end
                [V(j,i,u,n_h-n),Ind(j,i,u,n_h-n)]=max(V1m(j,:,i,u,n_h-n));
            end
        end
    end
end
%value function at 65
for u=1:lee;
    for j=1:lea;
        for i=1:lez;
            for jj=1:lea;
                V1m(j,jj,i,u,n_ret)=(1/(1-sigma))*(R*a(j)+exp(z(i)+f(n_ret)+eps(u))-a(jj))^(1-sigma)+beta*V(jj,i,u,n_ret+1);
                if R*a(j)+exp(z(i)+f(n_ret)+eps(u))-a(jj)<0;
                   V1m(j,jj,i,u,n_ret)=-999999-99*abs(R*a(j)+exp(z(i)+f(n_ret)+eps(u))-a(jj));
                end
            end
            [V(j,i,u,n_ret),Ind(j,i,u,n_ret)]=max(V1m(j,:,i,u,n_ret));
        end
    end
end
%value function from 25 to 64
for u=1:lee;
    for n=1:(n_ret-n_l);
        for j=1:lea;
            for i=1:lez;
                for jj=1:lea;
                    V1m(j,jj,i,u,n_ret-n)=(1/(1-sigma))*(R*a(j)+exp(z(i)+f(n_ret-n)+eps(u))-a(jj))^(1-sigma)+beta*(Gamma(i,1)*(0.5*V(jj,1,1,n_ret-(n-1))+0.5*V(jj,1,2,n_ret-(n-1)))+Gamma(i,2)*(0.5*V(jj,2,1,n_ret-(n-1))+0.5*V(jj,2,2,n_ret-(n-1)))+Gamma(i,3)*(0.5*V(jj,3,1,n_ret-(n-1))+0.5*V(jj,3,2,n_ret-(n-1))));
                    if R*a(j)+exp(z(i)+f(n_ret-n)+eps(u))-a(jj)<0;
                       V1m(j,jj,i,u,n_ret-n)=-999999-99*abs(R*a(j)+exp(z(i)+f(n_ret-n)+eps(u))-a(jj));
                    end
                end
                [V(j,i,u,n_ret-n),Ind(j,i,u,n_ret-n)]=max(V1m(j,:,i,u,n_ret-n));
            end
        end
    end
end

%consumption and policy function
C=zeros(lea,lez,lee,n_h);
Ap=zeros(lea,lez,lee,n_h);
%at age 80
for j=1:lea;
    for i=1:lez;
        C(j,i,:,n_h)=R*a(j)+Rz(i);
        Ap(j,i,:,n_h)=0;
    end
end
%at age 66 to 79
for n=1:(n_h-n_ret-1);
    for j=1:lea;
        for i=1:lez;
            C(j,i,:,n_h-n)=R*a(j)+Rz(i)-a(Ind(j,i,:,n_h-n));
            Ap(j,i,:,n_h-n)=a(Ind(j,i,:,n_h-n));
        end
    end
end
%at age 65
for j=1:lea;
    for i=1:lez;
        for u=1:lee;
            C(j,i,u,n_ret)=R*a(j)+exp(z(i)+f(n_ret)+eps(u))-a(Ind(j,i,u,n_ret));
            Ap(j,i,u,n_ret)=a(Ind(j,i,u,n_ret));
        end
    end
end
%at age 25 to 64
for n=1:(n_ret-n_l);
    for j=1:lea;
        for i=1:lez;
            for u=1:lee;
                C(j,i,u,n_ret-n)=R*a(j)+exp(z(i)+f(n_ret-n)+eps(u))-a(Ind(j,i,u,n_ret-n));
                Ap(j,i,u,n_ret-n)=a(Ind(j,i,u,n_ret-n));
            end
        end
    end
end

figure(1)
plot(a,Ap(:,1,1,25),'linewidth',2)
hold on
plot(a,Ap(:,1,1,40),'linewidth',2)
hold on
plot(a,Ap(:,1,1,60),'linewidth',2)
xlabel('a')
ylabel('Ap')
legend('Age=25','Age=40','Age=60','location','southeast')
title('Policy function at different ages')
saveas(figure(1),'Figure1.jpg')

%part b
CGamma=[0.6 0.95 1;0.1 0.9 1;0.05 0.4 1];
za=zeros(n_h,1);
za(n_l)=z(1);
for t=1:length(za);
    dioporco(t)=rand;
end
for t=(n_l+1):n_ret;
    if za(t-1)==z(1);
        if dioporco(t)<=CGamma(1,1);
            za(t)=z(1);
        elseif dioporco(t)>CGamma(1,1) & dioporco(t)<=CGamma(1,2);
            za(t)=z(2);
        else za(t)=z(3);
        end
    end
    if za(t-1)==z(2);
        if dioporco(t)<=CGamma(2,1);
            za(t)=z(1);
        elseif dioporco(t)>CGamma(2,1) & dioporco(t)<=CGamma(2,2);
            za(t)=z(2);
        else za(t)=z(3);
        end
    end
    if za(t-1)==z(3);
        if dioporco(t)<=CGamma(3,1);
            za(t)=z(1);
        elseif dioporco(t)>CGamma(3,1) & dioporco(t)<=CGamma(3,2);
            za(t)=z(2);
        else za(t)=z(3);
        end
    end
end
for t=n_ret+1:n_h;
    za(t)=za(n_ret);
end
epsa=zeros(n_h,1);
for t=1:length(epsa);
    porcamadonna(t)=rand;
end
for t=1:n_h;
    if porcamadonna(t)<=0.5;
    epsa(t)=eps(1);
    else epsa(t)=eps(2);
    end
end
ysim=zeros(n_h,1);
for t=n_l:n_ret;
    ysim(t)=exp(za(t)+f(t)+epsa(t));
end
for t=(n_ret+1):n_h;
    ysim(t)=0.7*exp(za(n_ret)+f(n_ret));
end
age=[n_l:n_h]';

Csim=zeros(n_h,1);
Apsim=zeros(n_h,1);
for v=1:length(za);
    if za(v)==z(1);
        indy(v)=1;
    elseif za(v)==z(2);
        indy(v)=2;
    else indy(v)=3;
    end
end
for q=1:length(epsa);
    if epsa(q)==eps(1);
        indeps(q)=1;
    else indeps(q)=2;
    end
end

Apsim(24)=a(1);
for k=n_l:n_h;
    Csim(k)=R*Apsim(k-1)+ysim(k)-Ap(find(a==Apsim(k-1)),indy(k),indeps(k),k);
    Apsim(k)=Ap(find(a==Apsim(k-1)),indy(k),indeps(k),k);
end
Asim=zeros(n_h,1);
for t=1:(n_h-1);
    Asim(n_h-(t-1))=Apsim(n_h-t);
end


figure(2)
plot(age,ysim(n_l:end),'linewidth',2)
hold on
plot(age,Csim(n_l:end),'linewidth',2)
hold on
plot(age,Asim(n_l:end),'linewidth',2)
set(gca, 'XLim', [age(1) age(end)])
set(gca, 'XTick', [age(1):5:age(end)])
xlabel('Age')
ylabel('Levels')
legend('Income','Consumption','Assets','location','northwest')
saveas(figure(2),'Figure2.jpg')

%part c
NN=10000;
zas=zeros(n_h,NN);
zas(n_l,:)=z(1);
for nn=1:NN;
    for t=1:size(zas,1);
        dioporcos(t,nn)=rand;
    end
end
for nn=1:NN;
for t=(n_l+1):n_ret;
    if zas(t-1,nn)==z(1);
        if dioporcos(t,nn)<=CGamma(1,1);
            zas(t,nn)=z(1);
        elseif dioporcos(t,nn)>CGamma(1,1) & dioporcos(t,nn)<=CGamma(1,2);
            zas(t,nn)=z(2);
        else zas(t,nn)=z(3);
        end
    end
    if zas(t-1,nn)==z(2);
        if dioporcos(t,nn)<=CGamma(2,1);
            zas(t,nn)=z(1);
        elseif dioporcos(t,nn)>CGamma(2,1) & dioporcos(t,nn)<=CGamma(2,2);
            zas(t,nn)=z(2);
        else zas(t,nn)=z(3);
        end
    end
    if zas(t-1,nn)==z(3);
        if dioporcos(t,nn)<=CGamma(3,1);
            zas(t,nn)=z(1);
        elseif dioporcos(t,nn)>CGamma(3,1) & dioporcos(t,nn)<=CGamma(3,2);
            zas(t,nn)=z(2);
        else zas(t,nn)=z(3);
        end
    end
end
end
for nn=1:NN;
for t=(n_ret+1):n_h;
    zas(t,nn)=zas(n_ret,nn);
end
end

epsas=zeros(n_h,NN);
for nn=1:NN;
for t=1:size(epsas,1);
    porcamadonnas(t,nn)=rand;
end
end

for nn=1:NN;
for t=1:n_h;
    if porcamadonnas(t,nn)<=0.5;
    epsas(t,nn)=eps(1);
    else epsas(t,nn)=eps(2);
    end
end
end
ysims=zeros(n_h,NN);
for nn=1:NN;
for t=n_l:n_ret;
    ysims(t,nn)=exp(zas(t,nn)+f(t)+epsas(t,nn));
end
for t=(n_ret+1):n_h;
    ysims(t,nn)=0.7*exp(zas(n_ret,nn)+f(n_ret));
end
end

Csims=zeros(n_h,NN);
Apsims=zeros(n_h,NN);
for nn=1:NN;
for v=1:size(zas,1);
    if zas(v,nn)==z(1);
        indys(v,nn)=1;
    elseif zas(v,nn)==z(2);
        indys(v,nn)=2;
    else indys(v,nn)=3;
    end
end
for q=1:size(epsas,1);
    if epsas(q,nn)==eps(1);
        indepss(q,nn)=1;
    else indepss(q,nn)=2;
    end
end
end

Apsims(24,:)=a(1);
Asims=zeros(n_h,NN);
for nn=1:NN;
for k=n_l:n_h;
    Csims(k,nn)=R*Apsims(k-1,nn)+ysims(k,nn)-Ap(find(a==Apsims(k-1,nn)),indys(k,nn),indepss(k,nn),k);
    Apsims(k,nn)=Ap(find(a==Apsims(k-1,nn)),indys(k,nn),indepss(k,nn),k);
end
for t=1:(n_h-1);
    Asims(n_h-(t-1),nn)=Apsims(n_h-t,nn);
end
end

for t=1:n_h;
    AveCsims(t)=(1/NN)*ones(1,NN)*(Csims(t,:)');
    Aveysims(t)=(1/NN)*ones(1,NN)*(ysims(t,:)');
    AveAsims(t)=(1/NN)*ones(1,NN)*(Asims(t,:)');
end
for nn=1:NN;
    for t=1:n_h;
        loCsims(t,nn)=log(Csims(t,nn));
        loysims(t,nn)=log(ysims(t,nn));
    end
end
for t=1:n_h;
    lo1Csims(t)=(1/NN)*ones(1,NN)*loCsims(t,:)';
    lo1ysims(t)=(1/NN)*ones(1,NN)*loysims(t,:)';
    lo2Csims(t)=(1/NN)*ones(1,NN)*(loCsims(t,:).^(2))';
    lo2ysims(t)=(1/NN)*ones(1,NN)*(loysims(t,:).^(2))';
    VarloCsims(t)=lo2Csims(t)-lo1Csims(t);
    Varloysims(t)=lo2ysims(t)-lo1ysims(t);
end


figure(3)
plot(age,Aveysims(n_l:end),'linewidth',2)
hold on
plot(age,AveCsims(n_l:end),'linewidth',2)
hold on
plot(age,AveAsims(n_l:end),'linewidth',2)
set(gca, 'XLim', [age(1) age(end)])
set(gca, 'XTick', [age(1):5:age(end)])
xlabel('Age')
ylabel('Levels')
legend('Income','Consumption','Assets','location','northwest')
saveas(figure(3),'Figure3.jpg')

figure(4)
plot(age,Varloysims(n_l:end),'linewidth',2)
hold on
plot(age,VarloCsims(n_l:end),'linewidth',2)
set(gca, 'XLim', [age(1) age(end)])
set(gca, 'XTick', [age(1):5:age(end)])
xlabel('Age')
ylabel('Logs')
legend('Income','Consumption','location','southeast')
saveas(figure(4),'Figure4.jpg')

%insurance coefficients for part e
deltaysims=zeros(n_h,1);
deltaCsims=zeros(n_h,1);
deltaysims(n_l)=Aveysims(n_l);
deltaCsims(n_l)=AveCsims(n_l);
for t=(n_l+1):n_h;
    deltaysims(t)=Aveysims(t)-Aveysims(t-1);
    deltaCsims(t)=AveCsims(t)-AveCsims(t-1);
end
deltaysims=deltaysims(n_l:end);
deltaCsims=deltaCsims(n_l:end);
Edeltaysims=(1/(n_h-n_l+1))*ones(1,(n_h-n_l+1))*deltaysims;
EdeltaCsims=(1/(n_h-n_l+1))*ones(1,(n_h-n_l+1))*deltaCsims;
E2deltaysims=(1/(n_h-n_l+1))*ones(1,(n_h-n_l+1))*(deltaysims.^(2));
Vardeltaysims=E2deltaysims-Edeltaysims;
for t=1:(n_h-n_l+1);
    deltayCsims(t)=deltaysims(t)*deltaCsims(t);
end
deltayCsims=deltayCsims';
CovdeltaysimsdeltaCsims=(1/(n_h-n_l+1))*ones(1,(n_h-n_l+1))*deltayCsims-(Edeltaysims*EdeltaCsims);
insurance_coeff=1-(CovdeltaysimsdeltaCsims/Vardeltaysims);

Aveysims0=Aveysims;
AveCsims0=AveCsims;
AveAsims0=AveAsims;
Valoysims0=Varloysims;
VarloCsims0=VarloCsims;
insurance_coeff0=insurance_coeff;

%% part d
% clear a_low V V1m C Ap Ind Apsims Asims Csims Aveysims AveCsims AveAsims Varloysims VarloCsims
clearvars -except r NN z eps age Aveysims0 AveCsims0 AveAsims0 Varloysims0 VarloCsims0 insurance_coeff0 zas epsas
n_l=25;
n_h=80;
n_ret=65;
sigma=2;
beta=0.92;
R=1+r;
lez=length(z);
Gamma=[0.6 0.35 0.05;0.1 0.8 0.1;0.05 0.35 0.6];
for n=n_l:n_h;
    f(n)=log(-1.0135+0.1086*n-0.001122*n^(2));
end
lee=length(eps);
pr_eps=[0.5 0.5];
stat_z=[1/3 1/3 1/3]*Gamma^(1000);
a_low=-50;
a_high=5000;
lea=460;
a=[linspace(a_low,a_high,lea-1) 0];
a=sort(a);
V=zeros(lea,lez,lee,n_h);
Ind=zeros(lea,lez,lee,n_h);
%value function at the age of 80
for u=1:lee;
    for j=1:lea;
        for i=1:lez;
            Rz(i)=0.7*exp(z(i)+f(n_ret));
            V(j,i,u,n_h)=(1/(1-sigma))*(R*a(j)+Rz(i))^(1-sigma);
            if R*a(j)+Rz(i)<0;
               V(j,i,u,n_h)=-999999-99*abs(R*a(j)+Rz(i));
            end
        end
    end
end

%value function from 66 to 79
for u=1:lee;
    for n=1:(n_h-n_ret-1);
        for j=1:lea;
            for i=1:lez;
                for jj=1:lea;
                    V1m(j,jj,i,u,n_h-n)=(1/(1-sigma))*(R*a(j)+Rz(i)-a(jj))^(1-sigma)+beta*V(jj,i,u,n_h-(n-1));
                    if R*a(j)+Rz(i)-a(jj)<0;
                       V1m(j,jj,i,u,n_h-n)=-999999-99*abs(R*a(j)+Rz(i)-a(jj));
                    end
                end
                [V(j,i,u,n_h-n),Ind(j,i,u,n_h-n)]=max(V1m(j,:,i,u,n_h-n));
            end
        end
    end
end
%value function at 65
for u=1:lee;
    for j=1:lea;
        for i=1:lez;
            for jj=1:lea;
                V1m(j,jj,i,u,n_ret)=(1/(1-sigma))*(R*a(j)+exp(z(i)+f(n_ret)+eps(u))-a(jj))^(1-sigma)+beta*V(jj,i,u,n_ret+1);
                if R*a(j)+exp(z(i)+f(n_ret)+eps(u))-a(jj)<0;
                   V1m(j,jj,i,u,n_ret)=-999999-99*abs(R*a(j)+exp(z(i)+f(n_ret)+eps(u))-a(jj));
                end
            end
            [V(j,i,u,n_ret),Ind(j,i,u,n_ret)]=max(V1m(j,:,i,u,n_ret));
        end
    end
end
%value function from 25 to 64
for u=1:lee;
    for n=1:(n_ret-n_l);
        for j=1:lea;
            for i=1:lez;
                for jj=1:lea;
                    V1m(j,jj,i,u,n_ret-n)=(1/(1-sigma))*(R*a(j)+exp(z(i)+f(n_ret-n)+eps(u))-a(jj))^(1-sigma)+beta*(Gamma(i,1)*(0.5*V(jj,1,1,n_ret-(n-1))+0.5*V(jj,1,2,n_ret-(n-1)))+Gamma(i,2)*(0.5*V(jj,2,1,n_ret-(n-1))+0.5*V(jj,2,2,n_ret-(n-1)))+Gamma(i,3)*(0.5*V(jj,3,1,n_ret-(n-1))+0.5*V(jj,3,2,n_ret-(n-1))));
                    if R*a(j)+exp(z(i)+f(n_ret-n)+eps(u))-a(jj)<0;
                       V1m(j,jj,i,u,n_ret-n)=-999999-99*abs(R*a(j)+exp(z(i)+f(n_ret-n)+eps(u))-a(jj));
                    end
                end
                [V(j,i,u,n_ret-n),Ind(j,i,u,n_ret-n)]=max(V1m(j,:,i,u,n_ret-n));
            end
        end
    end
end

%consumption and policy function
C=zeros(lea,lez,lee,n_h);
Ap=zeros(lea,lez,lee,n_h);
%at age 80
for j=1:lea;
    for i=1:lez;
        C(j,i,:,n_h)=R*a(j)+Rz(i);
        Ap(j,i,:,n_h)=0;
    end
end
%at age 66 to 79
for n=1:(n_h-n_ret-1);
    for j=1:lea;
        for i=1:lez;
            C(j,i,:,n_h-n)=R*a(j)+Rz(i)-a(Ind(j,i,:,n_h-n));
            Ap(j,i,:,n_h-n)=a(Ind(j,i,:,n_h-n));
        end
    end
end
%at age 65
for j=1:lea;
    for i=1:lez;
        for u=1:lee;
            C(j,i,u,n_ret)=R*a(j)+exp(z(i)+f(n_ret)+eps(u))-a(Ind(j,i,u,n_ret));
            Ap(j,i,u,n_ret)=a(Ind(j,i,u,n_ret));
        end
    end
end
%at age 25 to 64
for n=1:(n_ret-n_l);
    for j=1:lea;
        for i=1:lez;
            for u=1:lee;
                C(j,i,u,n_ret-n)=R*a(j)+exp(z(i)+f(n_ret-n)+eps(u))-a(Ind(j,i,u,n_ret-n));
                Ap(j,i,u,n_ret-n)=a(Ind(j,i,u,n_ret-n));
            end
        end
    end
end

ysims=zeros(n_h,NN);
for nn=1:NN;
for t=n_l:n_ret;
    ysims(t,nn)=exp(zas(t,nn)+f(t)+epsas(t,nn));
end
for t=(n_ret+1):n_h;
    ysims(t,nn)=0.7*exp(zas(n_ret,nn)+f(n_ret));
end
end

Csims=zeros(n_h,NN);
Apsims=zeros(n_h,NN);
for nn=1:NN;
for v=1:size(zas,1);
    if zas(v,nn)==z(1);
        indys(v,nn)=1;
    elseif zas(v,nn)==z(2);
        indys(v,nn)=2;
    else indys(v,nn)=3;
    end
end
for q=1:size(epsas,1);
    if epsas(q,nn)==eps(1);
        indepss(q,nn)=1;
    else indepss(q,nn)=2;
    end
end
end

Apsims(24,:)=0;
Asims=zeros(n_h,NN);
for nn=1:NN;
for k=n_l:n_h;
    Csims(k,nn)=R*Apsims(k-1,nn)+ysims(k,nn)-Ap(find(a==Apsims(k-1,nn)),indys(k,nn),indepss(k,nn),k);
    Apsims(k,nn)=Ap(find(a==Apsims(k-1,nn)),indys(k,nn),indepss(k,nn),k);
end
for t=1:(n_h-1);
    Asims(n_h-(t-1),nn)=Apsims(n_h-t,nn);
end
end

for t=1:n_h;
    AveCsims(t)=(1/NN)*ones(1,NN)*(Csims(t,:)');
    Aveysims(t)=(1/NN)*ones(1,NN)*(ysims(t,:)');
    AveAsims(t)=(1/NN)*ones(1,NN)*(Asims(t,:)');
end
for nn=1:NN;
    for t=1:n_h;
        loCsims(t,nn)=log(Csims(t,nn));
        loysims(t,nn)=log(ysims(t,nn));
    end
end
for t=1:n_h;
    lo1Csims(t)=(1/NN)*ones(1,NN)*loCsims(t,:)';
    lo1ysims(t)=(1/NN)*ones(1,NN)*loysims(t,:)';
    lo2Csims(t)=(1/NN)*ones(1,NN)*(loCsims(t,:).^(2))';
    lo2ysims(t)=(1/NN)*ones(1,NN)*(loysims(t,:).^(2))';
    VarloCsims(t)=lo2Csims(t)-lo1Csims(t);
    Varloysims(t)=lo2ysims(t)-lo1ysims(t);
end

%insurance coefficients for part e
deltaysims=zeros(n_h,1);
deltaCsims=zeros(n_h,1);
deltaysims(n_l)=Aveysims(n_l);
deltaCsims(n_l)=AveCsims(n_l);
for t=(n_l+1):n_h;
    deltaysims(t)=Aveysims(t)-Aveysims(t-1);
    deltaCsims(t)=AveCsims(t)-AveCsims(t-1);
end
deltaysims=deltaysims(n_l:end);
deltaCsims=deltaCsims(n_l:end);
Edeltaysims=(1/(n_h-n_l+1))*ones(1,(n_h-n_l+1))*deltaysims;
EdeltaCsims=(1/(n_h-n_l+1))*ones(1,(n_h-n_l+1))*deltaCsims;
E2deltaysims=(1/(n_h-n_l+1))*ones(1,(n_h-n_l+1))*(deltaysims.^(2));
Vardeltaysims=E2deltaysims-Edeltaysims;
for t=1:(n_h-n_l+1);
    deltayCsims(t)=deltaysims(t)*deltaCsims(t);
end
deltayCsims=deltayCsims';
CovdeltaysimsdeltaCsims=(1/(n_h-n_l+1))*ones(1,(n_h-n_l+1))*deltayCsims-(Edeltaysims*EdeltaCsims);
insurance_coeff=1-(CovdeltaysimsdeltaCsims/Vardeltaysims);

figure(5)
subplot(1,2,1)
plot(age,Aveysims0(n_l:end),'linewidth',1,'color','blue')
hold on
plot(age,AveCsims0(n_l:end),'linewidth',2,'color','green','linestyle',':')
hold on
plot(age,AveAsims0(n_l:end),'linewidth',2,'color','red','linestyle',':')
hold on
plot(age,AveCsims(n_l:end),'linewidth',1,'linestyle','-','color','green')
hold on
plot(age,AveAsims(n_l:end),'linewidth',1,'linestyle','-','color','red')
set(gca, 'XLim', [age(1) age(end)])
set(gca, 'XTick', [age(1):5:age(end)])
xlabel('Age')
ylabel('Levels')
legend('Income','Consumption_{0}','Assets_{0}','Consumption_{1}','Assets_{1}','location','northwest')
title('Means')

subplot(1,2,2)
plot(age,Varloysims(n_l:end),'linewidth',1,'color','blue')
hold on
plot(age,VarloCsims0(n_l:end),'linewidth',2,'color','red','linestyle',':')
hold on
plot(age,VarloCsims(n_l:end),'linewidth',1,'color','red','linestyle','-')
set(gca, 'XLim', [age(1) age(end)])
set(gca, 'XTick', [age(1):5:age(end)])
xlabel('Age')
ylabel('Logs')
legend('Income','Consumption_{0}','Consumption_{1}','location','southeast')
title('Variances')
saveas(figure(5),'Figure5.jpg')

%%
% To compute the insurance coefficients for permanent set:
% eps=[log(0.8) log(0.8)]; or log(1.2)
% To compute the insurance coefficients for transitory set:
% z=[log(15) log(15) log(15)]; or log(80) or log(500)

figure(6)
subplot(2,3,1)
histogram(Asims(26,:))
xlabel('Assets')
ylabel('Frequency')
title('Age=26')
subplot(2,3,2)
histogram(Asims(35,:))
xlabel('Assets')
ylabel('Frequency')
title('Age=35')
subplot(2,3,3)
histogram(Asims(45,:))
xlabel('Assets')
ylabel('Frequency')
title('Age=45')
subplot(2,3,4)
histogram(Asims(60,:))
xlabel('Assets')
ylabel('Frequency')
title('Age=60')
subplot(2,3,5)
histogram(Asims(65,:))
xlabel('Assets')
ylabel('Frequency')
title('Age=65')
subplot(2,3,6)
histogram(Asims(80,:))
xlabel('Assets')
ylabel('Frequency')
title('Age=80')

figure(7)
subplot(2,3,1)
histogram(zas(26,:))
xlabel('z')
ylabel('Frequency')
title('Age=26')
subplot(2,3,2)
histogram(zas(35,:))
xlabel('z')
ylabel('Frequency')
title('Age=35')
subplot(2,3,3)
histogram(zas(45,:))
xlabel('z')
ylabel('Frequency')
title('Age=45')
subplot(2,3,4)
histogram(zas(60,:))
xlabel('z')
ylabel('Frequency')
title('Age=60')
subplot(2,3,5)
histogram(zas(65,:))
xlabel('z')
ylabel('Frequency')
title('Age=65')
subplot(2,3,6)
histogram(zas(80,:))
xlabel('z')
ylabel('Frequency')
title('Age=80')

figure(8)
subplot(2,3,1)
histogram(epsas(26,:))
xlabel('epsilon')
ylabel('Frequency')
title('Age=26')
subplot(2,3,2)
histogram(epsas(35,:))
xlabel('epsilon')
ylabel('Frequency')
title('Age=35')
subplot(2,3,3)
histogram(epsas(45,:))
xlabel('epsilon')
ylabel('Frequency')
title('Age=45')
subplot(2,3,4)
histogram(epsas(60,:))
xlabel('epsilon')
ylabel('Frequency')
title('Age=60')
subplot(2,3,5)
histogram(epsas(65,:))
xlabel('epsilon')
ylabel('Frequency')
title('Age=65')
subplot(2,3,6)
histogram(epsas(80,:))
xlabel('epsilon')
ylabel('Frequency')
title('Age=80')

AAveAsims=(1/n_h)*ones(1,n_h)*AveAsims' %aggregate stock of assets;