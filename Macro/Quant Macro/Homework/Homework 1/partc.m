function [k]=partc(kk)
theta=0.3;
sigma=1.5;
delta=0.1;
beta=0.95;
A=(1/(theta*beta*((1-theta)^((1-theta)/(1+theta))))^((1+theta)/2));
kss=(1-beta*(1-delta))^(-((1+theta)/(1-theta))); %Capital of steady state
%lss=(1-beta*(1-delta))^(1/(1-theta))/((beta*theta*A)^(1/(1-theta)))*kss; % Labor of steady state
% kk0 = linspace(k0,kss,201);
k(1)=kk(1);
k(201)=kk(201);
options = optimset('Display','off');
for j=2:200;
    KKt=k(j-1);   
    KKt2=kk(j+1);
    k(j)=fsolve(@PARTC,KKt,options);
end
function [z]=PARTC(KKt1)
z=(((1/(beta*theta))*KKt^((2*theta)/(1+theta))+KKt*(1-theta)-KKt1-((1-theta)/(2*beta*theta))*KKt1^((2*theta)/(1+theta)))^(-sigma)-beta*(((1/(beta*theta))*KKt1^((2*theta)/(1+theta))+KKt1*(1-theta)-KKt2-((1-theta)/(2*beta*theta))*KKt1^((2*theta)/(1+theta)))^(-sigma))*((1/beta)*KKt1^(-(1-theta)/(1+theta))+1-delta));
end
end






