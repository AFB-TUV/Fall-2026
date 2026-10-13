% clear all
% clc
% 
% rho=0.75;
% sigma_eps=0.2;
% mu_eps=0;
% zbar=exp(1.39);
% mu_z=(1/(1-rho))*mu_eps;
% sigma2_eps=(sigma_eps)^(2);
% sigma2_z=sigma2_eps/(1-rho^(2));
% lambda=1;
% z(1)=mu_z-lambda*sigma_eps/(sqrt(1-rho^(2)));
% m=5;
% step=-2*z(1)/(m-1);
% for i=1:m;
%     z(i)=z(1)+(i-1)*step;
% end
% 
% for i=1:m;
%     Gamma(i,1)=normcdf(((z(1)-rho*z(i))/sigma_eps)+step/(2*sigma_eps));
%     for j=2:(m-1);
%         Gamma(i,j)=normcdf(((z(j)-rho*z(i))/sigma_eps)+step/(2*sigma_eps))-normcdf(((z(j)-rho*z(i))/sigma_eps)-step/(2*sigma_eps));
%     end
%     Gamma(i,m)=1-ones(1,(m-1))*Gamma(i,1:(m-1))';
% end
% 
% M0=(1/m)*ones(1,m);
% M1=M0*Gamma;
% Mstar=M0*Gamma^(1000);
% alpha=2/3;
% r=0.2;
% delta=0.1;
% cf=20;
% ce=40;  %substitute with 30 and 50 to see the behavior of p (1.1931 1.2617 1.3141)
% theta=100;
% p=0.9;
% z=z+log(zbar);
% z=exp(z)';
% lez=length(z);
% nmin=0;
% nmax=(alpha*p*z(end))^(1/(1-alpha));
% len=50;
% n=linspace(nmin,nmax,len)';
% V0=zeros(lez,len);
% tol=10^(-6);
% diff=tol+1000;
% its=0;
% while diff>tol;
%     for i=1:lez;
%         for j=1:len;
%             for jj=1:len;
%                 Vp(i,j,jj)=p*z(i)*n(jj)^(alpha)-n(jj)-delta*max(n(j)-n(jj),0)-cf+(1/(1+r))*Gamma(i,:)*max(V0(:,jj),-delta*n(jj));
%             end
%             [V(i,j),Ind(i,j)]=max(Vp(i,j,:));
%         end
%     end
%     diff=norm(V-V0);
%     its=its+1;
%     V0=V;
% end
% N=n(Ind);  %policy function for employment
% X=ones(lez,len);
% for i=1:lez;
%     for j=1:len;
%         Vaux(i,j)=max(V(i,Ind(i,j)),-delta*n(Ind(i,j)));
%         if abs(Vaux(i,j)-V(i,j))<tol;
%             X(i,j)=0;
%         end
%     end
% end
% 
% figure(1)
% subplot(1,3,1)
% plot(n,V(1,:),'linewidth',2)
% hold on
% plot(n,V(3,:),'linewidth',2)
% hold on
% plot(n,V(5,:),'linewidth',2)
% legend('Low','Middle','High')
% xlabel('n')
% ylabel('V')
% xlim([n(1) n(end)])
% title('Value function')
% 
% subplot(1,3,2)
% plot(n,N(1,:),'linewidth',2)
% hold on
% plot(n,N(3,:),'linewidth',2)
% hold on
% plot(n,N(5,:),'linewidth',2)
% legend('Low','Middle','High')
% xlabel('n')
% ylabel('N')
% xlim([n(1) n(end)])
% title('Policy function N')
% 
% subplot(1,3,3)
% plot(n,X(1,:),'linewidth',2)
% hold on
% plot(n,X(3,:),'linewidth',2)
% hold on
% plot(n,X(5,:),'linewidth',2)
% legend('Low','Middle','High')
% xlabel('n')
% ylabel('X')
% xlim([n(1) n(end)])
% title('Policy function X')
% 
% G=1/lez*ones(lez,1)';
% for i=1:lez;
%     tre(i)=max(V(i,1),0);
% end
% Ve=-ce+(1/(1+r))*G*tre';
% %%
% clear all
% clc
% 
% rho=0.75;
% sigma_eps=0.2;
% mu_eps=0;
% zbar=exp(1.39);
% mu_z=(1/(1-rho))*mu_eps;
% sigma2_eps=(sigma_eps)^(2);
% sigma2_z=sigma2_eps/(1-rho^(2));
% lambda=1;
% z(1)=mu_z-lambda*sigma_eps/(sqrt(1-rho^(2)));
% m=5;
% step=-2*z(1)/(m-1);
% for i=1:m;
%     z(i)=z(1)+(i-1)*step;
% end
% 
% for i=1:m;
%     Gamma(i,1)=normcdf(((z(1)-rho*z(i))/sigma_eps)+step/(2*sigma_eps));
%     for j=2:(m-1);
%         Gamma(i,j)=normcdf(((z(j)-rho*z(i))/sigma_eps)+step/(2*sigma_eps))-normcdf(((z(j)-rho*z(i))/sigma_eps)-step/(2*sigma_eps));
%     end
%     Gamma(i,m)=1-ones(1,(m-1))*Gamma(i,1:(m-1))';
% end
% 
% M0=(1/m)*ones(1,m);
% M1=M0*Gamma;
% Mstar=M0*Gamma^(1000);
% alpha=2/3;
% r=0.2;
% delta=0.1;
% cf=20;
% ce=40;  %substitute with 30 and 50 to see the behavior of p (1.1931 1.2617 1.3141)
% theta=100;
% p=0.9;
% p=1.2617; %price for part c, comment this for getting part a
% z=z+log(zbar);
% z=exp(z)';
% lez=length(z);
% nmin=0;
% nmax=(alpha*p*z(end))^(1/(1-alpha));
% len=50;
% n=linspace(nmin,nmax,len)';
% tic
% % pp=linspace(1.16,1.34,50);
% pp=linspace(1.27,1.37,20);
% lep=length(pp);
% for u=1:lep;
% p=pp(u);
% V0=zeros(lez,len);
% V=zeros(lez,len);
% Vp=zeros(lez,len,len);
% Ind=zeros(lez,len);
% N=zeros(lez,len);
% X=zeros(lez,len);
% tol=10^(-6);
% diff=tol+1000;
% its=0;
% while diff>tol;
%     for i=1:lez;
%         for j=1:len;
%             for jj=1:len;
%                 Vp(i,j,jj)=p*z(i)*n(jj)^(alpha)-n(jj)-delta*max(n(j)-n(jj),0)-cf+(1/(1+r))*Gamma(i,:)*max(V0(:,jj),-delta*n(jj));
%             end
%             [V(i,j),Ind(i,j)]=max(Vp(i,j,:));
%         end
%     end
%     diff=norm(V-V0);
%     its=its+1;
%     V0=V;
% end
% N=n(Ind);  %policy function for employment
% X=ones(lez,len);
% for i=1:lez;
%     for j=1:len;
%         Vaux(i,j)=max(V(i,Ind(i,j)),-delta*n(Ind(i,j)));
%         if abs(Vaux(i,j)-V(i,j))<tol;
%             X(i,j)=0;
%         end
%     end
% end
% G=1/lez*ones(lez,1)';
% tre=0;
% Ve=0;
% for i=1:lez;
%     tre(i)=max(V(i,1),0);
% end
% Ve=-ce+(1/(1+r))*G*tre';
% F(u)=Ve;
% end
% disp(ce)
% pstar=interp1(F,pp,0);
% disp(pstar)%price s.t. the value of entry is 0
% toc
% 
% figure(2)
% plot(pp,F,'linewidth',2)
% grid on
% xlabel('p')
% ylabel('F(p)')
% xlim([pp(1) pp(end)])
% title('Function F(p)')

%%
clear all
clc

rho=0.75;
sigma_eps=0.2;
mu_eps=0;
zbar=exp(1.39);
mu_z=(1/(1-rho))*mu_eps;
sigma2_eps=(sigma_eps)^(2);
sigma2_z=sigma2_eps/(1-rho^(2));
lambda=1;
z(1)=mu_z-lambda*sigma_eps/(sqrt(1-rho^(2)));
m=5;
step=-2*z(1)/(m-1);
for i=1:m;
    z(i)=z(1)+(i-1)*step;
end

for i=1:m;
    Gamma(i,1)=normcdf(((z(1)-rho*z(i))/sigma_eps)+step/(2*sigma_eps));
    for j=2:(m-1);
        Gamma(i,j)=normcdf(((z(j)-rho*z(i))/sigma_eps)+step/(2*sigma_eps))-normcdf(((z(j)-rho*z(i))/sigma_eps)-step/(2*sigma_eps));
    end
    Gamma(i,m)=1-ones(1,(m-1))*Gamma(i,1:(m-1))';
end

M0=(1/m)*ones(1,m);
M1=M0*Gamma;
Mstar=M0*Gamma^(1000);
alpha=2/3;
r=0.2;
delta=0.1;
cf=20;
ce=40;  %substitute with 30 and 50 to see the behavior of p (1.1931 1.2617 1.3141)
theta=100;
p=1.3534; %price for part c, comment this for getting part a
z=z+log(zbar);
z=exp(z)';
lez=length(z);
nmin=0;
nmax=(alpha*p*z(end))^(1/(1-alpha));
len=50;
n=linspace(nmin,nmax,len)';
V0=zeros(lez,len);
tol=10^(-6);
diff=tol+1000;
its=0;
while diff>tol;
    for i=1:lez;
        for j=1:len;
            for jj=1:len;
                Vp(i,j,jj)=p*z(i)*n(jj)^(alpha)-n(jj)-delta*max(n(j)-n(jj),0)-cf+(1/(1+r))*Gamma(i,:)*max(V0(:,jj),-delta*n(jj));
            end
            [V(i,j),Ind(i,j)]=max(Vp(i,j,:));
        end
    end
    diff=norm(V-V0);
    its=its+1;
    V0=V;
end
N=n(Ind);  %policy function for employment
X=ones(lez,len);
for i=1:lez;
    for j=1:len;
        Vaux(i,j)=max(V(i,Ind(i,j)),-delta*n(Ind(i,j)));
        if abs(Vaux(i,j)-V(i,j))<tol;
            X(i,j)=0;
        end
    end
end

G=1/lez*ones(lez,1)';
for i=1:lez;
    tre(i)=max(V(i,1),0);
end
Ve=-ce+(1/(1+r))*G*tre';

for i=1:len;
    if n(i)==0;
        Inl(i)=1;
    else Inl(i)=0;
    end
end
for j=1:lez;
    for i=1:len;
        for ii=1:len;
        if N(j,i)==n(ii);
            INnl(j,i,ii)=1;
        else INnl(j,i,ii)=0;
        end
        end
    end
end

m=0.1;
m=0.00165;
diff1=tol+1000;
its1=0;
GAMMA0=ones(lez,len);
while diff1>tol;
for j=1:lez;
    for i=1:len;
        for jj=1:lez;
            for ii=1:len;
                GAMMA111(j,i,jj,ii)=GAMMA0(jj,ii)*(1-X(jj,ii))*Gamma(jj,j)*INnl(jj,ii,i)+m*G(j)*Inl(i);
            end
        end
    end
end
GAMMA11=sum(GAMMA111,4);
GAMMA1=sum(GAMMA11,3);
diff1=norm(GAMMA1-GAMMA0);
its1=its1+1;
GAMMA0=GAMMA1;
end
D=theta/p;
for j=1:lez;
    for i=1:len;
        Qpp(j,i)=GAMMA1(j,i)*z(j)*N(j,i)^(alpha)*(1-X(j,i));
        Qwpp(j,i)=GAMMA1(j,i)*z(j)*N(j,i)^(alpha-1)*(1-X(j,i));
    end
end
Q=sum(sum(Qpp,2));
Qw=sum(sum(Qwpp,2));
% disp(Q)
% disp(D)
TotNumFirms=sum(sum(GAMMA1,2));
exit=sum(sum(GAMMA1.*X,2))/TotNumFirms;
L=sum(sum(GAMMA1.*N.*(1-X),2));

[zz nn] = meshgrid(z, n);
stem3(zz',nn',GAMMA1,'Marker','s','MarkerEdgeColor','m','MarkerFaceColor','g');

for j=1:lez;
    for i=1:len;
        TACp(j,i)=GAMMA1(j,i)*delta*max(n(i)-N(j,i),0)*(1-X(j,i));
    end
end
TAC=sum(sum(TACp,2));