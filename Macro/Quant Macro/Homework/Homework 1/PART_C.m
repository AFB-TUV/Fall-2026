K(:,1)=[linspace(0.75,36.096,201)]';
m=1;
tol=10^-3
tic
while m<10000
    m=m+1;
    K(:,m)=partc(K(:,m-1));
    if abs(K(:,m)-K(:,m-1))<tol
        break
    end
end
toc

figure(1);
plot(K(:,1));
title('Extende path: Capital')
hold on
plot(K(:,2));
hold on
plot(K(:,3));
hold on
plot(K(:,6));
hold on
plot(K(:,11));
hold on
plot(K(:,200)); %the last iteration is the optimal one
legend('K0','K1','K2','K5','K10')

%Doing the same for l since we know the path of K
theta=0.3;
sigma=1.5;
delta=0.1;
beta=0.95;
A=(1/(theta*beta*((1-theta)^((1-theta)/(1+theta))))^((1+theta)/2));
L=zeros(201,500);
for g=1:500
L(:,g)=((A*(K(:,g).^(theta))).^(1/(1-theta)));
end

figure(2);
plot(L(:,1));
title('Extended path: Labor')
hold on
plot(L(:,2));
hold on
plot(L(:,3));
hold on
plot(L(:,6));
hold on
plot(L(:,11));
hold on
plot(L(:,200)); %the last iteration is the optimal one
legend('L0','L1','L2','L5','L10')

%Part d. Plot output, employment, investment, consumption, and the real
%interest rate from t=0 to t=200.
%Output
Y=zeros(201,500);
for g=1:500
Y(:,g)=(A.*(K(:,g).^(theta)).*(L(:,g).^(1-theta)))+K(:,g).*(1-delta);
end

figure(3);
plot(Y);
title('Extended path: Output')


C=zeros(201,499);
for g=1:499
C(:,g)=Y(:,g)-K(:,g+1);
end

figure(4);
plot(C);
title('Extended path: Consumption')

I=zeros(201,499);
for g=1:499
I(:,g)=Y(:,g)-C(:,g);
end

figure(5);
plot(I);
title('Extended path: Investment')


