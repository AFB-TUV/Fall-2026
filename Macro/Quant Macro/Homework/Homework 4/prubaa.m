%Quantitative Macro
%Professor: Juan Sanchez
%Studen: Rocio Eugenia Suarez
%Homework 4

%% PART A
%%

clear all
clc

sigma=2;
rho=0.975;
beta=0.93;
R=1.025;
y1=16000;
y2=80;
y3=15;
y = [y1, y2, y3];
Gamma = [ 0.02, 0.98, 0; 0.001, 0.78, 0.219; 0, 0.036, 0.964 ];
abar = 0;
ymax = 3*y1; 
m1 = 30; 
a= linspace(abar, ymax, m1);
A0 = a(m1)*ones(m1,length(y));
%Initial guess of the iteration
a0=zeros(length(a),3);  %guess
for i=1:length(a);
a0(i,:)=a(i);
end

Diff = 100;
tol = 1e-3;
%
%
tic
while Diff > tol
    
    for i = 1:length(y)
        for h = 1:length(a)
            c1 = R*a(h) + y(i) - a(1);
            c2 = R*a(1) + y - a0(1,:);
            E1 = c1.^(-sigma) - rho*beta*R*sum( Gamma(i,:).*c2.^(-sigma) ); 

            if E1>0
                A0(h,i) = a(1);
            else 
                
                delta = -100;
                DELTA = NaN(length(a),1); 
                k = 0;
                
                while delta<0 & k<length(a)
                    k = k+1;
                    c1 = R*a(h) + y(i) - a(k);
                    c2 = R*a(k) + y - a0(k,:);
                    delta = max(0,c1).^(-sigma) - rho*beta*R*sum( Gamma(i,:).*max(0,c2).^(-sigma));
                    DELTA(k)=delta;
                end
                
                if abs( DELTA(k) )>abs( DELTA(k-1) )
                    A0(h,i) = a(k-1);
                else 
                    A0(h,i) = a(k);
                end
                
            end
        end
    end

    Diff = max(max(abs(A0-a0)));
    a0 = A0;
    
end
toc
%
aprime=A0;   
line=a;
figure(1)
plot(aprime(:,2));
legend('Policy function - m=30');
saveas(figure(1),'figure1.jpg')

           
           