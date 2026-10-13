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
