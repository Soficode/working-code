
function [X, In,We] =  Symmetric (SS)

T=1000;
N= 20;

tau_D= 0.20;
tau_F= 1.5;
tau_m=.01;
U= zeros(N,1) + 0.2;
dt=.001;

meanw = 0; 
variancew = 4;
d = 0.10;
W  = sprandn (N,N,d)*(variancew^1/2) + meanw;
We = W - tril(W,-1) + tril(W,1)';

In=[10*ones(1,100),5*zeros(1,T-1)];
%In=1+[50*zeros(1,(50)),50*ones(1,150),50*zeros(1,T-200)];
epsilon = 0.005;
Sigma = randn(T,1)*epsilon^1/2; 
% 
X=zeros(3*N,T);
 
% r= zeros(N,1) +40;
% u= zeros(N,1);
% x= zeros(N,1) +1;
% % 
r= zeros(N,1) + SS(1:N);
u= zeros(N,1) + SS (N+1:2*N);
x= zeros(N,1)+ SS(2*N + 1: 3*N);

for t=1:T
    r = r + dt/tau_m*(-r + We*diag(u.*x)*r + In(t));% + Sigma(t)/sqrt(dt));
     u = u + dt*(-u/tau_F + diag(U.*(1-u))*r);
     x = x + dt*((1-x)/tau_D - diag(u.*x)*r);
     
    X(:,t)= [r; u; x];
end

% SS = X(1:3*N,T);
% 
% re_o = X(1:N,T);

figure(1)
plot(X(1,:));
figure(2)
plot(X(N+1,:));
figure(3)
plot(X(2*N+1,:));


evalues = eig(We);    % Get the eigenvalues of J

   figure(4)
   plot(real(evalues),imag(evalues),'r*') %   Plot real and imaginary parts
   xlabel('Real')
   ylabel('Imaginary')

end


   
 