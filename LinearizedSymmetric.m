
function [J_x, Z, ue_o, xe_o] = LinearizedSymmetric (re_o)

T=1000;
N= 20;
tau_m= 0.01;
tau_d=0.2;
tau_f=1.5;
U= zeros(N,1) + 0.2; 
dt=.001;
I = eye(N);


meanw = 0; 
variancew =4;
d = 0.10;
W  = sprandn (N,N,d)*(variancew^1/2) + meanw;
We = W - tril(W,-1) + tril(W,1)';

In=[10*ones(1,100),5*zeros(1,T-1)];
epsilon = 0.005;
Sigma = randn(T,1)*epsilon^1/2; 
    
    %Steady States
   
    ue_o = U.*(1+tau_f*re_o/1+U.*re_o*tau_f);
   
    xe_o = 1/(1+(ue_o.*re_o*tau_d));
    xe_o = xe_o';
 
    Ds_o = (ue_o.*xe_o);
    Ds_o = diag(Ds_o);
    
    Df_o = (ue_o.*re_o);
    Df_o = diag(Df_o);
   
    Dd_o = (re_o.*xe_o);
    Dd_o = diag(Dd_o);
        
     
    %Linearized System - Jacobian Matrix (effective connectivity matrix of
    %the linearized system)
  
    a1 = 1/tau_m*(-I + We*Ds_o);

  
    a2 = 1/tau_m*(We*(Dd_o));

    a3 =1/tau_m*( We*(Df_o)); 
 
    b1 = 1/tau_f*(U*ue_o');
 
    b2 = 1/tau_f*(-1/tau_f-U*re_o');

    b3 = 1/tau_f*(zeros(N,N));

    c1 =1/tau_d*(Ds_o);
   
    c2 = 1/tau_d*(Dd_o );

    c3 = 1/tau_d*(-1/tau_d+diag(U)*(Ds_o));

 
%     
%  %Effective Connectivity Matrices

J_x = [ a1 a2 a3; b1 b2 b3; c1 c2 c3];  
J_x = zeros(3*N,3*N) + (J_x/sqrt(3*N));
J_I = [ones(1,N),zeros(1,2*N)]';
r= zeros(N,1)+re_o;
u= zeros(N,1);
x= zeros(N,1)+1;
Z=zeros(3*N,T);
z= [r;x;u];

for t=1:T

   z = z + dt*(-z + J_x*z) +J_I*In(t);% +  Sigma(t)/sqrt(dt)));
  
    Z(:,t)= z;
end


figure (1)
plot(Z(1,:));
figure (2)
plot(Z(N+1,:));
figure (3)
plot(Z(2*N+1,:));

Identity = eye(3*N);

evalues = eig(J_x-Identity);    % Get the eigenvalues of J

   figure(4)
   plot(real(evalues),imag(evalues),'r*') %   Plot real and imaginary parts
   xlabel('Real')
   ylabel('Imaginary')
  
   

end