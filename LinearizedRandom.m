
function [J_x, Z, ue_o, xe_o, We, r, u, x] = LinearizedRandom (re_o)

T=1000;
N= 500;
tau_m= 0.01;
tau_d=0.2;
tau_f=1.5;
U= zeros(N,1) + 0.2; 
dt=.0001;
I = eye(N);

% 
meanw = 0; 
variancew = 4;
d = 0.10;
We  = sprandn(N,N,d)*(variancew^1/2)/(sqrt(N)) + meanw/N;
%We= zeros(N,N) + 1/N;

%In=[100*ones(1,1),0*ones(1,T-1)];
In=1+[50*zeros(1,(50)),10*ones(1,50),50*zeros(1,T-100)];
epsilon = 0.005;
Sigma = randn(T,1)*epsilon^1/2; 
    
    %Steady States
   
    ue_o = U.*(1+tau_f.*re_o./1+U.*tau_f.*re_o);
   
    xe_o = 1./(1+(ue_o*tau_d.*re_o));
 
    Ds_o = (ue_o.*xe_o);
    Ds_o = diag(Ds_o);
    
    Df_o = (ue_o.*re_o);
    Df_o = diag(Df_o);
   
    Dd_o = (re_o.*xe_o);
    Dd_o = diag(Dd_o);
        
     
    %Linearized System - Jacobian Matrix (effective connectivity matrix of
    %the linearized system)
  
    a1 = I + 1/tau_m.*(-I + 1/N*Ds_o);
  
   
    a2 = 1/tau_m.*(1/N*(Dd_o));

    a3 = 1/tau_m.*( 1/N*(Df_o)); 
 
    b1 = I + U.*(1-ue_o');
 
    b2 = (-1./tau_f-U*re_o');

    b3 = zeros(N,N);

    
    c1 = Ds_o;
   
   
    c2 = Dd_o ;
 
    c3 = (-(1./tau_d)-Df_o);

 
%     
%  %Effective Connectivity Matrices

J_x = [ a1 a2 a3; b1 b2 b3; c1 c2 c3];  
%J_x = zeros(3*N,3*N) + (J_x/sqrt(3*N));
J_I = [ones(1,N),zeros(1,2*N)]';

r= zeros(N,1);
u= zeros(N,1);
x= zeros(N,1);

Z=zeros(3*N,T);
z= [r;u;x];

for t=1:T
    
   z = z + dt*(-z + J_x*z + In(t));% +  Sigma(t)/sqrt(dt)));
  
    Z(:,t)= z;
end


figure (1)
plot(Z(1,:));
figure (2)
plot(Z(N+1,:));
figure (3)
plot(Z(2*N+1,:));

figure(6)
plot(Z(1:N,:)');


evalues = eig(J_x);    % Get the eigenvalues of J

   figure(4)
   plot(real(evalues),imag(evalues),'r*') %   Plot real and imaginary parts
   xlabel('Real')
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
   ylabel('Imaginary')
  
   

end