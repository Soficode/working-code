T=3000;
dt=.001;
tau_D=0.2;
tau_F=1.5;
tau_m=.01;
U=.2;
We=-10;

In=1+[zeros(1,(T-200)),50*ones(1,50),zeros(1,150)];
X=zeros(3,T);
r=0; x=0; u=0;
for t=1:T
    r = r + dt/tau_m*(-r + We*x*u*r + In(t));
    x = x + dt*((1-x)/tau_D - u*x*r);
    u = u + dt*((U-u)/tau_F + U*(1-u)*r);

    X(:,t)=[r; x; u];
end

figure (1)
plot(X(1,:));



evalues = eig(We);    % Get the eigenvalues of J

   figure(2)
   plot(real(evalues),imag(evalues),'r*') %   Plot real and imaginary parts
   xlabel('Real')
   ylabel('Imaginary')

