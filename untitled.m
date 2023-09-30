% Solves the 1D heat equation with an explicit finite difference scheme
% 
% From Thorsten Becker's webpage
%

clear all
close all

% Physical parameters
L = 100; % Length of modeled domain [m]
Tmagma = 1200; % Temperature of magma [C]
Trock = 300; % Temperature of country rock [C]
kappa = 1e-6; % Thermal diffusivity of rock [m2/s]
W = 5; % Width of dike [m]
day = 3600*24; % # seconds per day
dt = 1.45*day; % Timestep [s]

% Numerical parameters
nx = 201; % Number of gridpoints in x-direction
nt = 100; % Number of timesteps to compute
dx = L/(nx-1) ;% Spacing of grid
x = -L/2:dx:L/2;% Grid

% Setup initial temperature profile
T = ones(size(x))*Trock;
T(find(abs(x)<=W/2)) = Tmagma;
c=kappa*dt/(dx*dx);
g=zeros(nx,nx);
g(1,1)=1;
%for i=1:201 
 %   g(1,i)=1;
%end

for i=2:201 
    for j=1:201
        if j==i-1 
            g(i,j)= c;
        end
        if j==i
            g(i,j)= 1-2*c;
        end
        if j==i+1
            g(i,j)= c;
        end
    end
end
g(nx,nx)=1;

time = 0; 

for n=1:nt % Timestep loop
    % Compute new temperature
   % Tnew = zeros(1,nx);
  %  for i=2:nx-1
   %     Tnew(i) =  c*T(i-1)+c*T(i+1)+ (1-2*c)*T(i);
    %end
    % Set boundary conditions
    
    Tnew=g*T';
    
   Tnew(1) = T(1);
   Tnew(nx) = T(nx);
    
    % Update temperature and time
    T = Tnew';
    time = time+dt;
    (dt/(dx*dx))*kappa

    % Plot solution
    figure(1), clf
    plot(x,Tnew);
    Ymat(n,:) = n*dt*(ones(1,nx));
    Xmat(n,:) = x;
    Tmat(n,:) = Tnew;
    axis([-50 50 200 1300])
    xlabel('x [m]')
    ylabel('Temperature [^oC]')
    title(['Temperature evolution after ',num2str(time/day),' days'])
    drawnow
end
figure
surf(Xmat,Ymat,Tmat)