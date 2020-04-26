close all; clear all; clc;
% Upper bound Problem by using CNAB2 algorithm.
% CN method for linear terms & AB2 nonlinear terms.
% Chebyshev (cheb) spectral collocation.
% we have 4 equations with 4 unknows (tau,theta,gamma, and w).
%% define the varibles
Ra= 50;   %4*pi*pi; % Ra number 998.
% ra = Ra/a so if a=0.5 then ra = 2*Ra. if a =1 we do not have a bound.
% check for the bast value for a! between (0,1).
ra = 2*Ra; 
%% time stuff:
Nt   = 3;       % save time  % 3200
t    = 0;
dt   = 2e-2;       % time step
nsave = 10;       % save energy
tmax = dt*Nt;    % t final
%% get cheb matrices on z:
Nz = 63;  % Nz graid points in phy space
[D0,D1,D2,zc]= cheb(Nz);           % check cheb.m
% Nz = Nz +1;                      % I can change Nx to Nz.
%% grid stuff:
% want the grid to go from -1 to 1 in order to apply D1, D2, etc
% domain discretization
Lz = 1;                       % Domain length
%x = dx*(0:Nx-1);             % the grid
%zc=x;                        % from cheb
z =(zc+1)/2;                  % zc is x from cheb. 
% zc = 2*x - 1;               % this goes from [-1, 1] 
%% x grid stuff: x is Lx-periodic
Nx = 32;   % number of grid points in x 
Lx = 2.0; 
dx = (Lx-0)/(Nx-1); 

%%
kvector = (2.0*pi/Lx)*linspace(1, Nx, Nx);%define a vector that holds k values

%% Set the initial condition, we have two:

%% I.C

% declare matrices - tau is x-averaged temperature - so doesn't depend on x

theta0 = zeros(Nz+1, Nx);
tau0   = zeros(Nz+1, 1);
W0     = zeros(Nz+1, Nx);
gamma0 = zeros(Nz+1, Nx);

theta1 = zeros(Nz+1, Nx);
tau1   = zeros(Nz+1, 1);
W1     = zeros(Nz+1, Nx);
gamma1 = zeros(Nz+1, Nx);

% I.C
for i=1:Nz+1
   tau0(i) = 1 - z(i);
   
   for j = 1:Nx
       theta0(i, j) = sin(j*pi*z(i));
   end
end


tau1 = tau0;
theata1 = theta0;
%%
Ainv = zeros(Nz+1, Nz+1, Nx);

for j=1:Nx
    k        = kvector(j);
    A        = (4.0*D2-k*k*D0); % the 4 comes due to cheb grid- coordinate transform from z to zc
    
    % Dirichlet BC
    A(1, 1)  = 1.0; A(1, 2:end) = 0.0;
    A(Nz+1, 1:Nz) = 0.0; A(Nz+1, Nz+1) = 1.0;
    
    Ainv(:, :, j) =  inv(A);
end
%%
for j = 1:Nx % n in the paper is j in code 
        k = kvector(j);
        W0(:, j) = (-ra * k * k)* Ainv(:, :, j) *  theta0(:, j);
        W1(:, j) = (-ra * k * k)* Ainv(:, :, j) *  theta1(:, j);
        
        gamma0(:, j) = Ainv(:, :, j) * ((-ra * 2.0 * D1 * tau0).* theta0(:, j));
        gamma1(:, j) = Ainv(:, :, j) * ((-ra * 2.0 * D1 * tau1).* theta1(:, j));
end

%% declare rhs matrices, vecs, etc
% A = zeros(Nz+1, Nz+1); overwrite the above A for each k as required
b = zeros(Nz+1, 1);

%% time loop CNAB2 
% (i, j) == ith row and jth column
% i corrosponds to zc
% j corrosponds to x
for n = 1:Nt  
    
    for j = 1:Nx % n in the paper is j in code 
        k = kvector(j);
        
        %% step # 1: update theta
        Ltheta = 2.0 * (4.0 * D2 - k * k * D0);
        Ntheta0 = -(2.0 * D1 * tau0) .* W0(:, j) - k * k * gamma0 (:, j);
        Ntheta1 = -(2.0 * D1 * tau1) .* W1(:, j) - k * k * gamma1 (:, j);
        
        % solve A * theta2 = b problem to update theta
        A = D0 - (0.5 * dt)*Ltheta;
        b = (D0 + (0.5 * dt)*Ltheta) * theta1(:, j) + 0.5 * dt*(3.0* Ntheta1 - Ntheta0);
        
        % apply BC's
        % Dirichlet BC
        A(1, 1)  = 1.0; A(1, 2:end) = 0.0;
        A(Nz+1, 1:Nz) = 0.0; A(Nz+1, Nz+1) = 1.0;
        
        b(1) = 0.0; b(Nz+1) = 0.0;
        
        theta2(:, j) = (A\b);
        
        %% step # 2: update W
        % solve  W2 = Ainv(:, :, j)*b problem to update W2
        
        b = - (ra * k * k) * theta2(:, j);

        % apply BC's
        % Dirichlet BC
        b(1) = 0.0; b(Nz+1) = 0.0;
        
        W2(:, j) = Ainv(:, :, j)*b;
        
        %% step # 3: solve for tau
        % get Wj* theta j here, sum columns later
        Wj_thetaj_0(:, j) = W1(:, j).*theta1(:, j);
        Wj_thetaj_1(:, j) = W2(:, j).*theta2(:, j);
        
    end
    %% step # 3: solve for tau cntd

    sum_Wj_thetaj0 = (sum(Wj_thetaj_0'))';
    sum_Wj_thetaj1 = (sum(Wj_thetaj_1'))';
%     sum_Wj_thetaj0 = zeros(Nz+1, 1);
%     sum_Wj_thetaj1 = zeros(Nz+1, 1);
%     for i = 1:Nz+1
%         for j = 1:Nx
%             sum_Wj_thetaj0(i) =  sum_Wj_thetaj0(i) + Wj_thetaj_0(i, j);
%             sum_Wj_thetaj1(i) =  sum_Wj_thetaj1(i) + Wj_thetaj_1(i, j);
%         end
%     end

    Ltau = 4.0 * D2;
    A = (D0 - 0.5 * dt * Ltau);   
    b = (D0 + 0.5 * dt * Ltau)*tau1 - 0.125*dt*(3.0*2.0*D1*sum_Wj_thetaj1 - 2.0*D1*sum_Wj_thetaj0);

    % apply BC's
    % Dirichlet BC
    A(1, 1)  = 1.0; A(1, 2:end) = 0.0;
    A(Nz+1, 1:Nz) = 0.0; A(Nz+1, Nz+1) = 1.0;

    b(1) = 0.0; b(Nz+1) = 1.0;

    tau2 = (A\b);

    %% solve for gamma
    for j = 1:Nx % n in the paper is j in code, solve for gamma
        % solve A * gamma2 = b problem to update gamma2

        b = -ra * (2.0 * D1 * tau2) .* theta2(:, j);
        % apply BC's
        % Dirichlet BC

        b(1) = 0.0; b(Nz+1) = 0.0;

        gamma2(:, j) = Ainv(:, :, j)*b;

    end
    %% update all variables for time-marching
    theta0 = theta1;
    theta1 = theta2;
    
    W0 = W1;
    W1 = W2;
    
    tau0 = tau1;
    tau1 = tau2;
    
    gamma0 = gamma1;
    gamma1 = gamma2;
end
%% ploting
figure(1)
plot(tau1, z, '-o');

%%%%
% After plotting, I will change Tb in my code in
% Energy Stability and replace it with tau >> plot(kvector,lambda1)  
      xlabel('tau','Interpreter','latex');
      ylabel('z','Interpreter','latex');
      %title(Check for values of tau);
%      grid on; 
%%
