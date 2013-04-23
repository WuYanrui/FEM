%% Project #2
% Authors: Blake Levy and Adedayo Lawal
% Quadratic Elements & Gauss-Legendre quadrature
clc;clear;
close all;
tic
%% set up constants from constants.m
constants

%% Setup geometry
M = 100; % Number of elements
N = 2*M+1; % Number of nodes
L = 5*lamb0; % Length of slab is 5x free space wavelength
y = L/(N-1):2*L/(N-1):L-L/(N-1); % discretize dielectric slab with M elements
x = 0:L/(N-1):L; % discretize dielectric slab into N points
phi = 0;%:pi/(2*M):pi/2;
%% Analytical solution
e_r = 4*ones(1,length(y));%+(2-0.1j)*((1-y/L).^2);
e_r = padarray(e_r, [0 1],1,'post'); % pad array with '1' for free space
mu_r = ones(1,length(e_r));
mu_r(1:end-1) = 1;
% compute kx
kx = k0*sqrt(mu_r.*e_r - sin(phi).^2);
% compute eta
eta = (kx(:,[2 3:end]) - kx(:,[1 2:end-1]))./(kx(:,[2 3:end]) + kx(:,[1 2:end-1]));
eta(:,end) = (mu_r(end-1).*kx(:,end) - mu_r(end).*kx(:,end-1))./(mu_r(end-1).*kx(:,end) + mu_r(end).*kx(:,end-1));
R = zeros(length(phi),M+1);
% initialize R at the PEC
R(:,1) = -1;
% compute R
for i = 2:M+1
    j = 2*(i-1)+1;
    R(:,i) = exp(2i*kx(:,i)*x(j)).*((eta(:,i-1) + R(:,i-1).*exp(-2i.*kx(:,i-1)*x(j))) ./ ...
             (1 + eta(:,i-1).*R(:,i-1).*exp(-2i*kx(:,i-1)*x(j))));
end
%% Numerical simulation
Rn = zeros(1,N);

for k = 1:numel(phi)
    % initialize f
    f = zeros(N,1);
    % initialize the K matrix
    K = zeros(N,N);
    alpha = 1/(2-0.1j);
    THETA = phi(k);
    le = x(3)-x(1);
    beta = -(k0^2)*(e_r - alpha.*sin(THETA).^2);    
    
    for i = 1:M
        j = 2*i-1;
        le = x(j+2)-x(j);
        Ke = compute_ke_quad(alpha, beta(i), le); 
        K(j,j) = K(j,j)+Ke(1,1);
        K(j,j+1) = Ke(1,2);
        K(j,j+2) = Ke(1,3);
        K(j+1,j) = K(j,j+1);
        K(j+1,j+1) = Ke(2,2);
        K(j+1,j+2) = Ke(2,3);
        K(j+2,j) = K(j,j+2);
        K(j+2,j+1) = K(j+1,j+2);
        K(j+2,j+2) = Ke(3,3);
    end
    
    % initialize the b vector
    b = zeros(N,1);
    for i = 1:M
        j = 2*i-1;
        le = x(j+2)-x(j);
        be = compute_be_quad(le,f(j:j+2));
        b(i) =  b(i) + be(1);
        b(i+1) = be(2);
        b(i+2) = be(3);
    end
    
    % apply boundary conditions
    K(1,1) = 1;
    K(1,2:end) = 0;
    K(2:end,1) = 0;
    E0 = 1; % 1V/m
    q = 2*1j*k0*cos(THETA)*E0*exp(1j*k0*L*cos(THETA));
    gamma = 1j*k0*cos(THETA);
    K(N,N) = K(N,N) + gamma;
    b(N) = b(N) + q;
    
    % Compute numerical solution
    Ez = K\b;
    E_inc = E0*exp(1j*k0*L*cos(THETA));
    Rn = (Ez - E_inc)/conj(E_inc);
    
end
%% Error result
ks = k0*sqrt(4);
Z_in = j*sqrt(mu0/4/eps0)*tan(k0*sqrt(4)*x);
R_exact = (Z_in - sqrt(mu0/eps0))./(Z_in + sqrt(mu0/eps0));
Ez_exact = (1+R_exact(end))*exp(-1j*sqrt(4)*k0.*x);
error = abs(Ez-transpose(Ez_exact));
semilogy(x(2:end)/lamb0,error(2:end))
toc