function [ R, phi ] = compute_fem(M)
%SOLVE_FEM Summary of this function goes here
%   M = % Number of elements
%   Detailed explanation goes here
c = 299792458;
f = 1e6; % initial frequency
w = 2*pi*f;
k0 = w/c;
lamb0 = c/f;
N = M+1; % Number of nodes
L = 5*lamb0; % Length of slab is 5x free space wavelength
x = 0:L/M:L; % discretize dielectric slab with N 
y = 0:L/(M-1):L; % discretize dielectric slab with M elements
phi = 0:pi/(2*M):pi/2;

mu_r = ones(1,length(x));
mu_r(1:end-1) = 2-0.1j;
e_r = 4+(2-0.1j)*((1-y/L).^2);
e_r = padarray(e_r, [0 1],1,'post'); % pad array with '1' for free space
l = L/M;
alpha = 1./mu_r;
R = zeros(1,length(phi));

for k = 1:numel(phi)
    % initialize f
    f = zeros(N,1);
    % initialize the K matrix
    K = zeros(N,N);
    THETA = phi(k);
    beta = -(k0^2)*(e_r - alpha.*sin(THETA).^2);
    K(1,1) = alpha(1)/l + beta(1)*l/3;
    K(N,N) = alpha(M)/l + beta(M)*l/3;
    for i = 2:N-1
        K(i,i) = alpha(i-1)/l + beta(i-1)*l/3 + alpha(i)/l + beta(i)*l/3;
        K(i+1,i) = -alpha(i)/l + beta(i)*l/6;
        K(i,i+1) = K(i+1,i);
    end
    % initialize the b vector
    b = zeros(N,1);
    b(1) = f(1)*l/2;
    b(N) = f(M)*l/2;
    for i = 2:N-1
        b(i) = f(i-1)*l/2 + f(i)*l/2;
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
    %% Compute numerical solution
    Ez = K\b;
    E_inc = E0*exp(1j*k0*L*cos(THETA));
    R(k) = (Ez(end) - E_inc)/conj(E_inc);
end

end

