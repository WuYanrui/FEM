%% Project #2
% Authors: Blake Levy and Adedayo Lawal
% Vectorized approach
clc;clear;
tic
%% set up constants from constants.m
constants

%% Setup geometry
M = 100; % Number of elements
N = M+1; % Number of nodes
L = 5*lamb0; % Length of slab is 5x free space wavelength
x = 0:L/M:L; % discretize dielectric slab with N 
y = 0:L/(M-1):L; % discretize dielectric slab with M elements
phi = 0:pi/(2*M):pi/2;
%% Analytical solution
mu_r = ones(1,length(x));
mu_r(1:end-1) = 2-0.1j;
e_r = 4+(2-0.1j)*((1-y/L).^2);
e_r = padarray(e_r, [0 1],1,'post'); % pad array with '1' for free space
%plot(1:n,abs(e_r));
% compute kx
kx = k0*sqrt(repmat(mu_r.*e_r,length(x),1) - (ones(length(phi),1)*sin(phi).^2)');
%plot(x,abs(kx(end,:)),x,abs(kx(1,:)));
% compute eta
eta = (kx(:,[2 3:end]) - kx(:,[1 2:end-1]))./(kx(:,[2 3:end]) + kx(:,[1 2:end-1]));
eta(:,end) = (mu_r(end-1).*kx(:,end) - mu_r(end).*kx(:,end-1))./(mu_r(end-1).*kx(:,end) + mu_r(end).*kx(:,end-1));
R = zeros(length(phi),N);
% initialize R at the PEC
R(:,1) = -1;
% compute R
for i = 2:N
    R(:,i) = exp(2i*kx(:,i)*x(i)).*((eta(:,i-1) + R(:,i-1).*exp(-2i.*kx(:,i-1)*x(i))) ./ ...
             (1 + eta(:,i-1).*R(:,i-1).*exp(-2i*kx(:,i-1)*x(i))));
end
%% Numerical simulation
Rn = zeros(1,length(phi));
for k = 1:numel(phi)
    % initialize f
    f = zeros(N,1);
    % initialize the K matrix
    K = zeros(N,N);
    l = L/M;
    alpha = 1./mu_r;
    THETA = phi(k);
    beta = -(k0^2)*(e_r - alpha.*sin(THETA).^2);
    K(1,1) = alpha(1)/l + beta(1)*l/3;
    K(N,N) = alpha(M)/l + beta(M)*l/3;
    for i = 2:N-1
        K(i,i) = alpha(i-1)/l + beta(i-1)/3 + alpha(i)/l + beta(i)*l/3;
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
    reflection(k) = (Ez(end) - E_inc)/conj(E_inc);
end
%% Plot solution
subplot(1,2,2)
plot(phi*180/pi,abs(R(:,end)));
xlabel('\theta (degrees)');
ylabel('Reflection coefficient');
legend({'Analytical'});
subplot(1,2,1)
plot(phi*180/pi,abs(reflection))
xlabel('\theta (degrees)');
ylabel('Reflection coefficient');
legend({'Numerical'});
figure
plot(phi*180/pi,abs(R(:,end)),'-',phi*180/pi,abs(reflection),'--')
xlabel('\theta (degrees)');
ylabel('Reflection coefficient');
s = sprintf('FEM, %d cells',M);
legend({'Analytical',s});
toc