%% Project #2
% Authors: Blake Levy and Adedayo Lawal
% (Quadratic Elements) 
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
phi = 0:pi/(2*M):pi/2;
%% Analytical solution
e_r = 4+(2-0.1j)*((1-y/L).^2);
e_r = padarray(e_r, [0 1],1,'post'); % pad array with '1' for free space
mu_r = ones(1,length(e_r));
mu_r(1:end-1) = 2-0.1j;
% compute kx
kx = k0*sqrt(repmat(mu_r.*e_r,length(e_r),1) - (ones(length(phi),1)*sin(phi).^2)');
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
Rn = zeros(1,length(phi));
for k = 1:numel(phi)
    % initialize f
    f = zeros(N,1);
    % initialize the K matrix
    K = zeros(N,N);
    alpha = 1./mu_r;
    THETA = phi(k);
    beta = -(k0^2)*(e_r - alpha.*sin(THETA).^2);
    le = x(3)-x(1);
    K(1,1) = alpha(1)*7/(3*le) + beta(1)*2*le/15;
    K(1,2) = -alpha(1)*8/(3*le) + beta(1)*le/15;    
    K(1,3) = alpha(1)/(3*le) - beta(1)*le/30;
    K(2,1) = K(1,2);    
    K(2,2) = alpha(1)*16/(3*le) + beta(1)*8*le/15;
    K(2,3) = -alpha(1)*8/(3*le) + beta(1)*le/15;    
    K(3,1) = K(1,3);
    K(3,2) = K(2,3);
    K(N,N) = alpha(M)*7/(3*le) + beta(M)*2*le/15;    
    for i = 3:2:N-2
        j = (i-1)/2+1;
        le = x(2*j+1)-x(2*j-1);
        K(i,i) = alpha(j)*7/(3*le) + beta(j)*2*le/15 + alpha(j-1)*7/(3*le) + beta(j-1)*2*le/15;
        K(i,i+1) =  -alpha(j)*8/(3*le) + beta(j)*le/15;
        K(i,i+2) = alpha(j)/(3*le) - beta(j)*le/30;
        K(i+1,i) =  K(i,i+1);
        K(i+1,i+1) = alpha(j)*16/(3*le) + beta(j)*8*le/15; 
        K(i+1,i+2) = -alpha(j)*8/(3*le) + beta(j)*le/15; 
        K(i+2,i) = K(i,i+2);
        K(i+2,i+1) = K(i+1,i+2);
    end
    % initialize the b vector
    b = zeros(N,1);
    b(1) = f(1)*le/6;
    b(2) = f(1)*2*le/3;
    b(N) = f(M)*le/6;
    for i =  3:2:N-2
        j = (i-1)/2+1;
        b(i) = f(j-1)*le/6 + f(j)*le/6;
        b(i+1) = f(j)*2*le/3;
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
    Rn(k) = (Ez(end) - E_inc)/conj(E_inc);
end
%% Plot solution
subplot(1,2,2)
plot(phi*180/pi,abs(R(:,end)));
xlabel('\theta (degrees)');
ylabel('Reflection coefficient');
legend({'Analytical'});
subplot(1,2,1)
plot(phi*180/pi,abs(Rn))
xlabel('\theta (degrees)');
ylabel('Reflection coefficient');
legend({'Numerical'});
figure
plot(phi*180/pi,abs(R(:,end)),'-',phi*180/pi,abs(Rn),'--')
xlabel('\theta (degrees)');
ylabel('Reflection coefficient');
s = sprintf('FEM (Quadratic Elements), %d cells',M);
legend({'Analytical',s});
toc