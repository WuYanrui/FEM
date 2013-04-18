%% Project #2
% Authors: Blake Levy and Adedayo Lawal
% Vectorized approach
clc;clear;
tic
%% set up constants from constants.m
constants

%% Set up geometry
M = 100; % Number of elements
n = M+1; % Number of nodes
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
R = zeros(length(phi),n);
% initialize R at the PEC
R(:,1) = -1;
% compute R
for i = 2:n
    R(:,i) = exp(2i*kx(:,i)*x(i)).*((eta(:,i-1) + R(:,i-1).*exp(-2i.*kx(:,i-1)*x(i))) ./ ...
             (1 + eta(:,i-1).*R(:,i-1).*exp(-2i*kx(:,i-1)*x(i))));
end
plot(phi*180/pi,abs(R(:,end)));
xlabel('\theta (degrees)');
ylabel('Reflection coefficient');
legend({'Analytical'});
toc