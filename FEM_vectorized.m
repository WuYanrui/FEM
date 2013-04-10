%% Project #2
% Authors: Blake Levy and Adedayo Lawal
clc;clear;
%% set up constants from constants.m
constants
%% Set up geometry
m = 50; % Number of elements
n = m+1; % Number of nodes
L = 5*lamb0; % Length of slab is 5x free space wavelength
x = 0:L/m:L; % discretize dielectric slab with M elements
theta = 0:pi/2/m:pi/2; % incident angle range from 0 to 90 degrees
er_slab = 4 + (2-1j*.1)*((1-x/L).^2);
e_slab = eps0*er_slab;
mur_slab = 2 - 1j*.1;
mu_slab = mu0*mur_slab;
kx = zeros(length(theta),length(e_slab));
for i = 1:length(theta)
    for j = 1:length(e_slab)
        kx(i,j) = k0*sqrt(mu_slab*e_slab(j) - sin(theta(i))^2);
    end
end
% kx = k0*sqrt(mu_slab*e_slab - sin(theta).^2);
%% FEM
eta = zeros(1,m+1);
R = zeros(length(theta),length(e_slab));
eta = zeros(length(theta),length(e_slab)+1);
eta(:,end) = sqrt(mu0/eps0);
R(:,1) = -1;





eta = (mu_slab*kx(:,[2:end]) - mu_slab*kx(:,[1:end-1]))...
    ./(mu_slab*kx(:,[2:end]) + mu_slab*kx(:,[1:end-1]));
for i = 2:length(e_slab)+1
   R(:,i) = exp(2*1j*kx(:,i)*x(i)).*(eta(:,i) + R(:,i-1).*exp(-2*1j*kx(:,i-1).*x(i)))...
       ./(1 + eta(:,i).*R(:,i-1).*exp(-2*1j*kx(:,i-1).*x(i)));
end










% for i = 2:m+1
%     if i == m+1
%         eta(i) = sqrt(mu0/eps0);
%     else
%         eta(i) = (mu_slab*kx(i) - mu_slab*kx(i-1))/(mu_slab*kx(i)...
%             + mu_slab*kx(i-1));
%     end
%         R(i) = exp(2*1j*kx(i)*x(i))*((eta(i) + R(i-1)*exp(-2*1j*kx(i-1)*x(i)))...
%             /(1+(eta(i)*R(i-1)*exp(-2*1j*kx(i-1)*x(i)))));
% 
% end

% eta = mu_slab*kx - mu_slab*kx(2:end-1);

subplot(1,2,1)
plot(x/L,abs(e_slab)); title('Permitivity profile in slab');
ylabel('permitivity');xlabel('distance from PEC (x/L)');
subplot(1,2,2)
plot(theta*180/pi,abs(R))




