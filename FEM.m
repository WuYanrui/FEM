%% Project #2
% Authors: Blake Levy and Adedayo Lawal
clc;clear;
%% set up constants from constants.m
constants
%% Set up geometry
m = 50; % Number of elements
n = m+1; % Number of nodes
L = 5*lamb0; % Length of slab is 5x free space wavelength
x = 0:L/m:L; % discretize dielectric slab with N nodes
y = 0:L/(m-1):L; % discretize dielectric slab with M elements
theta = 0:pi/2/m:pi/2; % incident angle range from 0 to 90 degrees
er_slab = 4 + (2-1i*0.1)*((1-y/L).^2);
er_slab = padarray(er_slab, [0 1],1,'post');
e_slab = eps0*er_slab;
mur_slab = 2 - 1j*.1;
mu_slab = mu0*mur_slab;
%% create wavenumber for each slab AND free space
kx = zeros(length(theta),length(e_slab));
for i = 1:length(theta)
    for j = 1:length(e_slab)
        if j == length(e_slab)
            kx(i,j) = k0*sqrt(mu0*e_slab(j) - sin(theta(i))^2);        
        else
            kx(i,j) = k0*sqrt(mu_slab*e_slab(j) - sin(theta(i))^2);
        end
    end
end
%% FEM
R = zeros(length(theta),n);
eta = zeros(length(theta),n);
R(:,1) = -1;
% Compute ETA_m 
for i = 1:length(theta)
    for j = 2:n
        if j == length(e_slab)
            eta(i,j) = (mu_slab*kx(i,j) - mu0*kx(i,j-1))...
                /(mu_slab*kx(i,j) + mu0*kx(i,j-1));         
        else
            eta(i,j) = (kx(i,j) - kx(i,j-1))...
                /(kx(i,j) + kx(i,j-1));   
        end
%         eta(i,j) = (mu_slab*kx(i,j) - mu_slab*kx(i,j-1))...
%             /(mu_slab*kx(i,j) + mu_slab*kx(i,j-1));       
    end
end
% Compute R_m
for i = 1:length(theta)
    for j = 2:length(e_slab)
        
           R(i,j) = exp(2*1j*kx(i,j)*x(j))*...
            ((eta(i,j) + R(i,j-1)*exp(-2*1j*kx(i,j-1)*x(j)))...
            /(1+(eta(i,j)*R(i,j-1)*exp(-2*1j*kx(i,j-1)*x(j)))));     
        
    end
end

subplot(1,2,1)
plot(padarray(y,[0 1],L+L/m,'post')/L,abs(e_slab)); title('Permitivity profile in slab');
ylabel('permitivity');xlabel('distance from PEC (x/L)');
subplot(1,2,2)
plot(theta*180/pi,abs(R(:,end)))




