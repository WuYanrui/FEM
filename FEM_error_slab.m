%% Project #2
% Authors: Blake Levy and Adedayo Lawal
clc;clear;
tic
%% set up constants from constants.m
constants
%% Set up geometry
m = 99; % Number of elements
n = m+1; % Number of nodes
L = 5*lamb0; % Length of slab is 5x free space wavelength
x_L = 0:L/m:L; % discretize dielectric slab with N nodes
y = 0:L/(m-1):L; % discretize dielectric slab with M elements
theta = 0:pi/(2*m):pi/2; % incident angle range from 0 to 90 degrees
degree = theta*180/pi;
er_slab = 4 + (2-1i*0.1)*((1-(y/L)).^2);
er_slab = padarray(er_slab, [0 1],1,'post'); % pad array with '1' for free space
e_slab = eps0*er_slab;
mur_slab = 2 - 1j*.1;
mu_slab = mu0*mur_slab;
%% create wavenumber for each slab AND free space
kx = zeros(length(theta),length(e_slab));
for i = 1:length(theta)
    for j = 1:length(e_slab)
        if j == length(e_slab)
            kx(i,j) = k0*sqrt(1 - sin(theta(i))^2);        
        else
            kx(i,j) = k0*sqrt(mur_slab*er_slab(j) - sin(theta(i))^2);
        end
    end
end
%% FEM
R = compute_analytic(m);
reflection = zeros(1,length(theta));
alpha_e = 1/mur_slab;
E0 = 1;
phi = zeros(length(theta),length(theta));
for T = 1:length(theta)
    %% Formation of the elemental K-matrix
    K_e = zeros(2,2);
    K = zeros(n,n); % This is a K-matrix for a SPECIFIED incidence angle    
    THETA = theta(T);    
    beta_e = -(k0^2)*(er_slab - alpha_e*sin(THETA)^2);    
    for e = 1:m
    l_e = x_L(e+1)-x_L(e);
        % Fill elemental K-matrix
        for i = 1:2
            for j = 1:2
                if i == j
                    K_e(i,j) = alpha_e/l_e + beta_e(e)*l_e/3;
                else
                    K_e(i,j) = -alpha_e/l_e + beta_e(e)*l_e/6;
                end
            end

        end

    % fill general K-matrix
    K(e,e) = K(e,e)+ K_e(1,1);
    K(e+1,e+1) = K(e+1,e+1)+ K_e(2,2);    
    K(e+1,e) = K(e+1,e)+ K_e(2,1);
    K(e,e+1) = K(e,e+1)+ K_e(1,2);

    
    
    end
    K(1,1) = 1;
    K(1,[2 end]) = 0;
    K([2 end], 1) = 0;
    %% Formation of the global b
    % form the global b vector
    q = 2*1j*k0*cos(THETA)*E0*exp(1j*k0*L*cos(THETA));
    gamma = 1j*k0*cos(THETA);
    b = zeros(n,1);
    % add the boundary conditions at L to the b vector and K matrix
    % respectively
    b(end) = q;
    K(end,end) = K(end,end)+ gamma;
    phi(T,:) = K\b;
    E_inc = E0*exp(1j*k0*L*cos(THETA));
    reflection(T) = (phi(T,end) - E_inc)/conj(E_inc);
end
% Plot Analytical solution
ks = k0*sqrt(4);
Z_in = j*sqrt(mu0/4/eps0)*tan(k0*sqrt(4)*x_L);
R_exact = (Z_in - sqrt(mu0/eps0))./(Z_in + sqrt(mu0/eps0));
Ez_exact = (1+R_exact(end))*exp(-1j*sqrt(4)*k0.*x_L);
error_L = abs(phi(1,:)-Ez_exact);
semilogy(x_L(2:end)/lamb0,error_L(2:end))
save('error.mat','error_L','x_L');
toc

