%% Project #2
% Authors: Blake Levy and Adedayo Lawal
clc;clear;
tic
%% set up constants from constants.m
constants
%% Set up geometry
m = 200; % Number of elements
n = m+1; % Number of nodes
Q = m+2;
L = 5*lamb0; % Length of slab is 5x free space wavelength
x = 0:L/(2*m):L; % discretize dielectric slab with (2N-1) nodes
y = 0:L/(m-1):L; % discretize dielectric slab with M elements
theta = 0:pi/(2*m):pi/2; % incident angle range from 0 to 90 degrees
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
R = zeros(length(theta),n);
eta = zeros(length(theta),m);
R(:,1) = -1;
% Compute ETA_m 
for i = 1:length(theta)
    for j = 1:size(eta,2)
        if j == size(eta,2)
            eta(i,j) = (mur_slab*kx(i,j+1) - kx(i,j))...
                /(mur_slab*kx(i,j+1) + kx(i,j));         
        else
            eta(i,j) = (kx(i,j+1) - kx(i,j))...
                /(kx(i,j+1) + kx(i,j));   
        end
%         eta(i,j) = (mu_slab*kx(i,j) - mu_slab*kx(i,j-1))...
%             /(mu_slab*kx(i,j) + mu_slab*kx(i,j-1));       
    end
end
% Compute R_m
for i = 1:length(theta)
    for j = 1:size(eta,2)
        
           R(i,j+1) = exp(2*1j*kx(i,j+1)*x(2*j+1))*...
            ((eta(i,j) + R(i,j)*exp(-2*1j*kx(i,j)*x(2*j+1)))...
            /(1+(eta(i,j)*R(i,j)*exp(-2*1j*kx(i,j)*x(2*j+1)))));     
        
    end
end


reflection = zeros(1,length(theta));
alpha_e = 1/mur_slab;

E0 = 1;
phi = zeros(length(theta),Q);
for T = 1:length(theta)
    %% Formation of the elemental K-matrix
    K_e = zeros(3,3);
    K = zeros(Q,Q); % This is a K-matrix for a SPECIFIED incidence angle    
    THETA = theta(T);    
    beta_e = -(k0^2)*(er_slab - alpha_e*sin(THETA)^2);    
    for e = 1:m
        l_e = x(2*e+1)-x(2*e-1);
            % Fill elemental K-matrix
        K_e(1,1) = alpha_e*7/(3*l_e) + beta_e(e)*2*l_e/15;
        K_e(1,2) = -alpha_e*8/(3*l_e) + beta_e(e)*l_e/15;
        K_e(1,3) = alpha_e/(3*l_e) - beta_e(e)*l_e/30;    
        K_e(2,2) = alpha_e*16/(3*l_e) + beta_e(e)*8*l_e/15;
        K_e(2,1) = K_e(1,2);   
        K_e(2,3) = -alpha_e*8/(3*l_e) + beta_e(e)*l_e/15;
        K_e(3,2) = K_e(2,3);
        K_e(3,1) = K_e(1,3);
        K_e(3,3) = K_e(1,1);

        % fill general K-matrix
        K(e,e) = K(e,e)+ K_e(1,1);
        K(e,e+1) = K(e,e+1)+ K_e(1,2);
        K(e,e+2) = K(e,e+2)+ K_e(1,3);   
        K(e+1,e) = K(e+1,e)+ K_e(2,1);
        K(e+1,e+1) = K(e+1,e+1)+ K_e(2,2);      
        K(e+1,e+2) = K(e+1,e+2)+ K_e(2,3);  
        K(e+2,e) = K(e+2,e)+ K_e(3,1);  
        K(e+2,e+1) = K(e+2,e+1)+ K_e(3,2);      
        K(e+2,e+2) = K(e+2,e+2)+ K_e(3,3);    
    end
    K(1,1) = 1;
    K(1,[2 end]) = 0;
    K([2 end], 1) = 0;
    %% Formation of the global b
    % form the global b vector
    q = 2*1j*k0*cos(THETA)*E0*exp(1j*k0*L*cos(THETA));
    gamma = 1j*k0*cos(THETA);
    b = zeros(Q,1);
    % add the boundary conditions at L to the b vector and K matrix
    % respectively
    b(end) = q;
    K(end,end) = K(end,end)+ gamma;
    phi(T,:) = K\b;
    E_inc = E0*exp(1j*k0*L*cos(THETA));
    reflection(T) = (phi(T,end) - E_inc)/conj(E_inc);
end
% Plot Analytical solution

% subplot(1,2,1)
% plot(padarray(y,[0 1],L+L/m,'post')/L,abs(e_slab)); title('Permitivity profile in slab');
% ylabel('permitivity');xlabel('distance from PEC (x/L)');
% subplot(1,2,2)
figure(2)
plot(theta*180/pi,abs(R(:,end)),'k',theta*180/pi,abs(reflection),'k--')
% subplot(1,2,1)
% plot(theta*180/pi,abs(reflection))
legend('Analytical','simulated')
toc

