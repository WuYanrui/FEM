%% Project #2
% Authors: Blake Levy and Adedayo Lawal
clc;clear;
tic
%% set up constants from constants.m
constants
trial = [50 100];
test = zeros(1,length(trial));
for i = 1: length(trial)
%% Set up geometry
m = trial(i); % Number of elements
n = m+1; % Number of nodes
    L = 5*lamb0; % Length of slab is 5x free space wavelength
    x = 0:L/m:L; % discretize dielectric slab with N nodes
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
        l_e = x(e+1)-x(e);
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

    % subplot(1,2,1)
    % plot(padarray(y,[0 1],L+L/m,'post')/L,abs(e_slab)); title('Permitivity profile in slab');
    % ylabel('permitivity');xlabel('distance from PEC (x/L)');
    % subplot(1,2,2)
    if m ==50
        reflection50 = reflection;
        phi50 = theta;
    end
        

end
    figure(1)
    hold on
    plot(degree,abs(R(:,end)),'k');
    plot(phi50*180/pi,abs(reflection50),'k-.');
    plot(degree,abs(reflection),'k--');


% figure(1)
subplot(1,2,1)
% Plot Analytic AND simulated results as a function of theta
plot(degree,abs(R(:,end)),'k',degree,abs(reflection),'k--')
legend('Analytical','simulated')
% Calculate Error (or difference)
diff = R(:,end) - transpose(reflection);
subplot(1,2,2)
plot(degree, 10*log10(abs(diff)));
title('Error (Difference) (dB) as function of Theta');
xlabel('Theta (\theta)');
ylabel('Error (dB)');
save('linear.mat','degree','R','reflection');
toc

