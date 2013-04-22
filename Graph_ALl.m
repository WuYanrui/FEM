%% Plot all functions on one graph
clc;clear;
Rn_cube = load('cube.mat','Rn');
Rn_linear = load('linear.mat','reflection');
Rn_quad = load('quad.mat','Rn');
phi = load('quad.mat','phi');
R = load('quad.mat','R');
plot(phi.phi*180/pi,abs(R.R(:,end)),'k-',phi.phi*180/pi,abs(Rn_linear.reflection),'b:',...
    phi.phi*180/pi,abs(Rn_quad.Rn),'g-.',phi.phi*180/pi,abs(Rn_cube.Rn),'r--');
xlabel('Degrees (\theta)'); ylabel('Reflection Coefficient');
s = sprintf('Reflection Coefficient vs. Incident Angle for %d cells FEM ',length(phi.phi)-1);
title(s);
legend('Analytic','Linear','Quadratic','Cubic');