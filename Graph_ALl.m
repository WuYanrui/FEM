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
diff_L = abs(R.R(:,end)-transpose(Rn_linear.reflection));
diff_q = abs(R.R(:,end)-transpose(Rn_quad.Rn));
diff_c = abs(R.R(:,end)-transpose(Rn_cube.Rn));
figure(2)
hold on
plot(phi.phi*180/pi,log10(diff_L),'k');
plot(phi.phi*180/pi,log10(diff_q),'r');
plot(phi.phi*180/pi,log10(diff_c),'g');
s = sprintf('Error (dB) vs. Incident Angle for %d cells FEM ',length(phi.phi)-1);
title(s);
xlabel('Degrees (\theta)'); ylabel('Error');
legend('Linear','Quadratic','Cubic');