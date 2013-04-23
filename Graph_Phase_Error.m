%% Project #2
% Authors: Blake Levy and Adedayo Lawal
% Compute the Phase Error
% Generate the lamda/h log space vector
%x = logspace(1,2,1000);
x = linspace(10,100,100);
%% Linear Elements
pe1 = (x.*acos((6-2*((2*pi)./x).^2)./(6+((2*pi)./x).^2))/(2*pi))-1;
loglog(x,abs(pe1)*360)
xlabel('\lambda/h');
ylabel('Phase error (degrees)')
%% Quadratic Elements
hold on
pe2 = ((x./2).*acos((15-26*((2*pi)./x).^2+3*((2*pi)./x).^4)./(15+4*((2*pi)./x).^2+((2*pi)./x).^4))/(2*pi))-1;
% plot the phase error
loglog(x,abs(pe2)*360)
%% Cubic Elements
pe3 = ((x./3).*acos((2800-11520*((2*pi)./x).^2+4860*((2*pi)./x).^4-324*((2*pi)./x).^6)./(2800+1080*((2*pi)./x).^2+270*((2*pi)./x).^4+81*((2*pi)./x).^6))/(2*pi))-1;
% plot the phase error
loglog(x,abs(pe3)*360)
hold off