%% convergence
x = [50 100 200 500];
linear = [-1.163 -1.717 -2.33  -2.88];
quad = [-1.983 -3.274 -4.481 -6.071];
cubic = [-2.825 -3.08 -3.379 -3.776];
plot(x,linear,'k',x,quad,'r',x,cubic,'g');
xlabel('Number of Nodes');
ylabel('Log(Error)');
title('Log(Error) as a function of nodes');
legend('Linear','Quadratic','Cubic');