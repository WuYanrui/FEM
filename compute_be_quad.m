function [ b ] = compute_be_quad(le,f)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    b = zeros(3,1);
    
    X = [-0.7745966692, 0.0000000000, 0.7745966692];
    W = [0.5555555556, 0.8888888889, 0.5555555556]; 
    
    N1 = X.*(X-1)/2;
    N2 = (1+X).*(1-X);
    N3 = X.*(X+1)/2;
    
    b(1) = (le/2)*f(1)*W*N1';
    b(2) = (le/2)*f(2)*W*N2';
    b(3) = (le/2)*f(3)*W*N3';
    
end

