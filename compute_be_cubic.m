function [ b ] = compute_be_cubic(le,f)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    b = zeros(4,1);
    
    X = [-0.8611363116, -0.3399810436,0.3399810436,0.8611363116];
    W = [0.3478548451,0.6521451549,0.6521451549,0.3478548451]; 
    
    N1 = -(9/16)*(X+1/3).*(X-1/3).*(X-1);
    N2 = (27/16)*(X+1).*(X-1/3).*(X-1);
    N3 = -(27/16)*(X+1).*(X+1/3).*(X-1);
    N4 = (9/16)*(X+1).*(X+1/3).*(X-1/3);
    
    b(1) = (le/2)*f(1)*W*N1';
    b(2) = (le/2)*f(2)*W*N2';
    b(3) = (le/2)*f(3)*W*N3';
    b(4) = (le/2)*f(4)*W*N4';    
    
end

