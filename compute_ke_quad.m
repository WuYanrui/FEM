function [ K ] = compute_ke_quad(alpha, beta, le)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    K = zeros(3,3);
    
    X = [-0.7745966692, 0.0000000000, 0.7745966692];
    W = [0.5555555556, 0.8888888889, 0.5555555556]; 
    
    N1 = X.*(X-1)/2;
    N2 = (1+X).*(1-X);
    N3 = X.*(X+1)/2;
    DN1 = X-1/2;
    DN2 = -2*X;
    DN3 = X+1/2; 
    
    K(1,1) = alpha*W*((2/le)*DN1.*DN1)'+ beta*W*((le/2)*N1.*N1)';
    K(1,2) = alpha*W*((2/le)*DN1.*DN2)'+beta*W*((le/2)*N1.*N2)';
    K(1,3) = alpha*W*((2/le)*DN1.*DN3)'+beta*W*((le/2)*N1.*N3)';
    K(2,1) = K(1,2);
    K(2,2) = alpha*W*((2/le)*DN2.*DN2)'+beta*W*((le/2)*N2.*N2)';
    K(2,3) = alpha*W*((2/le)*DN2.*DN3)'+beta*W*((le/2)*N2.*N3)';
    K(3,1) = K(1,3);
    K(3,2) = K(2,3);
    K(3,3) = alpha*W*((2/le)*DN3.*DN3)'+beta*W*((le/2)*N3.*N3)';
    
end

