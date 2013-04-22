function [ K ] = compute_ke_cubic(alpha, beta, le)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    K = zeros(3,3);
    
    X = [-0.8611363116, -0.3399810436,0.3399810436,0.8611363116];
    W = [0.3478548451,0.6521451549,0.6521451549,0.3478548451]; 
    
    N1 = -(9/16)*(X+1/3).*(X-1/3).*(X-1);
    N2 = (27/16)*(X+1).*(X-1/3).*(X-1);
    N3 = -(27/16)*(X+1).*(X+1/3).*(X-1);
    N4 = (9/16)*(X+1).*(X+1/3).*(X-1/3);
    DN1 = 1/16 + (9/8)*X - (27*X.^2)/16;
    DN2 = -(27/16) - (9*X)/8 + (81*X.^2)/16;
    DN3 = 27/16 - (9*X)/8 - (81*X.^2)/16;
    DN4 = -(1/16) + (9*X)/8 + (27*X.^2)/16;
    
    K(1,1) = alpha*W*((2/le)*DN1.*DN1)'+ beta*W*((le/2)*N1.*N1)';
    K(1,2) = alpha*W*((2/le)*DN1.*DN2)'+beta*W*((le/2)*N1.*N2)';
    K(1,3) = alpha*W*((2/le)*DN1.*DN3)'+beta*W*((le/2)*N1.*N3)';
    K(1,4) = alpha*W*((2/le)*DN1.*DN4)'+beta*W*((le/2)*N1.*N4)';
    K(2,1) = K(1,2);
    K(2,2) = alpha*W*((2/le)*DN2.*DN2)'+beta*W*((le/2)*N2.*N2)';
    K(2,3) = alpha*W*((2/le)*DN2.*DN3)'+beta*W*((le/2)*N2.*N3)';
    K(2,4) = alpha*W*((2/le)*DN2.*DN4)'+beta*W*((le/2)*N2.*N4)';
    K(3,1) = K(1,3);
    K(3,2) = K(2,3);
    K(3,3) = alpha*W*((2/le)*DN3.*DN3)'+beta*W*((le/2)*N3.*N3)';
    K(3,4) = alpha*W*((2/le)*DN3.*DN4)'+beta*W*((le/2)*N3.*N4)';
    K(4,1) = K(1,4);
    K(4,2) = K(2,4);
    K(4,3) = K(3,4);
    K(4,4) = alpha*W*((2/le)*DN4.*DN4)'+beta*W*((le/2)*N4.*N4)';    
    
end

