function [ dxdt ] = sysDyn( t, x, u )
% System dynamics: x(1) = x, x(2) = v, x(3) = phi, x(4) = \dot phi              
dxdt = [x(2); 
       (cos(x(3)) * ( -1.8*u + 11.5*x(2) + 9.8*sin( x(3) ) ) - 10.9*u + 68.4*x(2) - 1.2*x(4)^2*sin(x(3)) ) / ( cos(x(3)) - 24.7 );
       x(4);
       ( ( 9.3*u - 58.8*x(2) )*cos(x(3)) + 38.6*u - 234.5*x(2) - sin(x(3)) * (208.3 + x(4)^2*cos(x(3))) ) / (cos(x(3))^2 - 24.7);];

end

