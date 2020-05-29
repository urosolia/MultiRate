function [ uBarrier,h, Vout ] = CBF( X, Xn, uMPC, currentTime, vmax, dotphimax, xmax, phimax)
%UNTITLED6 Summary of this CBF goes here
%   Detailed explanation goes here

% Exstrac States
x = X(1);      xn = Xn(1);
v = X(2);      vn = Xn(2);
phi = X(3);    phin = Xn(3);
dotphi = X(4); dotphin = Xn(4);


% Lyap function parameters
vcost = 1000.0;
dotphicost = 1000.0;
xcost = 10.0;
phicost = 10.0;



% Build QP
u = sdpvar(1);
gamma = sdpvar(1);

% Sys dynamcis 
dxdt = [v;
       (cos(phi) * ( -1.8*(u+uMPC) + 11.5*v + 9.8*sin( phi ) ) - 10.9*(u+uMPC) + 68.4*v - 1.2*dotphi^2*sin(phi) ) / ( cos(phi) - 24.7 );
       dotphi;
       ( ( 9.3*(u+uMPC) - 58.8*v )*cos(phi) + 38.6*(u+uMPC) - 234.5*v - sin(phi) * (208.3 + dotphi^2*cos(phi)) ) / (cos(phi)^2 - 24.7);
       vn;
       (11.5+68.4)/(-23.7)*vn + 9.8/(-23.7)*phin + 12.8/23.7*uMPC;
       dotphin;
       (-58.8-234.5)/(-23.7)*vn - 208.3/(-23.7)*phin - 47.9/23.7*uMPC];

% Barrier
h = 1 - (phi - phin)^2/phimax^2 - (v - vn)^2/vmax^2 - (x - xn)^2/xmax^2 - (dotphi - dotphin)^2/dotphimax^2;
 
 
dhdx =[         -(2*x - 2*xn)/xmax^2
                -(2*v - 2*vn)/vmax^2
          -(2*phi - 2*phin)/phimax^2
 -(2*dotphi - 2*dotphin)/dotphimax^2
                 (2*x - 2*xn)/xmax^2
                 (2*v - 2*vn)/vmax^2
           (2*phi - 2*phin)/phimax^2
  (2*dotphi - 2*dotphin)/dotphimax^2];

% Lyap
V = [dotphicost*(dotphi - dotphin)^2 + phicost*(phi - phin)^2 + vcost*(v - vn)^2 + xcost*(x - xn)^2];
 
 
dvdx =[          xcost*(2*x - 2*xn)
                 vcost*(2*v - 2*vn)
           phicost*(2*phi - 2*phin)
  dotphicost*(2*dotphi - 2*dotphin)
                -xcost*(2*x - 2*xn)
                -vcost*(2*v - 2*vn)
          -phicost*(2*phi - 2*phin)
 -dotphicost*(2*dotphi - 2*dotphin)];

% cost
cost = 100*u*u + 1*gamma^2;

% constraints
constraint = [dvdx'*dxdt - gamma <= -10.0*V;
              dhdx'*dxdt >= -10*h];

% solve QP
if V > 10^-4 || h< 1-10^-4
    ops = sdpsettings('verbose',0,'solver','gurobi');
    solution = solvesdp(constraint,cost,ops);

    if solution.problem ~= 0
        errorNotFeasible = 1
    end
    uBarrier = double(u);
else
    uBarrier = 0;
end
% output variables
Vout = V;

% plot to screen
[h, V, uBarrier, double(gamma), double(dvdx'*dxdt), currentTime]
here = 1;
end

