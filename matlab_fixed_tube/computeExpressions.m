%% This file is used to compute the dV/dt and dh/dt.

clc
clear all
close all
syms x xn xmax xcost
syms v vn vmax vcost
syms phi phin phimax phicost
syms dotphi dotphin dotphimax dotphicost
syms u uMPC
syms cross_h cross_V


dxdt = [v;
       (cos(phi) * ( -1.8*(u) + 11.5*v + 9.8*sin( phi ) ) - 10.9*(u) + 68.4*v - 1.2*dotphi^2*sin(phi) ) / ( cos(phi) - 24.7 );
       dotphi;
       ( ( 9.3*(u) - 58.8*v )*cos(phi) + 38.6*(u) - 234.5*v - sin(phi) * (208.3 + dotphi^2*cos(phi)) ) / (cos(phi)^2 - 24.7);
       vn;
       (11.5+68.4)/(-23.7)*v + 9.8/(-23.7)*dotphin + 12.7/23.7*uMPC;
       dotphin;
       (-58.8-234.5)/(-23.7)*v - 208.3/(-23.7)*dotphin - 47.9/23.7*uMPC];
   
A = [diff(dxdt,x), diff(dxdt,v), diff(dxdt,phi),diff(dxdt,dotphi)]

B = [diff(dxdt,u)]

h = 1 - (phi - phin)^2/phimax^2 - (x - xn)^2/xmax^2 - (v - vn)^2/vmax^2 - (dotphi - dotphin)^2/dotphimax^2

dhdx = [diff(h,x);
        diff(h,v);
        diff(h,phi);
        diff(h,dotphi);
        diff(h, xn);
        diff(h, vn);
        diff(h, phin);
        diff(h, dotphin)]
    
    

V = (v - vn)^2*vcost + (dotphi - dotphin)^2*dotphicost + (phi - phin)^2*phicost + (x - xn)^2*xcost


dvdx = [diff(V,x);
        diff(V,v);
        diff(V,phi);
        diff(V,dotphi);
        diff(V, xn);
        diff(V, vn);
        diff(V, phin);
        diff(V, dotphin)]