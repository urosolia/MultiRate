clc
clear all
close all

%% Initialize Parameters

T  = 5;         % Simulation time
dt = 0.001;      % Discretization time
dt_mpc = 0.5;    % MPC discretization time

useLowLevel  = 0; % 1 for low level active and 0 for low level non-active
constrActive = 0; % 1 for constraint on \theta and 0 for unconstrained case 

% Pick Initial Condition
x0    =  -1.0;                
v0    =   0.0;    
phi0  =   0.0;
dotphi0= -0.0;  
%% Initialize Linear Dyanmics (continous time from small angle approzimation)
Alinear =[0                     1              0  0;
          0   (11.5+68.4)/(-23.7)    9.8/(-23.7)  0;
          0                     0              0  1;
          0 (-58.8-234.5)/(-23.7)  -208.3/(-23.7) 0];

Blinear =[0;
          12.8/23.7;
          0;
         -47.9/23.7;];
     

 %% Initialize MPC and Barrrier parameters
Q  = diag([1,0.001,0.0001,0.0001]);
Qf = 10*diag([100,100,100,200]);
Qe = 100;
R = 0.001;
N = 10;

n = size(Q,2);
d = size(R,2);

%% System Matrices Converted in Discrete Time (Used in the mpc)
sys = ss(Alinear,Blinear,eye(n),0);
sysd=c2d(sys,dt_mpc); % Discretization dt_mpc
[Alin,Blin,Clin,Dlin]=ssdata(sysd);
              
xPred = zeros(size(Q,2),N+1);
%% Compute Constraint Tightening
% Q_e = diag(xmax, vmax, phimax, dotphimax)
xmax      = 0.2;
vmax      = 0.1;
phimax    = 0.05;
dotphimax = 0.1;
% Assuming the the error set is an ellispe with the above semi-axis
error_max = [xmax, vmax, phimax, dotphimax]'; 
% Input used at the low-level
input_low = [5];

% State error set (box over approximation which is all we need to computing the difference)
X_constrTightening = Polyhedron([eye(n);-eye(n)], [error_max;error_max]);
U_constrTightening = Polyhedron([eye(d);-eye(d)], [input_low;input_low]); % Input used by the low-level
%% Initialize Constraint Sets
if constrActive == 1
    state_cnstr = [10;
                   5;
                   0.78;
                   10*pi];
else
    state_cnstr = [10;
                   5;
                   2*pi;
                   10*pi];
end
           
input_cnstr = [25];

X = Polyhedron([eye(n); -eye(n)], [state_cnstr; state_cnstr]);
U = Polyhedron([eye(d); -eye(d)], [input_cnstr; input_cnstr]);


%% Matrices Used To Compute Evolution of the Nominal Model
% These matrices are used to compute the evolution of the nominal model
sys = ss(Alinear,Blinear,eye(n),0);
sysd=c2d(sys,dt);
[A_int,B_int,C_int,D_int]=ssdata(sysd);

%% Simulation Loop
time = 0; % simulation time
x_cl  = [x0; v0; phi0; dotphi0]; % set initial conditions

timeIndex = 1;
h(1) = 1; % Initialize h 
V(1) = 0; % Initialize V

mpcCounter = dt_mpc/dt;
while (time(end) <= T)
    
    if mpcCounter == dt_mpc/dt % every dt_mpc compute MPC solution
        mpcCounter = 0;
        [uMPC(timeIndex), xPred, uPred] = FTOCP( x_cl(:,timeIndex), Alin, Blin, X, U, N, Q, R, Qf, zeros(4,1), X_constrTightening, U_constrTightening, error_max);
        % Pick the nominal state based on MPC solution
        xn_cl(:,timeIndex) = xPred(:,1);
    else
        uMPC(timeIndex) = uMPC(timeIndex - 1);
    end
    
    
    % compute low-level control
    if useLowLevel == 1
        [uBarrier(timeIndex), h(timeIndex), V(timeIndex)] = CBF(x_cl(:,timeIndex), xn_cl(:,timeIndex), uMPC(timeIndex), time(timeIndex), vmax, dotphimax, xmax, phimax);
    else
        uBarrier(timeIndex) = 0;
    end
    
    % control action (summation low-level + mid-level, if low-level active)
    u(timeIndex) = useLowLevel * uBarrier(timeIndex) + uMPC(timeIndex);

    % interate non-linear system
    tspan=[0 dt];
    [~,NextStates]=ode45(@(t, x) sysDyn(t,x,u(timeIndex)),tspan, x_cl(:,timeIndex));
    x_cl(:,timeIndex + 1) = NextStates(end,:)';
    
    % integrate nominal system (Exact integration using matrices computed previously)
    xn_cl(:,timeIndex + 1) = A_int * xn_cl(:,timeIndex) + B_int * uMPC(timeIndex);

    % print to screen true, nominal and error states
    display(['Time: ', num2str(time(end)),'/', num2str(T)])
    display('True State | Nominal State | Difference ')
    display([x_cl(:,end), xn_cl(:,end), x_cl(:,end)-xn_cl(:,end)])
 
    % update time step
    time(timeIndex+1) = time(timeIndex) + dt;
    timeIndex = timeIndex + 1;
    mpcCounter = mpcCounter + 1;
    
    if x_cl(3,end)>1.25
        break
    end

end
display('Done with sim: plotting')
%% Plotting
figure 
title('states')
subplot(5,1,1)
plot(time, x_cl(1,:),'-o')
hold on
plot(time, xn_cl(1,:),'-rs')
xlim([0 T])
ylabel('x')
subplot(5,1,2)
plot(time, x_cl(2,:),'-o')
hold on
plot(time, xn_cl(2,:),'-rs')
xlim([0 T])
ylabel('v')
subplot(5,1,3)
plot(time, x_cl(3,:),'-o')
hold on
plot(time, xn_cl(3,:),'-rs')
xlim([0 T])
ylabel('phi')
subplot(5,1,4)
plot(time, x_cl(4,:),'-o')
hold on
plot(time, xn_cl(4,:),'-rs')
xlim([0 T])
ylabel('dotphi')
subplot(5,1,5)
hold on
plot(time(1:end-1), uBarrier,'-o')
plot(time(1:end-1), uMPC,'-rs')
plot(time(1:end-1), u,'-k*')
xlim([0 T ])
ylabel('u')
xlabel('time')

%%
figure 
subplot(2,1,1)
plot(time(1:end-1), h,'-o')
xlim([0 T])
ylabel('h')
subplot(2,1,2)
plot(time(1:end-1), V,'-o')
xlim([0 T])
ylabel('V')
xlabel('time')

%%
figure
hold on
pltUll  = plot(time(1:end-1), uBarrier,'-o');
pltUml  = plot(time(1:end-1), uMPC,'-rs');
pltUtot = plot(time(1:end-1), u,'-k*');
xlim([0 T ])
legend([pltUll, pltUml, pltUtot], {'u low-level', 'u mid-level', 'u total'})
ylabel('u')
xlabel('time')
