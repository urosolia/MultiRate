clc
clear all
close all

%% Initialize Parameters

T  = 10;         % Simulation time
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

%% System Matrices Converted in Discrete Time
 sys = ss(Alinear,Blinear,eye(n),0);
 sysd=c2d(sys,dt_mpc);
[Alin,Blin,Clin,Dlin]=ssdata(sysd);
              
xPred = zeros(size(Q,2),N+1);


%% Compute Constraint Tightening
% Q_e = diag(xmax, vmax, phimax, dotphimax)
xmax      = 0.2;
vmax      = 0.1;
phimax    = 0.01;
dotphimax = 0.1;

% Compute the propagation for the vx, phi, dotphi (the position is
% independent!)
[K_small,~,~] = dlqr(Alin(2:end,2:end),Blin(2:end),eye(3),1);
K = [0, K_small];
Astable = Alin(2:end,2:end) - Blin(2:end)*K_small;
Alin_small = Astable; % pick 3X3 matrix

bVector = [vmax; phimax; dotphimax];
Se = Polyhedron([eye(3);-eye(3)], [bVector;bVector]);
Se.computeVRep()
Se.computeHRep()
Se.computeVRep()

% Initialize propagation for the position (Independent!)
Sex = Polyhedron([eye(1);-eye(1)], [xmax;xmax]);

E{1} = 0*Se;
E{1}.computeVRep()
Ex{1} = 0*Sex;
Ex{1}.computeVRep()

X_constrTightening{1} = zeros(4,1);
U_constrTightening{1} = 0;

for i = 1:N
    E{i+1}  = Alin_small*E{i} + Se;
    E{i+1}.computeVRep()
    E{i+1}.computeHRep()
    E{i+1}.computeVRep()
    
    Ex{i+1} = Ex{i} + Alin(1,2:end)*E{i} + Sex
    Ex{i+1}.computeVRep()
    Ex{i+1}.computeHRep()
    Ex{i+1}.computeVRep()
    
    U_constrTightening{i+1} = constrTighteningFunc( K_small, E{i+1} );
    X_constrTightening{i+1} = zeros(4,1);
    X_constrTightening{i+1}(1,1) = max(Ex{i+1}.V);
    X_constrTightening{i+1}(2,1) = constrTighteningFunc( [1,0,0], E{i+1} );
    X_constrTightening{i+1}(3,1) = constrTighteningFunc( [0,1,0], E{i+1} );
    X_constrTightening{i+1}(4,1) = constrTighteningFunc( [0,0,1], E{i+1} );
end

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
 sys = ss(Alinear,Blinear,eye(n),0);
 sysd=c2d(sys,dt);
[A_int,B_int,C_int,D_int]=ssdata(sysd);

%% Save quantities to txt file for python implementation
% state dim
stateInputDim =[n, d, N];
dlmwrite(strcat('../python/parameter_dt_mpc_',num2str(dt_mpc),'/stateInputDim.txt'),stateInputDim,'delimiter','\t')

% mpc matrices
MPC_A_B = [Alin, Blin];
dlmwrite(strcat('../python/parameter_dt_mpc_',num2str(dt_mpc),'/MPC_A_B.txt'),MPC_A_B,'delimiter','\t')

% int matrices
int_A_B = [A_int, B_int];
dlmwrite(strcat('../python/parameter_dt_mpc_',num2str(dt_mpc),'/int_A_B.txt'),int_A_B,'delimiter','\t')

MPC_X_cntr = [X.A, X.b]; 
MPC_U_cntr = [U.A, U.b];
dlmwrite(strcat('../python/parameter_dt_mpc_',num2str(dt_mpc),'/MPC_X_cntr.txt'),MPC_X_cntr,'delimiter','\t')
dlmwrite(strcat('../python/parameter_dt_mpc_',num2str(dt_mpc),'/MPC_U_cntr.txt'),MPC_U_cntr,'delimiter','\t')

% cost matrices
Qdiag = diag(Q);
Qfdiag = diag(Qf);
Rdiag = diag(R);
costMatrices=[Qdiag', Qfdiag', Rdiag'];
dlmwrite(strcat('../python/parameter_dt_mpc_',num2str(dt_mpc),'/costMatrices.txt'),costMatrices,'delimiter','\t')

% constraint tightening and K
x_tightening =[];
u_tightening =[];
for i=1:size(X_constrTightening,2)
    x_tightening = [x_tightening, X_constrTightening{i}];
    u_tightening = [u_tightening, U_constrTightening{i}];    
end
dlmwrite(strcat('../python/parameter_dt_mpc_',num2str(dt_mpc),'/x_tightening.txt'),x_tightening,'delimiter','\t')
dlmwrite(strcat('../python/parameter_dt_mpc_',num2str(dt_mpc),'/u_tightening.txt'),u_tightening,'delimiter','\t')
dlmwrite(strcat('../python/parameter_dt_mpc_',num2str(dt_mpc),'/K.txt'),K,'delimiter','\t')

% save CLF_CBF parameters
CLF_CBF_parameters = [xmax, vmax, phimax, dotphimax];
dlmwrite(strcat('../python/parameter_dt_mpc_',num2str(dt_mpc),'/CLF_CBF_parameters.txt'),CLF_CBF_parameters,'delimiter','\t')


% save parameters
parameters = [T, dt, dt_mpc];
dlmwrite(strcat('../python/parameter_dt_mpc_',num2str(dt_mpc),'/parameters.txt'),parameters,'delimiter','\t')

%% Simulation Loop
time = 0; % simulation time
x_cl  = [x0; v0; phi0; dotphi0]; % set initial conditions
% xn_cl = [x0; v0; phi0; dotphi0]; % set initial conditions for nominal model
timeIndex = 1;
h(1) = 1; % Initialize h 
V(1) = 0; % Initialize V

mpcCounter = dt_mpc/dt;
while (time(end) <= T)
    
    if mpcCounter == dt_mpc/dt; % every dt_mpc compute MPC solution
        mpcCounter = 0;
        [uMPC(timeIndex), xPred, uPred] = FTOCP( x_cl(:,timeIndex), Alin, Blin, X, U, N, Q, R,Qf, zeros(4,1), K, X_constrTightening, U_constrTightening);
        
        xn_cl(:,timeIndex) = xPred(:,1);
    else
        uMPC(timeIndex) = uMPC(timeIndex - 1);
    end
    
    
    % compute barrier
    [uBarrier(timeIndex), h(timeIndex), V(timeIndex)] = CBF(x_cl(:,timeIndex), xn_cl(:,timeIndex), uMPC(timeIndex), time(timeIndex), vmax, dotphimax, xmax, phimax);
    
    % control action
    u(timeIndex) = useLowLevel * uBarrier(timeIndex) + uMPC(timeIndex);

    % interate system
    tspan=[0 dt];
    [~,NextStates]=ode45(@(t, x) sysDyn(t,x,u(timeIndex)),tspan, x_cl(:,timeIndex));
    x_cl(:,timeIndex + 1) = NextStates(end,:)';
    
    % integrate nominal system
    xn_cl(:,timeIndex + 1) = A_int * xn_cl(:,timeIndex) + B_int * uMPC(timeIndex);

    % print to screen true, nominal and error states
    [x_cl(:,end), xn_cl(:,end), x_cl(:,end)-xn_cl(:,end)]
 
    % update time
    time(timeIndex+1) = time(timeIndex) + dt;
    timeIndex = timeIndex + 1;
    mpcCounter = mpcCounter + 1;
    
    if x_cl(3,end)>1.25
        break
    end

end

%% Plotting
figure 
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
%%
figure
hold on
plot(time(1:end-1), uBarrier,'-o')
plot(time(1:end-1), uMPC,'-rs')
plot(time(1:end-1), u,'-k*')
xlim([0 T ])
ylabel('u')

%%
if constrActive == 1;
    cnstrString = '_withConstr.mat';
else
    cnstrString = '.mat';
end;

if useLowLevel == 1
    save(strcat('dt_mpc',num2str(dt_mpc*10),'_withBarrier',cnstrString))
else
    save(strcat('dt_mpc',num2str(dt_mpc*10),'_withoutBarrier',cnstrString))
end