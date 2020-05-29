clc
clear all
close all

%% Initialize Parameters

T  = 10;         % Simulation time
dt = 0.001;      % Discretization time
dt_mpc = 0.001;    % MPC discretization time

useLowLevel  = 0; % 1 for low level active and 0 for low level non-active
constrActive = 1; % 1 for constraint on \theta and 0 for unconstrained case 

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
N_constraintTightening = 10;
N = 1000;

n = size(Q,2);
d = size(R,2);

%% System Matrices Converted in Discrete Time
 sys = ss(Alinear,Blinear,eye(n),0);
 sysd=c2d(sys,dt_mpc);
[Alin,Blin,Clin,Dlin]=ssdata(sysd);

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

for i = 1:N_constraintTightening
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
    for j = 1:(N/N_constraintTightening)
        x_tightening = [x_tightening, X_constrTightening{i}];
        u_tightening = [u_tightening, U_constrTightening{i}];    
    end
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

