function [ uMPC, xPred, uPred ] = FTOCP( x0, Alin, Blin, X, U, N, Q, R, Qf, x_goal, K, X_constrTightening, U_constrTightening )

% define variable
n = size(Alin, 2);
d = size(Blin, 2);
x = sdpvar(n, N+1);
u = sdpvar(d, N);

Fx = X.A; bx = X.b;
Fu = U.A; bu = U.b;

% constraints 
constraints = [ x(:,1) ==x0];

for i = 1:N
    X_constrTighteningVector = [X_constrTightening{i};X_constrTightening{i}];
    U_constrTighteningVector = [U_constrTightening{i};U_constrTightening{i}];
          
    constraints = [constraints;
                   x(:,i + 1) == (Alin-Blin*K)*x(:,i) + Blin*u(:,i);
                   Fx*x(:,i) <= bx + X_constrTighteningVector;
                   Fu*(u(:,i)- K*x(:,i)) <= bu + U_constrTighteningVector;];
end

% cost 
cost = 0;
for i = 1:N
    cost = cost + (x(:,i)-x_goal)'*Q*(x(:,i)-x_goal) + (u(:,i)-K*x(:,i))'*R*(u(:,i)-K*x(:,i));
end
cost = cost + (x(:,N+1)-x_goal)'*Qf*(x(:,N+1)-x_goal);

% solve QP        
ops = sdpsettings('verbose',1,'solver','gurobi');
solution = solvesdp(constraints,cost,ops);

if solution.problem ~= 0
    errorNotFeasible = 1
    xPred = zeros(n,N+1);
    uPred = zeros(d, N);
    uMPC = 0;
else
    xPred = double(x);
    uPred = double(u);
    uMPC = double(u(:,1)-K*x(:,1));
end
end

