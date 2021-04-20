function [ constrTightening ] = constrTighteningFunc( vector, Set )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
error = sdpvar(size(vector,2),1);

cost = -vector*error;
constraint = Set.A*error <= Set.b;

ops = sdpsettings('verbose',0,'solver','gurobi');
solution = solvesdp(constraint,cost,ops);

constrTightening = double(cost);
end

