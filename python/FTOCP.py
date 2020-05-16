import numpy as np
import pdb 
import scipy
from cvxpy import *
import time

class FTOCP(object):
	""" Finite Time Optimal Control Problem (FTOCP)
	Methods:
		- solve: solves the FTOCP given the initial condition x0, terminal contraints (optinal) and terminal cost (optional)
		- model: given x_t and u_t computes x_{t+1} = Ax_t + Bu_t

	"""
	def __init__(self, N, A, B, Q, Qf, R, K, Fx, bx, x_tightening, Fu, bu, u_tightening):
		# Define variables
		self.N = N # Horizon Length

		# System Dynamics (x_{k+1} = A x_k + Bu_k)
		self.A = A 
		self.B = B 
		self.n = A.shape[1]
		self.d = B.shape[1]

		# Cost (h(x,u) = x^TQx +u^TRu)
		self.Q = Q
		self.Qf = Qf
		self.R = R

		# Constraints
		self.Fx = Fx
		self.Fu = Fu
		self.bx = bx
		self.bu = bu
		self.x_tightening = x_tightening
		self.u_tightening = u_tightening
		self.K = K

		# Initialize Predicted Trajectory
		self.xPred = []
		self.uPred = []
		self.solverTime = []

	def solve(self, x0, verbose = False):
		"""This methos solve a FTOCP given:
			- x0: initial condition
			- SS: (optional) contains a set of state and the terminal constraint is ConvHull(SS)
			- Qfun: (optional) cost associtated with the state stored in SS. Terminal cost is BarycentrcInterpolation(SS, Qfun)
		""" 
		# Initialize Variables
		x = Variable((self.n, self.N+1))
		u = Variable((self.d, self.N))

		# State Constraints
		constr = [x[:,0] == x0[:]]
		for i in range(0, self.N):
			vec_X = np.squeeze(self.x_tightening[:,i])
			X_constrTighteningVector = np.concatenate((vec_X,vec_X));

			vec_U = np.squeeze(self.u_tightening[i])
			U_constrTighteningVector = np.array([vec_U, vec_U]);
			
			constr += [x[:,i+1] == (self.A - self.B*self.K)*x[:,i] + self.B*u[:,i],
						self.Fu * (u[:,i] - self.K*x[:,i]) <=  self.bu + U_constrTighteningVector,
						self.Fx * x[:,i] <=  self.bx + X_constrTighteningVector,]

		# Cost Function
		cost = 0
		for i in range(0, self.N):
			# Running cost h(x,u) = x^TQx + u^TRu
			cost += quad_form(x[:,i], self.Q) + norm(self.R**0.5*(u[:,i] - self.K*x[:,i]))**2
			# cost += norm(self.Q**0.5*x[:,i])**2 + norm(self.R**0.5*u[:,i])**2

		cost += quad_form(x[:,self.N], self.Qf) # If SS is not given terminal cost is quadratic


		# Solve the Finite Time Optimal Control Problem
		problem = Problem(Minimize(cost), constr)
		# if CVX == True:
		# 	problem.solve(verbose=verbose, solver=ECOS) # I find that ECOS is better please use it when solving QPs
		# else:
		# 	problem.solve(verbose=verbose)
		start = time.time()
		problem.solve(verbose=verbose, solver=ECOS) # I find that ECOS is better please use it when solving QPs
		end = time.time()
		self.solverTime = end - start 

		# Store the open-loop predicted trajectory


		if problem.status == 'optimal':
			self.feasible = 1
			self.xPred = x.value
			self.uPred = u.value	
		else:
			self.feasible = 0
			self.xPred = np.zeros((self.n, self.N+1))
			self.uPred = np.zeros((self.d, self.N+1))	
		
		self.mpcInput = self.uPred[:,0][0] - np.dot(self.K,self.xPred[:,0])

	def model(self, x, u):
		# Compute state evolution
		return (np.dot(self.A,x) + np.squeeze(np.dot(self.B,u))).tolist()





	

