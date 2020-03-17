from casadi import *
from numpy import *
import pdb
import itertools
import numpy as np
from cvxpy import *

##### MY CODE ######
class nonlinearFTOCP(object):
	""" Finite Time Optimal Control Problem (FTOCP)
	Methods:
		- solve: solves the FTOCP given the initial condition x0 and terminal contraints
		- buildNonlinearProgram: builds the nonlinear program solved by the above solve methos
		- model: given x_t and u_t computes x_{t+1} = Ax_t + Bu_t

	"""

	def __init__(self, N, dt_mpc, bx, x_tightening, bu, u_tightening, Q, R, Qf, fixInput):
		# Define variables
		self.N  = N
		self.n  = 4
		self.d  = 1
		self.dt = dt_mpc
		self.bx = bx[0:self.n]
		self.bu = bu[0:self.d]
		self.x_tightening = x_tightening
		self.u_tightening = u_tightening
		self.Q = Q
		self.Qf = Qf
		self.R = R
		self.fixInput = fixInput

		self.buildNonlinearProgram()

	def solve(self, x0, verbose=False):
			# Set box constraints on states (here we constraint the last i steps of the horizon to be xf)

			if self.x_tightening.shape[1] == self.N:
				self.lbx = x0
				self.ubx = x0

				for i in range(0, self.N):
					self.lbx = self.lbx + (-self.bx-self.x_tightening[:,i]).tolist()
					self.ubx = self.ubx + ( self.bx+self.x_tightening[:,i]).tolist()

				for i in range(0, self.N):
					self.lbx = self.lbx + (-self.bu-self.u_tightening[i]).tolist()
					self.ubx = self.ubx + ( self.bu+self.u_tightening[i]).tolist()
			else:
				self.lbx = x0 + (-self.bx).tolist()*(self.N) + (-self.bu).tolist()*self.N
				self.ubx = x0 + ( self.bx).tolist()*(self.N) + ( self.bu).tolist()*self.N


			# Solve nonlinear programm
			sol = self.solver(lbx=self.lbx, ubx=self.ubx, lbg=self.lbg_dyanmics, ubg=self.ubg_dyanmics)

			# Check if the solution is feasible
			if (self.solver.stats()['success']) and (np.max(x0 - self.bx) < 0):
				self.feasible = 1
				x = sol["x"]
				self.xPred = np.array(x[0:(self.N+1)*self.n].reshape((self.n,self.N+1)))
				self.uPred = np.array(x[(self.N+1)*self.n:((self.N+1)*self.n + self.d*self.N)].reshape((self.d,self.N)))

				self.mpcInput = self.uPred[0][0]

			else:
				self.xPred = np.zeros((self.n, self.N+1))
				self.uPred = np.zeros((self.d, self.N))
				self.mpcInput = 0
				self.feasible = 0

				print("Unfeasible")
				# pdb.set_trace()

	def buildNonlinearProgram(self):
		# Define variables
		n  = self.n
		d  = self.d
		N  = self.N
		Q  = self.Q
		Qf = self.Qf
		R  = self.R
		fixInput = self.fixInput

		X      = SX.sym('X', n*(N+1));
		U      = SX.sym('X', d*N);

		# Define dynamic constraints
		constraint = []
		for i in range(0, N):
			constraint = vertcat(constraint, X[n*(i+1)+0] - (X[n*i+0] + self.dt*X[n*i+1]) ) 
			constraint = vertcat(constraint, X[n*(i+1)+1] - (X[n*i+1] + self.dt* ( (np.cos(X[n*i+2]) * ( -1.8*U[d*i] + 11.5*X[n*i+1] + 9.8*np.sin( X[n*i+2] ) ) - 10.9*U[d*i] + 68.4*X[n*i+1] - 1.2*X[n*i+3]**2*np.sin(X[n*i+2]) ) / ( np.cos(X[n*i+2]) - 24.7 ) ) )) 
			constraint = vertcat(constraint, X[n*(i+1)+2] - (X[n*i+2] + self.dt*X[n*i+3])) 
			constraint = vertcat(constraint, X[n*(i+1)+3] - (X[n*i+3] + self.dt* ( ( ( 9.3*U[d*i] - 58.8*X[n*i+1] )*np.cos(X[n*i+2]) + 38.6*U[d*i] - 234.5*X[n*i+1] - sin(X[n*i+2]) * (208.3 + X[n*i+3]**2*np.cos(X[n*i+2])) ) / (np.cos(X[n*i+2])**2 - 24.7) ) )) 

		# When the discretization of the nonlinear MPC is different from the dicretizaion of the input we need to constraint the intermedies inputs to be the same
		counter = 0
		totCounter = 0
		if fixInput > 0:
			for i in range(0,N-1):
				if counter < fixInput:
					constraint = vertcat(constraint, U[d*i] - U[d*i + 1] ) 
					counter += 1
					totCounter += 1
				else:
					counter = 0


			
		# Defining Cost (We will add stage cost later)
		cost = 0
		for i in range(0, N):
			cost += Q[0,0]*X[n*i+0]**2 + Q[1,1]*X[n*i+1]**2 + Q[2,2]*X[n*i+2]**2 + Q[3,3]*X[n*i+3]**2
			cost += R[0,0]*U[d*i]**2

		cost += Qf[0,0]*X[n*N+0]**2 + Qf[1,1]*X[n*N+1]**2 + Qf[2,2]*X[n*N+2]**2 + Qf[3,3]*X[n*N+3]**2
			


		# Set IPOPT options
		# opts = {"verbose":False,"ipopt.print_level":0,"print_time":0,"ipopt.mu_strategy":"adaptive","ipopt.mu_init":1e-5,"ipopt.mu_min":1e-15,"ipopt.barrier_tol_factor":1}#, "ipopt.acceptable_constr_viol_tol":0.001}#,"ipopt.acceptable_tol":1e-4}#, "expand":True}
		opts = {"verbose":False,"ipopt.print_level":0,"print_time":0}#, "ipopt.acceptable_constr_viol_tol":0.001}#,"ipopt.acceptable_tol":1e-4}#, "expand":True}
		nlp = {'x':vertcat(X,U), 'f':cost, 'g':constraint}
		self.solver = nlpsol('solver', 'ipopt', nlp, opts)

		# Set lower bound of inequality constraint to zero to force n*N state dynamics
		if fixInput > 0:
			self.lbg_dyanmics = [0]*(n*N) + [0]*totCounter
			self.ubg_dyanmics = [0]*(n*N) + [0]*totCounter
		else:
			self.lbg_dyanmics = [0]*(n*N)
			self.ubg_dyanmics = [0]*(n*N)

