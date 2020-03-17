import numpy as np
import pdb 
import scipy
from cvxpy import *

class CBF_CLF(object):
	""" CBF_CLF
	"""
	def __init__(self, vmax, dotphimax, xmax, phimax):
		self.vmax      = vmax
		self.dotphimax = dotphimax
		self.xmax      = xmax
		self.phimax    = phimax

	def solve(self, X, Xn, uMPC, verbose=False):
		vmax      = self.vmax
		dotphimax = self.dotphimax
		xmax      = self.xmax
		phimax    = self.phimax

		# Exstrac States
		x = X[0]      
		v = X[1]      
		phi = X[2]    
		dotphi = X[3] 

		xn = Xn[0]
		vn = Xn[1]
		phin = Xn[2]
		dotphin = Xn[3]


		# Lyap function parameters
		vcost = 1000.0;
		dotphicost = 1000.0;
		xcost = 10.0;
		phicost = 10.0;

		# Build QP
		uCBF = Variable(1)
		gamma = Variable(1)

		# Sys dynamcis 
		dxdt = [v,
				(np.cos(phi) * ( -1.8*(uCBF+uMPC) + 11.5*v + 9.8*np.sin( phi ) ) - 10.9*(uCBF+uMPC) + 68.4*v - 1.2*dotphi**2*np.sin(phi) ) / ( np.cos(phi) - 24.7 ),
				dotphi,
				( ( 9.3*(uCBF+uMPC) - 58.8*v )*np.cos(phi) + 38.6*(uCBF+uMPC) - 234.5*v - np.sin(phi) * (208.3 + dotphi**2*np.cos(phi)) ) / (np.cos(phi)**2 - 24.7),
				vn,
				(11.5+68.4)/(-23.7)*v + 9.8/(-23.7)*phin + 12.8/23.7*uMPC,
				dotphin,
				(-58.8-234.5)/(-23.7)*v - 208.3/(-23.7)*phin - 47.9/23.7*uMPC]

		# Barrier
		h = 1 - (phi - phin)**2/phimax**2 - (v - vn)**2/vmax**2 - (x - xn)**2/xmax**2 - (dotphi - dotphin)**2/dotphimax**2
 
 
		dhdx =[	-(2*x - 2*xn)/xmax**2,
				-(2*v - 2*vn)/vmax**2,
				-(2*phi - 2*phin)/phimax**2,
				-(2*dotphi - 2*dotphin)/dotphimax**2,
				(2*x - 2*xn)/xmax**2,
				(2*v - 2*vn)/vmax**2,
				(2*phi - 2*phin)/phimax**2,
				(2*dotphi - 2*dotphin)/dotphimax**2]

		# Lyap
		V = dotphicost*(dotphi - dotphin)**2 + phicost*(phi - phin)**2 + vcost*(v - vn)**2 + xcost*(x - xn)**2
 
 
		dvdx =[xcost*(2*x - 2*xn),
				vcost*(2*v - 2*vn),
				phicost*(2*phi - 2*phin),
				dotphicost*(2*dotphi - 2*dotphin),
				-xcost*(2*x - 2*xn),
				-vcost*(2*v - 2*vn),
				-phicost*(2*phi - 2*phin),
				-dotphicost*(2*dotphi - 2*dotphin)]

		# cost
		cost = 100*uCBF**2 + gamma**2;

		# constraints
		# constraint = [np.dot(dvdx,dxdt) - gamma <= -10.0*V,
		# 				np.dot(dhdx,dxdt) >= -10*h]
		constraint = [np.dot(dvdx,dxdt) - gamma <= -10.0*V,
					  np.dot(dhdx,dxdt) >= -10*h]

		if (V > 1e-4) and (h<1-1e-4):
			problem = Problem(Minimize(cost), constraint)
			problem.solve(verbose=verbose) # I find that ECOS is better please use it when solving QPs

			self.uCBF = uCBF.value[0]	
		else:
			self.uCBF = 0.0

	def model(self, x, u):
		# Compute state evolution
		return (np.dot(self.A,x) + np.squeeze(np.dot(self.B,u))).tolist()





	

