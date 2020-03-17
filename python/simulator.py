import numpy as np
import pdb
import datetime

class Simulator():
    """ simulator"""
    def __init__(self, dt, steps):
        """Initialization
        dt:
        steps:
        """
        self.dt = dt
        self.steps = steps
    def sim(self, x, u):

        for i in range(0,self.steps):
            dxdt = [x[1],
                    (np.cos(x[2]) * ( -1.8*u + 11.5*x[1] + 9.8*np.sin( x[2] ) ) - 10.9*u + 68.4*x[1] - 1.2*x[3]**2*np.sin(x[2]) ) / ( np.cos(x[2]) - 24.7 ),
                    x[3],
                    ( ( 9.3*u - 58.8*x[1] )*np.cos(x[2]) + 38.6*u - 234.5*x[1] - np.sin(x[2]) * (208.3 + x[3]**2*np.cos(x[2])) ) / ((np.cos(x[2]))**2 - 24.7),];
            
            xNext = np.array(x) + self.dt/self.steps*np.array(dxdt)
            x = xNext

        return xNext.tolist()


