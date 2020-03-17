import pdb
import numpy as np
from simulator import Simulator
from FTOCP import FTOCP
from nonlinearFTOCP import nonlinearFTOCP
from CBF_CLF import CBF_CLF
import pdb
import dill
import matplotlib.pyplot as plt
from tempfile import TemporaryFile
import copy
import datetime
import os

# Check if a storedData folder exist.	
if not os.path.exists('storedData'):
	os.makedirs('storedData')

# Initialization Parameters
Example     = 1 # Pich between example 1 and 2 in the paper
useLowLevel = 0 # 1 low level active and 0 low level not active
linearMPC   = 0 # 1 linear highlevel MPC and 0 nonlinear high level MPC

# selet folder for reading data and set nonlinear MPC parameters
if Example == 1:
	folderToLoad = 'parameter_dt_mpc_0.5'
	
	# dt_nonlinear = 0.01 # Discretization time NMPC model ---> 1/0.01 model update rate
	# N_nonlinear = 500   # Horizon length 
	# fixInput = 49       # Number of consecutive inputs kept constant --> 1/0.5 input update rate 

	dt_nonlinear = 0.05 # Discretization time NMPC model ---> 1/0.05 model update rate
	N_nonlinear = 100   # Horizon length 
	fixInput = 9        # Number of consecutive inputs kept constant --> 1/0.5 input update rate 

else:
	folderToLoad = 'parameter_dt_mpc_0.1'
	dt_nonlinear = 0.1
	N_nonlinear = None
	fixInput = 0 # No need to fix inputs as the discretization time of the nonlinear mpc is the same as the linear one


if linearMPC == 0: # If using nonlinear model deactivate low level controller
	useLowLevel = 0
############################# Load Matlab Parameters #####################################

# Load Dim and MPC horizon
lines = np.loadtxt(folderToLoad+"/stateInputDim.txt", unpack=False)
n = int(lines[0])
d = int(lines[1])
if (linearMPC == 1) or (N_nonlinear == None):
	N = int(lines[2])
else:
	N = N_nonlinear

# Load MPC matrices
lines = np.loadtxt(folderToLoad+"/MPC_A_B.txt", unpack=False)
Alin  = lines[:,0:n]
Blin  = lines[:,n:(n+d)]

# Load int matrices
lines = np.loadtxt(folderToLoad+"/int_A_B.txt", unpack=False)
A_int  = lines[:,0:n]
B_int  = lines[:,n:(n+d)]

# Load Constraint set, constraint tightening and gain K
lines = np.loadtxt(folderToLoad+"/MPC_X_cntr.txt", unpack=False)
Fx    = lines[:,0:n]
bx    = np.squeeze(lines[:,n:(n+1)])
lines = np.loadtxt(folderToLoad+"/MPC_U_cntr.txt", unpack=False)
Fu    = lines[:,0:d]
bu    = np.squeeze(lines[:,(d):(d+1)])
lines = np.loadtxt(folderToLoad+"/u_tightening.txt", unpack=False)
u_tightening    = lines
lines = np.loadtxt(folderToLoad+"/x_tightening.txt", unpack=False)
x_tightening    = lines
lines = np.loadtxt(folderToLoad+"/K.txt", unpack=False)
K     = lines

# Load Cost Matrices
lines = np.loadtxt(folderToLoad+"/costMatrices.txt", unpack=False)
Q     = np.diag(lines[0:n])
Qf    = np.diag(lines[n:(2*n)])
R     = np.diag(lines[(2*n):(2*n+d)])

# Load CLF CBF parameters
lines = np.loadtxt(folderToLoad+"/CLF_CBF_parameters.txt", unpack=False)
xmax      = lines[0]
vmax      = lines[1]
phimax    = lines[2]
dotphimax = lines[3]

# Load Parameters
lines  = np.loadtxt(folderToLoad+"/parameters.txt", unpack=False)
T      = lines[0]
dt     = lines[1]
# if linearMPC == 1:
# 	dt_mpc = lines[2]
# else:
# 	dt_mpc = dt_nonlinear
dt_mpc = lines[2]

############################# Initialize Objects #####################################
steps = 10
simulator = Simulator(dt, steps)
xcl     = [[-1,0,0,0]]
xncl    = [[-1,0,0,0]]
ucl     = []
ucl_mpc = []
ucl_cbf = []

# Initialize Controller
if linearMPC == 1:
	ftocp = FTOCP(N, Alin, Blin, Q, Qf, R, K, Fx, bx, x_tightening, Fu, bu, u_tightening)
else:
	ftocp = nonlinearFTOCP(N, dt_nonlinear, bx, x_tightening, bu, u_tightening, Q, R, Qf, fixInput)

cbf_clf = CBF_CLF(vmax, dotphimax, xmax, phimax)

############################# Run main loop #####################################
time = [0]
mpcCounter = int(dt_mpc/dt)


while (time[-1] <= T ):

	# Solve FTOCP
	if mpcCounter == int(dt_mpc/dt):
		mpcCounter = 1
		ftocp.solve(xcl[-1], verbose = False) 
		# Read optimal input
		ut_Mpc = ftocp.mpcInput
		# Reset High Level Planning Model Internal State
		xncl[-1] = ftocp.xPred[:,0]

		print("u MPC: ", ut_Mpc)
	else:
		mpcCounter += 1

	# Solve CLF-CBF QP
	if useLowLevel == 1:
		cbf_clf.solve(xcl[-1], xncl[-1], ut_Mpc, verbose = False)
		u_CBF = cbf_clf.uCBF
	else:
		u_CBF = 0

	# Integreate nominal system
	xnNext = (np.dot(A_int,xncl[-1]) + np.squeeze(np.dot(B_int,ut_Mpc))).tolist()
	xncl.append(xnNext)

	# Compute total input
	if useLowLevel == 1:
		u_tot = ut_Mpc + u_CBF
	else:
		u_tot = ut_Mpc

	# Store inputs
	ucl_mpc.append(ut_Mpc)
	ucl_cbf.append(u_CBF)
	ucl.append(u_tot)

	# Apply input to the system
	xcl.append(simulator.sim(xcl[-1], ucl[-1]))

	time.append(time[-1]+dt)

	if xcl[-1][2] > 1.3:
		break
	# Increasi time counter
	print("Time: ", time[-1], "Closed-Loop: ", xcl[-1], "Input MPC: ",ut_Mpc, "Input CBF: ", u_CBF)



############################# Save closed loop data #####################################
outfile = TemporaryFile()

if linearMPC == 0:
	useLowLevel = 2
	dt_mpc = dt_nonlinear

np.savetxt('storedData/xcl_dt_mpc_'+str(dt_mpc)+'_lowLevelUsed_'+str(useLowLevel)+'.txt', np.round(np.array(xcl), decimals=5).T, fmt='%f' )
np.savetxt('storedData/xncl_dt_mpc_'+str(dt_mpc)+'_lowLevelUsed_'+str(useLowLevel)+'.txt', np.round(np.array(xncl), decimals=5).T, fmt='%f' )
np.savetxt('storedData/ucl_dt_mpc_'+str(dt_mpc)+'_lowLevelUsed_'+str(useLowLevel)+'.txt', np.round(np.array(ucl), decimals=5).T, fmt='%f' )
np.savetxt('storedData/ucl_mpc_dt_mpc_'+str(dt_mpc)+'_lowLevelUsed_'+str(useLowLevel)+'.txt', np.round(np.array(ucl_mpc), decimals=5).T, fmt='%f' )
np.savetxt('storedData/ucl_cbf_dt_mpc_'+str(dt_mpc)+'_lowLevelUsed_'+str(useLowLevel)+'.txt', np.round(np.array(ucl_cbf), decimals=5).T, fmt='%f' )
np.savetxt('storedData/time_'+str(dt_mpc)+'_lowLevelUsed_'+str(useLowLevel)+'.txt', np.round(np.array(time), decimals=5).T, fmt='%f' )
np.save(outfile, xcl)

############################# Plot Results #####################################

xcl = np.array(xcl).T
xncl = np.array(xncl).T
ucl = np.array(ucl).T
ucl_mpc = np.array(ucl_mpc).T
ucl_cbf = np.array(ucl_cbf).T
time = np.array(time)

plt.figure()
plt.subplot(4, 1, 1)	
plt.plot(time, xcl[0,:], '-o', color='C0')#, label="LMPC closed-loop for P = "+str(P))
plt.plot(time, xncl[0,:], '-s', color='C2')#, label="LMPC closed-loop for P = "+str(P))
plt.ylabel('$p$', fontsize=20)
plt.legend()

plt.subplot(4, 1, 2)
plt.plot(time, xcl[1,:], '-o', color='C0')#, label="LMPC closed-loop for P = "+str(P))
plt.plot(time, xncl[1,:], '-s', color='C2')#, label="LMPC closed-loop for P = "+str(P))
plt.ylabel('$v$', fontsize=20)
plt.legend()

plt.subplot(4, 1, 3)
plt.plot(time, xcl[2,:], '-o', color='C0')#, label="LMPC closed-loop for P = "+str(P))
plt.plot(time, xncl[2,:], '-s', color='C2')#, label="LMPC closed-loop for P = "+str(P))
plt.ylabel('$\\theta$', fontsize=20)
plt.legend()

plt.subplot(4, 1, 4)
plt.plot(time, xcl[3,:], '-o', color='C0')#, label="LMPC closed-loop for P = "+str(P))
plt.plot(time, xncl[3,:], '-s', color='C2')#, label="LMPC closed-loop for P = "+str(P))
plt.ylabel('$\omega$', fontsize=20)
plt.xlabel('$\mathrm{time}$', fontsize=20)
plt.legend()

plt.show()

plt.figure()
plt.plot(time[0:ucl.shape[0]], ucl, '-o', color='C0', label="Total input")
plt.plot(time[0:ucl_mpc.shape[0]], ucl_mpc, '-s', color='C3', label="High level control action")
plt.plot(time[0:ucl_cbf.shape[0]], ucl_cbf, '-s', color='C4', label="Low level control action")
plt.xlabel('$\mathrm{time}$', fontsize=20)
plt.legend()

plt.show()
