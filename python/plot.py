import numpy as np 
import matplotlib.pyplot as plt
import copy
import matplotlib
from matplotlib import rc
import pdb
from pylab import rcParams
from matplotlib.ticker import StrMethodFormatter

matplotlib.rc('xtick', labelsize=14) 
matplotlib.rc('ytick', labelsize=14) 

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

# Example to plot
Example = 1

# Parameter to pick
if Example == 1:
	dt_mpc = 0.5
	dt_nonlinear = [0.01, 0.05]
	plotNonlinearMPC = 1

elif Example == 2:
	dt_mpc = 0.1
	dt_nonlinear = [dt_mpc]
	plotNonlinearMPC = 1
	fixInput = 0


# Load data
useLowLevel = 1
xcl_withBarrier     = np.loadtxt('storedData/xcl_dt_mpc_'+str(dt_mpc)+'_lowLevelUsed_'+str(useLowLevel)+'.txt')
xncl_withBarrier    = np.loadtxt('storedData/xncl_dt_mpc_'+str(dt_mpc)+'_lowLevelUsed_'+str(useLowLevel)+'.txt')
ucl_withBarrier     = np.loadtxt('storedData/ucl_dt_mpc_'+str(dt_mpc)+'_lowLevelUsed_'+str(useLowLevel)+'.txt')
ucl_mpc_withBarrier = np.loadtxt('storedData/ucl_mpc_dt_mpc_'+str(dt_mpc)+'_lowLevelUsed_'+str(useLowLevel)+'.txt')
ucl_cbf_withBarrier = np.loadtxt('storedData/ucl_cbf_dt_mpc_'+str(dt_mpc)+'_lowLevelUsed_'+str(useLowLevel)+'.txt')
time_withBarrier    = np.loadtxt('storedData/time_'+str(dt_mpc)+'_lowLevelUsed_'+str(useLowLevel)+'.txt')

useLowLevel = 0
xcl_withoutBarrier     = np.loadtxt('storedData/xcl_dt_mpc_'+str(dt_mpc)+'_lowLevelUsed_'+str(useLowLevel)+'.txt')
xncl_withoutBarrier    = np.loadtxt('storedData/xncl_dt_mpc_'+str(dt_mpc)+'_lowLevelUsed_'+str(useLowLevel)+'.txt')
ucl_withoutBarrier     = np.loadtxt('storedData/ucl_dt_mpc_'+str(dt_mpc)+'_lowLevelUsed_'+str(useLowLevel)+'.txt')
ucl_mpc_withoutBarrier = np.loadtxt('storedData/ucl_mpc_dt_mpc_'+str(dt_mpc)+'_lowLevelUsed_'+str(useLowLevel)+'.txt')
ucl_cbf_withoutBarrier = np.loadtxt('storedData/ucl_cbf_dt_mpc_'+str(dt_mpc)+'_lowLevelUsed_'+str(useLowLevel)+'.txt')
time_withoutBarrier    = np.loadtxt('storedData/time_'+str(dt_mpc)+'_lowLevelUsed_'+str(useLowLevel)+'.txt')

if plotNonlinearMPC == 1:
	useLowLevel = 2
	xcl_nonLinear          = []
	xncl_nonLinear         = []
	ucl_nonLinear         = []
	ucl_mpc_nonLinear      = []
	ucl_cbf_nonLinear      = []
	time_nonLinear         = []

	for i in range(0,len(dt_nonlinear)):
		xcl_nonLinear.append(np.loadtxt('storedData/xcl_dt_mpc_'+str(dt_nonlinear[i])+'_lowLevelUsed_'+str(useLowLevel)+'.txt'))
		xncl_nonLinear.append(np.loadtxt('storedData/xncl_dt_mpc_'+str(dt_nonlinear[i])+'_lowLevelUsed_'+str(useLowLevel)+'.txt'))
		ucl_nonLinear.append(np.loadtxt('storedData/ucl_dt_mpc_'+str(dt_nonlinear[i])+'_lowLevelUsed_'+str(useLowLevel)+'.txt'))
		ucl_mpc_nonLinear.append(np.loadtxt('storedData/ucl_mpc_dt_mpc_'+str(dt_nonlinear[i])+'_lowLevelUsed_'+str(useLowLevel)+'.txt'))
		ucl_cbf_nonLinear.append(np.loadtxt('storedData/ucl_cbf_dt_mpc_'+str(dt_nonlinear[i])+'_lowLevelUsed_'+str(useLowLevel)+'.txt'))
		time_nonLinear.append(np.loadtxt('storedData/time_'+str(dt_nonlinear[i])+'_lowLevelUsed_'+str(useLowLevel)+'.txt'))

colors = ['C1', 'C2', 'C4']


rcParams['figure.figsize'] = 8, 9
matplotlib.rc('xtick', labelsize=10) 
matplotlib.rc('ytick', labelsize=10) 
# Plot Data
plt.figure()
plt.subplot(4, 1, 1)	
plt.plot([time_withoutBarrier[0],time_withoutBarrier[-1]], [0,0], '-', color='k', label="Goal")
plt.plot(time_withBarrier, xcl_withBarrier[0,:], '-', color='C0', label="Proposed strategy")
plt.plot(time_withoutBarrier, xcl_withoutBarrier[0,:], '-', color='C3', label="Linear MPC")
if plotNonlinearMPC == 1:
	for i in range(0,len(dt_nonlinear)):
		plt.plot(time_nonLinear[i], xcl_nonLinear[i][0,:], '--', color = colors[i], label="Nonlinaer MPC discretized at "+str(int(1/dt_nonlinear[i]))+"Hz" )

plt.ylabel('$p\mathrm{[m]}$', fontsize=20)
plt.legend(loc='upper right',fontsize=12)
plt.xlim(time_withBarrier[0],time_withBarrier[-1])
plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.1f}')) # 2 decimal places
plt.gca().xaxis.set_major_formatter(StrMethodFormatter('{x:,.1f}')) # 2 decimal places

plt.subplot(4, 1, 2)
plt.plot([time_withoutBarrier[0],time_withoutBarrier[-1]], [0,0], '-', color='k')
plt.plot(time_withBarrier, xcl_withBarrier[1,:], '-', color='C0')#, label="LMPC closed-loop for P = "+str(P))
plt.plot(time_withoutBarrier, xcl_withoutBarrier[1,:], '-', color='C3')#, label="LMPC closed-loop for P = "+str(P))
if plotNonlinearMPC == 1:
	for i in range(0,len(dt_nonlinear)):
		plt.plot(time_nonLinear[i], xcl_nonLinear[i][1,:], '--', color = colors[i])#, label="LMPC closed-loop for P = "+str(P))

plt.ylabel('$v\mathrm{[m/s]}$', fontsize=20)
plt.xlim(time_withBarrier[0],time_withBarrier[-1])
plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.1f}')) # 2 decimal places
plt.gca().xaxis.set_major_formatter(StrMethodFormatter('{x:,.1f}')) # 2 decimal places

plt.subplot(4, 1, 3)
plt.plot([time_withoutBarrier[0],time_withoutBarrier[-1]], [0,0], '-', color='k')
plt.plot(time_withBarrier, xcl_withBarrier[2,:], '-', color='C0')#, label="LMPC closed-loop for P = "+str(P))
plt.plot(time_withoutBarrier, xcl_withoutBarrier[2,:], '-', color='C3')#, label="LMPC closed-loop for P = "+str(P))
if plotNonlinearMPC == 1:
	for i in range(0,len(dt_nonlinear)):
		plt.plot(time_nonLinear[i], xcl_nonLinear[i][2,:], '--', color = colors[i])#, label="LMPC closed-loop for P = "+str(P))

plt.ylabel('$\\theta\mathrm{[rad]}$', fontsize=20)
plt.xlim(time_withBarrier[0],time_withBarrier[-1])
plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.1f}')) # 2 decimal places
plt.gca().xaxis.set_major_formatter(StrMethodFormatter('{x:,.1f}')) # 2 decimal places

plt.subplot(4, 1, 4)
plt.plot([time_withoutBarrier[0],time_withoutBarrier[-1]], [0,0], '-', color='k')
plt.plot(time_withBarrier, xcl_withBarrier[3,:], '-', color='C0')#, label="LMPC closed-loop for P = "+str(P))
plt.plot(time_withoutBarrier, xcl_withoutBarrier[3,:], '-', color='C3')#, label="LMPC closed-loop for P = "+str(P))
if plotNonlinearMPC == 1:
	for i in range(0,len(dt_nonlinear)):
		plt.plot(time_nonLinear[i], xcl_nonLinear[i][3,:], '--', color = colors[i])#, label="LMPC closed-loop for P = "+str(P))

plt.ylabel('$\omega\mathrm{[m/s]}$', fontsize=20)
plt.xlabel('$\mathrm{Time~[s]}$', fontsize=20)
plt.xlim(time_withBarrier[0],time_withBarrier[-1])
plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.1f}')) # 2 decimal places
plt.gca().xaxis.set_major_formatter(StrMethodFormatter('{x:,.1f}')) # 2 decimal places


rcParams['figure.figsize'] = 8, 7
matplotlib.rc('xtick', labelsize=14) 
matplotlib.rc('ytick', labelsize=14) 

plt.figure()
plt.plot(time_withBarrier, xcl_withBarrier[2,:], '-', color='C0', label="Proposed strategy",linewidth=2)
plt.plot(time_withoutBarrier, xcl_withoutBarrier[2,:], '-', color='C3', label="Linear MPC",linewidth=2)
if plotNonlinearMPC == 1:
	for i in range(0,len(dt_nonlinear)):
		plt.plot(time_nonLinear[i], xcl_nonLinear[i][2,:], '-', color='k', label="Nonlinear MPC",linewidth=2)

if Example == 2:
	plt.plot([time_withBarrier[0],time_withBarrier[-1]], [0.78, 0.78], '--', color='k', label="Constraint",linewidth=2)
plt.ylim(-0.7,1.5)
plt.yticks([-0.5,0.0,0.78,1.5])

plt.ylabel('$\\theta \mathrm{[rad]}$', fontsize=24)
plt.xlabel('$\mathrm{Time~[s]}$', fontsize=24)
plt.legend(fontsize=18)
plt.xlim(time_withBarrier[0],time_withBarrier[-1])
plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.1f}')) # 2 decimal places
plt.gca().xaxis.set_major_formatter(StrMethodFormatter('{x:,.1f}')) # 2 decimal places


plt.figure()
plt.plot(time_withBarrier[0:ucl_mpc_withBarrier.shape[0]], ucl_mpc_withBarrier, '-', color='C3', label="High level control action $v$",linewidth=2)
plt.plot(time_withBarrier[0:ucl_cbf_withBarrier.shape[0]], ucl_cbf_withBarrier, '-', color='C0', label="Low level control action $u$",linewidth=2)
plt.plot(time_withBarrier[0:ucl_withBarrier.shape[0]], ucl_withBarrier, '-', color='k', label="Total input $u+v$",linewidth=2)
plt.xlabel('$\mathrm{Time~[s]}$', fontsize=24)
plt.ylabel('$\mathrm{Voltage~[V]}$', fontsize=24)
plt.legend(fontsize=18)
plt.xlim(time_withBarrier[0],time_withBarrier[-1])
plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.1f}')) # 2 decimal places
plt.gca().xaxis.set_major_formatter(StrMethodFormatter('{x:,.1f}')) # 2 decimal places


plt.figure()
plt.plot(time_withBarrier[0:ucl_withBarrier.shape[0]], ucl_withBarrier, '-', color='C0', label="Total input with Barrier")
plt.plot(time_withoutBarrier[0:ucl_withoutBarrier.shape[0]], ucl_withoutBarrier, '-', color='C3', label="Total input withoutBarrier")
plt.xlabel('$\mathrm{time}$', fontsize=24)
plt.legend(fontsize=18)
plt.xlim(time_withBarrier[0],time_withBarrier[-1])
plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.1f}')) # 2 decimal places
plt.gca().xaxis.set_major_formatter(StrMethodFormatter('{x:,.1f}')) # 2 decimal places


plt.figure()
plt.plot(time_withBarrier[0:ucl_withBarrier.shape[0]], ucl_withBarrier, '-', color='C0', label="Total input proposed strategy")
if plotNonlinearMPC == 1:
	for i in range(0,len(dt_nonlinear)):
		plt.plot(time_nonLinear[i][0:ucl_nonLinear[i].shape[0]], ucl_nonLinear[i], '-', color='C3', label="Total input nonlinear MPC")
plt.xlabel('$\mathrm{time}$', fontsize=24)
plt.legend(fontsize=18)
plt.xlim(time_withBarrier[0],time_withBarrier[-1])
plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.1f}')) # 2 decimal places
plt.gca().xaxis.set_major_formatter(StrMethodFormatter('{x:,.1f}')) # 2 decimal places

plt.show()