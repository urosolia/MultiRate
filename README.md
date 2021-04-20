# Multi-Rate Control Design Leveraging Control Barrier Functions and Model Predictive Control Policies

### Prerequisites

The python packeges needed for running the code can be installed using pip

```
pip install cvxopt
```

To run the Matlab code is required the MPT toolbox [link](https://www.mpt3.org/Main/Installation)

## Description

### python folder
Run the main.py to generate the data contained in the 'storedData' folder.

You can modify the following paramters:
1) Line 20: choose if you want to run the first or second example from the paper.
2) Line 21: choose if you want to use the low level high freqeuncy feedback.
3) Line 22: choose if you want to use linear or nonlinaer MPC for high-level control

To plot the results run the file plot.py

### matlab
Run main.m if you want to

1) to update the QP python parameters stored in the folders "python/parameters_dt_mpc_0.1" and "python/parameters_dt_mpc_0.1" (not required to execute closed-loop pythong simulations)
2) Run the closed-loop simulation in matlab