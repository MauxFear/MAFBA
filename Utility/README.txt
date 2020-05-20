MAFBA (Membrane constrained Allocation Flux Balance Analysis) 
Some requirements to run MAFBA are to have a CPLEX solver added to the MATLAB Path.
Also is required the Cobra tool box package to read the model file (https://github.com/opencobra/cobratoolbox) 

The software version tested was:
MATLAB R2017b
IBM ILOG CPLEX Optimization Studio V12.8.0
COBRAtoolbox version 3
Model E coli: iJR904


Test:
Run the file setUpModel_iJR904_MAFBA_zy_v5.m change the model file path 
 ¨model_iJR_M3 = readCbModel('D:\Code\MATLAB\Projects\CAFBA\Models\Ec_iJR904_flux1.xml');¨

We can set a different weight for each protein sector or for all of them. 

Run MAFBA using the cplex solver. 
We can change the gamma ratio given the ¨gamma¨ value
sol_iJR_m3 = MAFBA_OptCbModel_cplex_v4(model_iJR_M3, 'gamma', 1, 'cSense_m', 'U');

With the file C_lim_plot_iJR_m_MAFBA_v3.m we can simulate carbon limitation changing the value of WCc. 

With the file Z_C_lim_plot_iJR_m_v3.m we can simulate the protein overexpression under carbon limitation conditions changing the fraction of a cytosolic protein Z or a membrane protein Y. 