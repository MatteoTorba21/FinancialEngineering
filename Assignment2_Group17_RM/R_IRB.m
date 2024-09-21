function [R] = R_IRB(PD)
% Function to compute IRB correlation (large corporate-sovereign) for
% exposure with given PD
% INPUTS:
% PD:       Probability of default
% OUTPUTS:
% R:        IRB correlation

R_min = 0.12;
R_max = 0.24;
% k-factor:
k = 50;
% Regulatory correlation function:
R = R_min*(1-exp(-k*PD))/(1-exp(-k)) + R_max*(exp(-k*PD))/(1-exp(-k));


end % functionR_IRB