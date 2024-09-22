function asset = simulateAsset(S0,sigma,r,timefrac,div,Nsim)
% Function which simulates the dynamics of an asset in the Black and
% Scholes framework using the plus and the minus sign at the exponential to
% apply the antithetic variables technique
% 
% INPUTS:
% S0:           First asset value at the initial time
% sigma:        Volatility of the asset
% r:            Forward rate        
% timefrac:     Time fractions
% div:          Dividend yield of the asset
% Nsim:         Number of simulations
% 
% OUPUTS:
% asset:        Vector of the simulations of the asset using the antithetic
%               variables technique

% Set the seed
rng(20)
% Stock inizialization
asset = zeros(Nsim,1);
% Yearly stock simulation
% Standard normal vector sampling
Z = randn(Nsim/2);
% Stock update
asset(1:Nsim/2) = S0.*exp((r - div - sigma^2/2)*timefrac + sqrt(timefrac) * sigma * Z);
asset(Nsim/2+1:end) = S0.*exp((r - div - sigma^2/2)*timefrac - sqrt(timefrac) * sigma * Z);

end