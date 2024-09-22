function [asset1,asset2,asset1_minus,asset2_minus] = simulateAsset(S1,S2,sigma1, sigma2,r,timefrac,div1,div2,correlation)
% Function which simulates the dynamics of two correlated asset in the Black and
% Scholes framework using the plus and the minus sign at the exponential to
% apply the antithetic variables technique
% INPUTS:
% S1:           First asset value at the initial time
% S2:           Second asset value at the initial time
% sigma1:       Volatility of the first asset
% sigma2:       Volatility of the second asset
% r:            Vector of the forward rates        
% timefrac:     Vector of time fractions
% div1:         Dividend yield of the first asset
% div2:         Dividend yield of the second asset
% correlation:  Correlation between the Brownian Motion of the dynamics
% OUPUTS:
% asset1:       Matrix of the simulations of the first asset: in each row
%               there is the simulation of the asset from the time 0 to the
%               time 4
% asset2:       Matrix of the simulations of the second asset:the i-th row
%               corresponds to the i-th simulation of the asset from time
%               t = 0 to time t = 4years
% asset1_minus: Matrix of the simulations of the first asset using the
%               minus sign to apply the antithetic variables technique
% asset2_minus: Matrix of the simulations of the second asset using the
%               minus sign to apply the antithetic variables technique

% Set the seed
rng(20)
% Number of simulations
Nsim = 1e5;
% Stock inizialization
asset1 = zeros(Nsim,length(timefrac)+1);
asset2 = zeros(Nsim,length(timefrac)+1);
asset1_minus = zeros(Nsim,length(timefrac)+1);
asset2_minus = zeros(Nsim,length(timefrac)+1);
asset1(:,1) = S1;
asset2(:,1) = S2;
asset1_minus(:,1) = S1;
asset2_minus(:,1) = S2;
% Yearly cicle for stocks simulation
for ii = 1:length(timefrac)
    % Correlated standard normal vector sampling for each simulation
    Z = mvnrnd([0;0],[1,correlation;correlation,1], Nsim);
    % Stock update
    asset1(:,ii+1) = asset1(:,ii).*exp((r(ii) - div1 - sigma1^2/2)*timefrac(ii) + sqrt(timefrac(ii)) * sigma1 * Z(:,1));
    asset2(:,ii+1) = asset2(:,ii).*exp((r(ii) - div2 - sigma2^2/2)*timefrac(ii) + sqrt(timefrac(ii)) * sigma2 * Z(:,2));
    asset1_minus(:,ii+1) = asset1_minus(:,ii).*exp((r(ii) - div1 - sigma1^2/2)*timefrac(ii) - sqrt(timefrac(ii)) * sigma1 * Z(:,1));
    asset2_minus(:,ii+1) = asset2_minus(:,ii).*exp((r(ii) - div2 - sigma2^2/2)*timefrac(ii) - sqrt(timefrac(ii)) * sigma2 * Z(:,2));
end

end