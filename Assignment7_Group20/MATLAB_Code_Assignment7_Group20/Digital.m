function [difference] = Digital(data,r,TTM,Notional,digitalPayoff)
% Function to compute the difference between digital option price per Black model
% and digital option price considering implied volatility curve
%
% INPUT:
%           data:       data struct given 
%           r:          zero interest rate at TTM
%           TTM:        Time To Maturity
%           Notional:   Notional of the digital option
% OUTPUT:    
%           difference: difference between the two digital prices

% extract strikes from dataset
strikes = data.strikes;
% extract volatilities from dataset
vol = data.surface;
% extract strike 
K = data.K;
% extract Stock value in 0
S0 = data.reference;
% extract the dvidend yield from dataset
d = data.dividend;

% extract index of the first strike bigger than K
t = find(strikes > K,1);

% compute the slope of the tangent at the volatility curve passing per K
slope = (vol(t) - vol(t-1))/(strikes(t) - strikes(t-1));

% interpolate to get the volatility relative to K
sigma = interp1([strikes(t-1),strikes(t)],[vol(t-1),vol(t)],K);
% compute d1
d1 = (r - d + sigma^2/2)*sqrt(TTM)/sigma;
d1 = (log(S0/K)+(r - d + sigma^2/2)*TTM)/(sigma*sqrt(TTM));
% compute vega
vega = S0*exp(-d*TTM)*sqrt(TTM)*normpdf(d1);
% compute the difference
difference = -vega*slope*Notional*digitalPayoff;
end