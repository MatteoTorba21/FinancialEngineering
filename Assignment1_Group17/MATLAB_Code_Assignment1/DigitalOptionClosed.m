function optionPrice=DigitalOptionClosed(F0,K,B,T,sigma)
%Digital option price with Closed formula
%
%INPUT
% F0:          forward price
% B:           discount factor
% K:           strike
% T:           time-to-maturity
% sigma:       volatility
%
%OUTPUT
% optionPrice: Digital price with Closed Formula

d2 = log(F0/K)/(sigma*sqrt(T)) - sigma*sqrt(T)/2;
optionPrice = B * normcdf(d2);

end % function EuropeanOptionClosed