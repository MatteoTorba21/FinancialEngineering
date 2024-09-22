function optionPrice=EuropeanOptionKIMC(F0,K,KI,B,T,sigma,N) 
%European Barrier option price with MonteCarlo method
%
%INPUT
% F0:          forward price
% B:           discount factor
% K:           strike
% KI:          barrier
% T:           time-to-maturity
% sigma:       volatility
% N:           number of simulations
%
%OUTPUT
% optionPrice: price of an European Call with European Barrier via Monte-Carlo method

g = randn(1,N);
callPayoff = max((F0*exp(-(sigma^2)*T/2+sigma*sqrt(T)*g)-K) .* ...
    (F0*exp(-(sigma^2)*T/2+sigma*sqrt(T)*g)> KI ),0);
optionPrice = B*mean(callPayoff);

end %function EuropeanOptionKIMC 

