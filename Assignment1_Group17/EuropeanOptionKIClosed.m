function optionPrice = EuropeanOptionKIClosed(F0,K, KI,B,T,sigma)
%European option price with exact formula. The value of the option is
%obtained by replicating strategy: it uses a plain vanilla having KI, the
%barrier value, as strike and a (KI - K) digital options with KI as the
%strike.
%
%INPUT
% F0:           forward price
% B:            discount factor
% K:            strike
% KI:           barrier
% T:            time-to-maturity
% sigma:        volatility
%  
%OUTPUT
% option price: price of an European Call with European Barrier via Closed Formula

flag = 1;
optionPrice = (KI - K) * DigitalOptionClosed(F0,KI,B,T,sigma) + ...
    EuropeanOptionClosed(F0,KI,B,T,sigma,flag);

end % function EuropeanOptionKIClosed