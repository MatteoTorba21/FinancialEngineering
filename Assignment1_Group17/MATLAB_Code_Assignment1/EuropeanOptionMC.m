function optionPrice=EuropeanOptionMC(F0,K,B,T,sigma,N,flag) 
%European option price with Monte-Carlo method
%
%INPUT
% F0:           forward price
% B:            discount factor
% K:            strike
% T:            time-to-maturity
% sigma:        volatility
% N:            number of simulation
% flag:         1 call, -1 put
%
%OUTPUT
% optionPrice:  Call/Put price with Monte-Carlo method

g = randn(1,N); % random samples

if flag==1
    % the payoff of a Call is max(F0*exp(-(sigma^2)*T/2+sigma*sqrt(T))-K)
    callPayoff = max(F0*exp(-(sigma^2)*T/2+sigma*sqrt(T)*g)-K,0);
    optionPrice = B*mean(callPayoff); % discounted Call payoff
end

if flag==-1
    % the payoff of a Put is max(K-F0*exp(-(sigma^2)*T/2+sigma*sqrt(T)))
    putPayoff = max(K-F0*exp(-(sigma^2)*T/2+sigma*sqrt(T)*g),0);
    optionPrice = B*mean(putPayoff); % discounted Put payoff
end

end % function EuropeanOptionMC


