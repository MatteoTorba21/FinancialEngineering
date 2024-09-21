function optionPrice=EuropeanOptionKICRR(F0,K, KI,B,T,sigma,N)
%European option price with CRR tree
%
%INPUT
% F0:          forward price
% B:           discount factor
% K:           strike
% KI:          barrier
% T:           time-to-maturity
% sigma:       volatility
% N:           number of CRR tree steps
%
%OUTPUT
% optionPrice: price of an European Call with European Barrier via CRR 

deltaT = T/N;  % length of each interval of the tree
u = exp(sigma*sqrt(deltaT)); % up coefficient
d = 1/u;  % down coefficient 
q = (1-d)/(u-d); % up probability

% costruction of the tree
EuropeanKITree = zeros(N+1,N+1);
for ii=0:N
    EuropeanKITree(ii+1,N+1) = max((F0*(u^(N-2*ii))-K) .* (F0*(u^(N-2*ii))>KI),0);
end

for jj=N:-1:1
    for ii=1:jj
        % for each node we take the expectation between the values at its branches
        EuropeanKITree(ii,jj) = q*EuropeanKITree(ii,jj+1)+(1-q)*EuropeanKITree(ii+1,jj+1);
    end
end
optionPrice = B*EuropeanKITree(1,1); % discoounted EuropeanKICRR price

end % function EuropeanOptionKICRR

