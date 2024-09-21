function BermudanOptionPrice = BermudanOptionCRR(F0,K,r,T,EEFreq,sigma,N)
%Bermudan option price with CRR tree
%
%INPUT
% F0:                  forward price
% r:                   zero-rate assumed constant
% K:                   strike
% T:                   time-to-maturity
% EEFreq:              early execise frequency
% sigma:               volatility
% N:                   number of intervals of one period in CRR tree
%
%OUTPUT
% BermudanOptionPrice: Bermudan price computed with CRR tree

deltaT = EEFreq/N;           % length of each interval of the tree
u = exp(sigma*sqrt(deltaT)); % up coefficient 
d = 1/u;                     % down coefficient 
q = (1-d)/(u-d);             % up probability

nPeriods = T/EEFreq;         % number of periods between value date, maturity or EE opportunities
div = r - log(F0)/T;         % dividend yield

% costruction of the tree
callTree = zeros(nPeriods*N+1,nPeriods*N+1);
for ii=0:N*nPeriods  
     callTree(ii+1,nPeriods*N+1) = max(F0*(u^(N*nPeriods-2*ii))-K,0);
end
tic;
for jj=nPeriods * N:-1:1
    for ii=1:jj
        % discounting is performed at each step
        callTree(ii,jj) = (q*callTree(ii,jj+1)+(1-q)*callTree(ii+1,jj+1)) * exp(-r * deltaT);
    end
    if rem(jj,N) == 0 && jj~=nPeriods*N % verified only if we are at one of early exercise dates
        for ii=1:jj
             % from the forward, the undelying is retrived and the intrinsict
             % value is calculated
             intrinsictValue =  max((F0*(u^(jj-2*(ii-1)))) * exp((div - r)*EEFreq)-K,0);
             % take the maximum between IV and CV
             callTree(ii,jj) = max(callTree(ii,jj), intrinsictValue);
        end
    end
end
BermudanOptionPrice = callTree(1,1); % Bermudian price
elapsed_time=toc;
disp(['Execution time to compute the price of the Bermudan via the CRR: ', ...
    num2str(elapsed_time)]);

end % function BermudanOptionCRR





