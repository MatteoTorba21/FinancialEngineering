function optionPrice=EuropeanOptionCRR(F0,K,B,T,sigma,N,flag)
% European option price with CRR tree
%
%INPUT
% F0:            forward price
% B:             discount factor
% K:             strike
% T:             time-to-maturity
% sigma:         volatility
% N:             number of intervals in CRR tree
% flag:          1 call, -1 put
%
%OUTPUT
% optionPrice:  Call/Put price with the CRR tree method

deltaT = T/N;                % length of each interval of the tree
u = exp(sigma*sqrt(deltaT)); % up coefficient
d = 1/u;                     % down coefficient 
q = (1-d)/(u-d);             % up probability

if flag==1 % Call case
    % costruction of the tree
    callTree = zeros(N+1,N+1);
    for ii=0:N
        % in this relation we use the fact that d=1/u
        callTree(ii+1,N+1) = max(F0*(u^(N-2*ii))-K,0);
    end
    for jj=N:-1:1
        for ii=1:jj
            % for each node we take the expectation between the values at its branches
            callTree(ii,jj) = q*callTree(ii,jj+1)+(1-q)*callTree(ii+1,jj+1);
        end
    end
    optionPrice = B*callTree(1,1); % discounted Call price
end

if flag==-1 % Put case 
    % costruction of the tree
    putTree = zeros(N+1,N+1);
    for ii=0:N
        % in this relation we use the fact that d=1/u
        putTree(ii+1,N+1) = max(K-F0*(u^(N-2*ii)),0);
    end
    for jj=N:-1:1
        for ii=1:jj
            % for each node we take the expectation between the values at its branches
            putTree(ii,jj) = q*putTree(ii,jj+1)+(1-q)*putTree(ii+1,jj+1);
        end
    end
    optionPrice = B*putTree(1,1); % discouted Put price
end

end % function EuropeanOptionCRR