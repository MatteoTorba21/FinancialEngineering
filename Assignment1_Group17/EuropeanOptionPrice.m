function OptionPrice=EuropeanOptionPrice(F0,K,B,T,sigma,pricingMode,N,flag)
% Option Price with different pricing methods
%
% INPUT:
% F0:           forward price
% B:            discount factor
% K:            strike
% T:            time-to-maturity
% sigma:        volatility
% pricingMode:  1 ClosedFormula, 2 CRR, 3 Monte Carlo
% N:            either number of time steps (knots for CRR tree)
%               or number of simulations in MC   
% flag:         1 call, -1 put

if (nargin < 7)
 N = 10000; % Default: N
end 

if (nargin < 8)
 flag = 1; % Default: Call price
end 

switch (pricingMode)
    % For each case we computed the execution time via the functions tic & toc
    case 1  % Closed Formula
        % Start timer
        tic;
        OptionPrice = EuropeanOptionClosed(F0,K,B,T,sigma,flag);
        % Stop timer
        elapsed_time = toc;
        % Print the elapsed time
        fprintf('Execution time to compute the price via the closed formula: %.6f seconds\n', elapsed_time);
    
    case 2  % CRR
        % Start timer
        tic;
        OptionPrice = EuropeanOptionCRR(F0,K,B,T,sigma,N,flag);
        % Stop timer
        elapsed_time = toc;
        % Print the elapsed time
        fprintf('Execution time to compute the price via CRR: %.6f seconds\n', elapsed_time);

    case 3  % Monte Carlo
        % Start timer
        tic;
        OptionPrice = EuropeanOptionMC(F0,K,B,T,sigma,N,flag);
        % Stop timer
        elapsed_time = toc;
        % Print the elapsed time
        fprintf('Execution time to compute the price via MC: %.6f seconds\n', elapsed_time);
    otherwise
end
return
end % function EuropeanOptionPrice