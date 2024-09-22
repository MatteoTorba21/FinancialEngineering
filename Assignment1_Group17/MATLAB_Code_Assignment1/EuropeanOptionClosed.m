function optionPrice=EuropeanOptionClosed(F0,K,B,T,sigma,flag)
% European option price with Closed formula
%
%INPUT
% F0:          forward price
% B:           discount factor
% K:           strike
% T:           time-to-maturity
% sigma:       volatility
% flag:        1 call, -1 put
%
%OUTPUT
% optionPrice: Call/Put price with Black Formula

% blkprice: built-in Matlab function to compute black price
[call, put] = blkprice(F0, K, 0, T, sigma); % 

if flag == 1 
   optionPrice= B*call; % discounted Call price
else
   optionPrice= B*put; % discounted Put price
end

end % function EuropeanOptionClosed