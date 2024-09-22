function  Delta = DeltaBucketSwap(dates, quotedSwaps, rates, swapMaturities, swapRates)
%Function that computes the Delta Buckets of a Swap, shocking one quoted
%rate at a time and compute the Delta leaving the rest unchanged
%
% INPUTS:
% dates:         Struct of dates of the financial instruments quoted in the market
% quotedSwaps:   Vector of the quoted swap expiries   
% rates:         Struct of rates of the financial instruments quoted in the market
% swapMaturities:Vector of swaps maturities
% swapRates:     Vector of swap rates
%
%OUTPUTS:
%Delta:          Matrix of Delta Buckets of the swaps, in which the rows
%represent the quoted dates and the coloumns are the different maturities
%of the Swaps 

%Depos Delta bucket initialization
Delta_depos = zeros(4,length(swapRates));
parfor ii = 1:4
    % Rates struct inizialization
    newrates = rates;
    % Increase each bucket iteratively
    newrates.depos(ii,:) = newrates.depos(ii,:) + 1e-4;
    % Store the X value after rate shift
    Delta_depos(ii,:) = (payerSwapNPV(dates, quotedSwaps, newrates, swapMaturities, swapRates)-payerSwapNPV(dates, quotedSwaps, rates, swapMaturities, swapRates))';
end

% Futures Delta bucket inizialization
Delta_futures = zeros(7,length(swapRates));
parfor ii = 1:7
    % Rates struct inizialization
    newrates = rates;
    % Increase each bucket iteratively
    newrates.futures(ii,:) = newrates.futures(ii,:) + 1e-4;
    % Store the X value after rate shift
    Delta_futures(ii,:) = (payerSwapNPV(dates, quotedSwaps, newrates, swapMaturities, swapRates)-payerSwapNPV(dates, quotedSwaps, rates, swapMaturities, swapRates))';
end

% Swaps Delta bucket inizialization
Delta_swaps = zeros(size(rates.swaps,1)-1, length(swapRates));
parfor ii = 2:size(rates.swaps,1)
    % Rates struct inizialization
    newrates = rates;
    % Increase each bucket iteratively
    newrates.swaps(ii,:) = newrates.swaps(ii,:) + 1e-4;
    % Store the X value after rate shift
    Delta_swaps(ii-1,:) = (payerSwapNPV(dates, quotedSwaps, newrates, swapMaturities, swapRates)-payerSwapNPV(dates, quotedSwaps, rates, swapMaturities, swapRates))';

end
%Aggregating the results
Delta = [Delta_depos;Delta_futures;Delta_swaps];
end