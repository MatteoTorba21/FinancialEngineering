function NPV = payerSwapNPV(dates, quotedSwaps, rates, swapMaturities, swapRates)
% PAYERSWAPNPV 
% Compute the NPVs of a fixed leg payer swaps, given maturities and swap
% rates
% INPUTS:
% dates:        Struct of dates of the financial instruments quoted in the market
% quotedSwaps:  Vector of the quoted swap expiries
% rates:        Struct of rates of the financial instruments quoted in the market
% swapMaturities:      Vector of swaps maturities
% swapRates:    Vector of swap rates
% OUTPUT:
% NPV:            NPVs of the swaps


% Yearly swap expiries up to 50 years:
FullSwapExpiries = [1:50]';
% Yearly swap dates up to 50 years:
FullSwDates = findSwapDates(dates.settlement, FullSwapExpiries, eurCalendar);
% Bootstrap discounts up to 50 years
[dates_bootstrap, discounts_bootstrap] = interpolateAndLaunchBootsrap(dates,rates,quotedSwaps,FullSwapExpiries);
% swap date convention
SwapsConvention = 6;
% NPV inizialization
NPV = zeros(size(swapMaturities));
% yearly discounts computation
yearlyDiscounts = InterpDFviaRates(dates_bootstrap,discounts_bootstrap,FullSwDates');
% deltas compèutation
deltas = yearfrac([dates.settlement; FullSwDates(1:end-1)'], FullSwDates(1:end)', SwapsConvention);
% Basis point value vector (discount times respective delta)
BPVvector = deltas.*yearlyDiscounts;
% Swaps cycle
for ii = 1:length(NPV)
    % NPV computation for each swap
    NPV(ii) = 1 - yearlyDiscounts(swapMaturities(ii)) - sum(BPVvector(1:swapMaturities(ii)))*swapRates(ii);
end

end

