function [dates_bootstrap, discounts_bootstrap] = interpolateAndLaunchBootsrap(dates,rates,SwapExpiries,FullSwapExpiries)
% Function to build a complete curve of swap rates via spline interpolation
% on the mid rates of quoted swap rates and to launch the bootstrap to
% obtain the whole discounts curve.
%
% INPUTS:
% dates:                Struct of dates from market data         
% rates:                Struct of rates from market data
% SwapExpiries:         Quoted swap expires (in years)
% FullSwapExpiries:     Swap expiries (in years) up to 50 years
%
% OUTPUTS:
% dates_bootstrap:      Vector of dates obtained via bootstrap   
% discounts_bootstrap:  Vector of discounts obtained via bootstrap 

% Dates of quoted swaps:
SwDates = findSwapDates(dates.settlement, SwapExpiries, eurCalendar);
% Dates of swaps up to 50 years:
FullSwDates = findSwapDates(dates.settlement, FullSwapExpiries, eurCalendar);
% Mid swap rates computed as the mean between bid and ask:
midSwapRates = mean(rates.swaps,2);
% Spline interpolation of the midrates:
FullSwapRates = interp1(SwDates,midSwapRates,FullSwDates,'spline')';
% Update dates and rates structs:
dates.swaps = FullSwDates';
rates.swaps = FullSwapRates;

% Bootstrap dates and discounts up to 50 years:
[dates_bootstrap, discounts_bootstrap] = bootstrap(dates, rates);

end