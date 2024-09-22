function df = InterpDFviaRates(dates,dfs,date)
% Interpolate/Extrapolate the discount factors from the settlement date to date
% by passing to the linear interpolation/flat extrapolation of the zero rate
%
% INPUT:
% dates:     vector of the dates for which the discount factors are known
% dfs:       vector of the discount factors from the settlement day to the corresponding date in the same 
%            position in the vector dates
% date:      vector of dates for which the interpolation/extrapolation of the discounts must be done
%
% OUTPUT:
% df:        vector of interpolated/extrapolated discount factors at the dates given in the vector date 
%

% Compute the rates from the known discounts
rates = zeroRates(dates, dfs);

% Interpolate/extrapolate the rates
rate = InterpZeroRate(dates(2:end,1), rates, date);

% Pass from the interpolated/extrapolated rates to the discounts:
df = DiscountsFromRates([dates(1); date],rate);

end % function InterpDFviaRates

