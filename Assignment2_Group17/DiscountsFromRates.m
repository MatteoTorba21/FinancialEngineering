function discounts = DiscountsFromRates(dates,rates)
% Obtain discount factors from the zero rates
%
%INPUT:
% dates:        vector of the dates for which we have the zero rates
% rates:        vector of the zero rates
%OUTPUT:
% discounts:    zero rates at the dates given in the vector dates

% basis of year fractions equal to 3 because time intervals for zero rates are measured with an Act/365 yearfrac
ACT_365 = 3;

% compute the vector of the year fractions needed
delta = yearfrac(dates(1)*ones(length(dates)-1,1),dates(2:end),ACT_365); 

% inversed discount formula to compute the rates from the discounts
discounts = exp(-delta.*rates);

end % function DiscountsFromRates

