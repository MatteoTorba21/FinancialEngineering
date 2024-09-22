function zRates = zeroRates(dates, discounts)
% Obtain zero Rates from the discounts and their dates
%
% INPUT:
% dates:        vector of the dates for which we have the discount
% discounts:    vector of the discounts from the settlement day to the
%               date in the same position in the vector dates
% OUTPUT:
% zRates:       zero rates at the dates given in the vector dates
%

% Basis of year fractions equal to 3 because time intervals for zero rates are measured with an Act/365 yearfrac
ACT_365 = 3;

% Compute the vector of the year fractions needed
delta = yearfrac(dates(1)*ones(length(dates)-1,1),dates(2:end),ACT_365); 

% Inversed discount formula to compute the rates from the discounts
zRates = -log(discounts(2:end))./delta;

end % function zeroRates

