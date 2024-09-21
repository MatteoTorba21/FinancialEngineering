function FLDates = findFloatingLegDates(settlementDate,Maturity, Holidays)
% Function to determine the dates for the floating leg payments every 3 months
%
% INPUTS:
% settlementDate:       Settlement date
% Maturity:             Time to maturity (in years) of the contract
% Holidays:             Vector of holidays dates
%
% OUTPUT:
% FLDates:              Vector of the floating leg dates 
%

% Create the the vector of floating leg payment dates (every 3 month)
FLDates = (arrayfun(@(x) addtodate(settlementDate, 3*x, 'month'), 1:4*Maturity)');

% Verify for each floating leg payment date if it is a business day
nonBusinessDays = ~isbusday(FLDates,Holidays);

% If the date obtained is not a business day, the next business day is taken, according to the 'modified follow' convenction
FLDates(nonBusinessDays) = busdate(FLDates(nonBusinessDays), 'modifiedfollow', Holidays);

end % function findFloatingLegDates