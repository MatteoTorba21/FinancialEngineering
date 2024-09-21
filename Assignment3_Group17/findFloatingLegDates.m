function FLDates = findFloatingLegDates(settlementDate,Maturity)
% Function to determine the floating leg dates for a float with payments every 3 months
%
% INPUTS:
% settlementDate:       Settlement date
% Maturity:             Time to maturity (in years) of the contract
%
% OUTPUT:
% FLDates:              Vector of the floating leg dates 
%

% Create the the vector of floating leg payment dates (every 3 month)
FLDates = (arrayfun(@(x) addtodate(settlementDate, 3*x, 'month'), 1:4*Maturity)');

% Verify for each floating leg payment date if it is a business day
nonBusinessDays = ~isbusday(FLDates);

% If the date obtained is not a business day, the next business day is taken, according to the 'modified follow' convenction
FLDates(nonBusinessDays) = busdate(FLDates(nonBusinessDays), 'modifiedfollow',[datenum('02-Apr-2010') datenum('24-Dec-2010') datenum('31-Dec-2010') datenum('22-Apr-2011')]);

end % function findFloatingLegDates