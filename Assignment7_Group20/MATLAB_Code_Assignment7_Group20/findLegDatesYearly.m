function FLDates = findLegDatesYearly(settlementDate,Maturity, Holidays)
% Function to determine the dates for the floating leg payments yearly
%
% INPUTS:
% settlementDate:       Settlement date
% Maturity:             Time to maturity (in years)
% Holidays:             Vector of holidays dates
%
% OUTPUT:
% FLDates:              Vector of the yearly dates 
%

% Create the the vector of yearly dates
FLDates = (arrayfun(@(x) addtodate(settlementDate, x, 'year'), 1:Maturity)');

% Verify for each yearly date if it is a business day
nonBusinessDays = ~isbusday(FLDates,Holidays);

% If the date obtained is not a business day, the next business day is taken, according to the 'modified follow' convenction
FLDates(nonBusinessDays) = busdate(FLDates(nonBusinessDays), 'modifiedfollow', Holidays);

end % function findFloatingLegDates