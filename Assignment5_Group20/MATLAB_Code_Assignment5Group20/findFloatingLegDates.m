function FLDates = findFloatingLegDates(settlementDate,Maturity, Holidays)
% Function to determine the floating leg dates for the payments every 3 months
%
% INPUTS:
% settlementDate:       Settlement date
% Maturity:             Time to maturity (in years) of the contract
% Holidays:             Vector of holidays dates of the Italian Stock Exchange
%
% OUTPUT:
% FLDates:              Vector of the floating leg dates 
%

% Create the the vector of floating leg payment dates (every 3 month)
FLDates = (arrayfun(@(x) addtodate(settlementDate, 3*x, 'month'), 1:4*Maturity)');

% Verify for each floating leg payment date if it is a business day
nonBusinessDays = ~isbusday(FLDates,Holidays);
% [datenum('21-Mar-2008') datenum('24-Mar-2008') datenum('15-Ago-2008') datenum('24-Dec-2008') datenum('25-Dec-2008') datenum('26-Dec-2008') datenum('31-Dec-2008') datenum('01-Jan-2009') datenum('10-Apr-2009') datenum('13-Apr-2009') datenum('01-May-2009') datenum('24-Dec-2009') datenum('25-Dec-2009') datenum('31-Dec-2009') datenum('01-Jan-2010') datenum('02-Apr-2010') datenum('05-Apr-2010') datenum('24-Dec-2010') datenum('31-Dec-2010') datenum('22-Apr-2011') datenum('25-Apr-2011') datenum('15-Ago-2011') datenum('26-Dec-2011')]

% If the date obtained is not a business day, the next business day is taken, according to the 'follow' convenction
FLDates(nonBusinessDays) = busdate(FLDates(nonBusinessDays), 'follow', Holidays);

end % function findFloatingLegDates