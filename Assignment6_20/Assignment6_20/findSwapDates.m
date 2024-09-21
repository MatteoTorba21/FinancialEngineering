function SwDates = findSwapDates(settlementDate, SwapExpiries, Holidays)
% Function to determine the swap dates for the bootstrap
%
% INPUTS:
% settlementDate:       Settlement date
% SwapExpiries:         Vector with the expiries of the swaps
% Holidays:             Vector of holidays dates of the Italian Stock Exchange
%
% OUTPUT:
% SwDates:              Vector of the swap dates 
%

% Create the the vector of swap dates
SwDates = (arrayfun(@(x) addtodate(settlementDate, x, 'year'), SwapExpiries)');

% Verify for each date if it is a business day
nonBusinessDays = ~isbusday(SwDates,Holidays);

% If the date obtained is not a business day, the next business day is taken, according to the 'modifiedfollow' convenction
SwDates(nonBusinessDays) = busdate(SwDates(nonBusinessDays), 'modifiedfollow', Holidays);

end % function findSwapDates