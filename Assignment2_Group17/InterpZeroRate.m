function rate = InterpZeroRate(dates,rates,date)
% Linear interpolate/flat extrapolate the zero rate at date
%
%INPUT:
% dates:      vector of the dates for which the rates are known
% rates:      vector of the zero rates from the settlement day to the
%             date in the same position in the vector dates
% date:       vector of dates for which the interpolation/extrapolation
%             of the rate has to be done
%OUTPUT:
% rate:       vector of rates at the dates given in the vector date   

% preallocate the memory for the rate vector
rate = zeros(size(date));

% for each element of the vector date, check if extrapolation or interpolation is needed
for ii = 1:length(date)
    if date(ii)<dates(1) % extrapolate before the first date
        r = rates(1);
    elseif date(ii)>dates(end) % extrapolate after the last date
        r = rates(end);
    else
        r = interp1(dates,rates,date(ii),'linear'); % linear interpolation
    end
    % update the vector rate:
    rate(ii,1) = r;
end

end % function InterpZeroRates