function Price = capPrice(discounts,dates,deltas,K,sigma,YTM)
% Function to compute the price of a Cap option via the Bachelier's formula
% INPUTS:
% discounts:    Vector of the discount factors interpolated every 3 months
%               from the bootstrap results
% dates:        Vector of caplet dates
% deltas:       Vector of year fractions
% K:            Vector of strikes
% sigma:        Volatility vector for the maturity indicated in YTM
% YTM:          Years to maturity of the Caps at the start date
% OUTPUTS:
% capPrice:     Vector of Cap prices computed via Bachelier's formula

% Initialization of the results vector:
Price = zeros(length(K),1);
% Number of caplets:
ncapl = 3+4*(YTM-1);
% yearfrac Act/360:
DayCount = 2;
% yearfrac Act/365:
ACT365 = 3;
% Select the needed year fractions:
deltas = deltas(2:ncapl+1);
% Compute the forward discounts:
fwd_discounts = discounts(2:end)./discounts(1:end-1);
fwd_discounts = fwd_discounts(2:ncapl+1); % select only the needed forward discounts
% Compute the forward Libor:
L = (1./fwd_discounts -1)./deltas;
for ii=1:length(K)
    % Compute dn:
    d = (L-K(ii))./(sigma(ii)*1e-4.*sqrt(yearfrac(dates(1),dates(2:ncapl+1),3)));
    % Computhe the caplets' vector:
    caplet = discounts(3:ncapl+2).*deltas.*...
        ((L-K(ii)).*normcdf(d)+sigma(ii)*1e-4.*sqrt(yearfrac(dates(1),dates(2:ncapl+1),ACT365)).*normpdf(d));
    % Compute the Cap price for the ii-th strike K:
    Price(ii) = sum(caplet);
end
end