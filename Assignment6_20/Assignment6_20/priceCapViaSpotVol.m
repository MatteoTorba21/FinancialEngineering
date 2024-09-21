function Cap = priceCapViaSpotVol(dates, SwapExpiries, rates, Strikes, K, CapMaturity, spotvol)
% PRICESTRUCTUREDCAPS 
% Price a cap with different strikes for a given matrix of spot volatilities. Price is found calling PriceStructuredCaps with a single period
% INPUTS:
% dates:        Struct of dates of the financial instruments quoted in the market
% SwapExpiries: Vector of the quoted swap expiries
% rates:        Struct of dates of the financial instruments quoted in the market
% Strikes:      Vector of strikes of quoted volatilities
% K:            Strike of the cap
% CapMaturity:  Maturity of the cap
% spotvol:      Matrix of spot volatilities
%
% OUTPUT:
% Cap:            Price of the cap

% Yearly swap expiries up to 50 years:
FullSwapExpiries = [1:50]';
% Bootstrap:
[dates_bootstrap, discounts_bootstrap] = interpolateAndLaunchBootsrap(dates,rates,SwapExpiries,FullSwapExpiries);


% Find the payment dates of the Caps (every 3 months)
FLDates = findFloatingLegDates(datenum(dates.settlement),50, eurCalendar);
FLDates = [dates.settlement; FLDates];
% Interpolate the discounts at payment dates:
discounts3m = InterpDFviaRates(dates_bootstrap,discounts_bootstrap,FLDates);
% yearfrac Act/360:
DayCount = 2;
% yearfrac Act/365:
ACT365 = 3;
% Compute the year fractions:
deltas = yearfrac(FLDates(1:end-1),FLDates(2:end),DayCount);

% Fwd discounts vector
fwd_discounts = discounts3m(2:end)./discounts3m(1:end-1); %tutti i fwd
% Compute the forward Libor:
L = (1./fwd_discounts -1)./deltas;
% Number of caplets computaion
ncaplets = CapMaturity*4 -1;
% Price of the cap calling PriceStructuredCaps with a single period
Cap = PriceStructuredCaps(K, Strikes, L, deltas, FLDates, spotvol, ncaplets, discounts3m);


end