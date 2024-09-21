function NPV_Flat = PriceXFlatVol(dates, SwapExpiries, rates, Strikes, Caps_vol, YTM)
% Function to determine the percentage value (wrt the notional) of the
% upfront in the contract given in the text via the flat volatilities: the
% upfront is determine by replicating the coupons with vanilla caps
%
% INPUTS:
% dates:        Struct of dates of the financial instruments quoted in the market
% SwapExpiries: Vector of the quoted swap expiries
% rates:        Struct of dates of the financial instruments quoted in the market
% Strikes:      Vector of strikes of quoted volatilities
% Caps_vol:     Matrix of volatilities
% YTM:          Maturities of the market volatilities
% PlotVolSmile: Flag variable
%               -If PlotVolSmile=1 plot the spot volatilities smile
%               -If PlotVolSmile=0 do not plot the volatilities smile
%
% OUTPUT:
% x:            Upfront's value of the contract

% Yearly swap expiries up to 50 years:
FullSwapExpiries = [1:50]';
% Bootstrap:
[dates_bootstrap, discounts_bootstrap] = interpolateAndLaunchBootsrap(dates,rates,SwapExpiries,FullSwapExpiries);

% Determine the upfront via the flat vols
% Start by computing the Caps price for the volatility grid
% Find the payment dates of the Caps (every 3 months)
FLDates = findFloatingLegDates(datenum(dates.settlement),50, eurCalendar);
FLDates = [dates.settlement; FLDates];  % add the settlement date
% Interpolate the discounts at payment dates:
discounts3m = InterpDFviaRates(dates_bootstrap,discounts_bootstrap,FLDates);
% yearfrac Act/360:
ACT360 = 2;
% yearfrac Act/365:
ACT365 = 3;
% Compute the year fractions:
deltas = yearfrac(FLDates(1:end-1),FLDates(2:end),ACT360);

% Strikes of the structured cap:
K = [0.032, 0.035, 0.04]; 
% Interpolation of the flat volatilities for the interesed strike
% Pricing of the used caps
% The behaviour of the coupons is replicated by a long and short position
% in casp with the respective strike of each period
sigma1 = interp1(Strikes, Caps_vol(5,:), K(1), "spline");
CAP(1) = capPrice(discounts3m,FLDates,deltas,K(1),sigma1,YTM(5));
sigma2 = interp1(Strikes, Caps_vol(5,:), K(2), "spline");
CAP(2) = -capPrice(discounts3m,FLDates,deltas,K(2),sigma2,YTM(5));
sigma3 = interp1(Strikes, Caps_vol(10,:), K(2), "spline");
CAP(3) = capPrice(discounts3m,FLDates,deltas,K(2),sigma3,YTM(10));
sigma4 = interp1(Strikes, Caps_vol(10,:), K(3), "spline");
CAP(4) = -capPrice(discounts3m,FLDates,deltas,K(3),sigma4,YTM(10));
sigma5 = interp1(Strikes, Caps_vol(12,:), K(3), "spline");
CAP(5) = capPrice(discounts3m,FLDates,deltas,K(3),sigma5,YTM(12));
% BPV of a swap neglecting the first coupon
BPV = discounts3m(3:61)'*deltas(2:60);
% NPV of the structured product
NPV_Flat = 0.009*BPV + 1 - discounts3m(2) - 0.01*discounts3m(2)'*deltas(1) + sum(CAP);
end

