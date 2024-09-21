function [Cap] = PriceStructuredCaps(Strikes, volStrikes, L, deltas, FLDates, spotvol, ncaplets, Discounts)
% PRICESTRUCTUREDCAPS 
% Price a cap with different strikes for each period by using the
% interpolated spot volatilities
% INPUTS:
% Strikes:      Vector of Strikes for the structured cap 
% volStrikes:   Vector of strikes of quoted volatilities
% L:            Vector of fwd Libor rates
% deltas:       Vector of deltas (ACT360 convention)
% FLDates:      Dates each 3 months includibg settlement 
% spotvol:      Interpolated spot vols
% ncaplets:     Number of same strike caplets for each period
% Discounts:    Discounts each 3 months
%
% OUTPUT:
% Cap:          Price of the structured cap
% Cap price inizialization
Cap = 0;
% Caplets number vectors used in computations
ncaplets = [0, ncaplets];
cumCaplets = cumsum(ncaplets);
% Different strikes cycle
for ii = 1:length(ncaplets)-1
    % volatilities inizialization
    sigma = zeros(ncaplets(ii), 1);
    % Volatilities at the required strikes obtained via spline interpolation on the strikes:
    sigma =  interp1(volStrikes, spotvol(cumCaplets(ii)+1:cumCaplets(ii+1), :)', Strikes(ii), "spline")';
    % d computation
    d = (L(cumCaplets(ii)+2:cumCaplets(ii+1)+1)-Strikes(ii))./(sigma*1e-4.*sqrt(yearfrac(FLDates(1),FLDates(cumCaplets(ii)+2:cumCaplets(ii+1)+1),3)));
    % Caplet prices vector for current period
    caplet = Discounts(cumCaplets(ii)+3:cumCaplets(ii+1)+2).*deltas(cumCaplets(ii)+2:cumCaplets(ii+1)+1).*...
        ((L(cumCaplets(ii)+2:cumCaplets(ii+1)+1)-Strikes(ii)).*normcdf(d)+sigma*1e-4.*sqrt(yearfrac(FLDates(1),FLDates(cumCaplets(ii)+2:cumCaplets(ii+1)+1),3)).*normpdf(d));
    % Cap price update
    Cap = Cap + sum(caplet);
end

end