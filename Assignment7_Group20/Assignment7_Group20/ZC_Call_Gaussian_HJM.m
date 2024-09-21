function Price = ZC_Call_Gaussian_HJM(sigma,discounts,fwd_discounts,Expiry,Maturity,K,t)
% Price a call option written on a zero coupon bond in the Gaussian HJM framework.
%
% INPUTS:
% sigma:            Function handle to compute sigma
% discounts:        Discount factors form t to Maturity
% fwd_discounts:    Forward discount factors computed in t from Maturity to
%                   Expiry
% Expiry:           Expiry of the call option
% Maturity:         Maturity of the ZCB
% K:                Strike
% t:                Instant in which the call price is computed
%
% OUTPUT:
% Price:            Price of the ZC call option

% Time to maturity
TTM = Expiry-t;
% Integrand to compute V^2
integrand = @(u) (sigma(u,Maturity)-sigma(u,Expiry)).^2;
% V^2 value:
V_squared = 1./(TTM) .*integral(integrand, t, Expiry,'ArrayValued',true,"AbsTol",1e-12,"RelTol",1e-12);
% d1 and d2 values:
d1 = log(fwd_discounts./K)./(sqrt(V_squared*TTM)) + 1/2*sqrt(V_squared*TTM);
d2 = log(fwd_discounts./K)./(sqrt(V_squared*TTM)) - 1/2*sqrt(V_squared*TTM);

% Price of the call option:
Price = discounts.*(fwd_discounts.*normcdf(d1)-K.*normcdf(d2));

end