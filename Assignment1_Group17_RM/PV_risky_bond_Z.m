function [ PV ] = PV_risky_bond_Z( z, cf_schedule, ZC_curve)
% Dirty price for a given risky bond (from scalar Z-spread)
%
% INPUT:
% z:                Z-spread of the bond
% cf_schedule:      Table of cash flows of corp. bonds with
%                   -column #1: cash flow date (year frac)
%                   -column #2: cash flow amount (US $)
% ZC_curve:         Table of zero-coupon rates (continuous compounding)
%                   -column #1: maturity (year frac)
%                   -column #2: MID rate
%
% OUTPUT:
% PV:               Price of the risky bond (dirty)

% Define a vector of cash flows dates for simplicity:
% payment_dates = cf_schedule(:,1);

% Spline interpolation of the zero coupon rates for all the cash flow dates:
rates = interp1(ZC_curve(:,1),ZC_curve(:,2),cf_schedule(:,1),'spline');
% Discounts for pure interest rates:
discounts = exp(-cf_schedule(:,1).*rates); 

% Computation of the survival probabilities:
P = exp(-z*cf_schedule(:,1));

% Computation of the (dirty) price of the risky bond:
PV = sum(discounts.*cf_schedule(:,2).*P);

end % function PV_risky_bond_Z