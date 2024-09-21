function [ PV ] = PV_risky_bond_h(cf_schedule, h_curve, ZC_curve, R)
% Dirty price for a given risky bond (from hazard rate curve- fixed recoovery rate R)
%
% INPUT:
% cf_schedule:      Table of cash flows of corp. bonds with
%                   -column #1: cash flow date (year frac)
%                   -column #2: cash flow amount (US $)
% h_curve:          Table of piece-wise hazard rates
% ZC_curve:         Table of zero-coupon rates (continuous compounding)
%                   -column #1: maturity (year frac)
%                   -column #2: MID rate
% R:                Recovery rate
%
% OUTPUT:
% PV:               Price of the risky bond (dirty)

% Spline interpolation of the zero coupon rates for all the cash flow dates:
rates = interp1(ZC_curve(:,1),ZC_curve(:,2),cf_schedule(:,1),'spline');
% Discoun cqxts for pure interest rates:
discounts = exp(-cf_schedule(:,1).*rates);

% Create a vector in which every element of the second column of h_curve is replicated twice
h = repelem(h_curve(:,2),2);
% If h is not a column, trasform it in a column vector
h = reshape(h,[],1);

% Computation of the vector of the survival probabilities:
P = exp(-cumsum(h*0.5));

% Vector of rate of default between two coupon payments T1 and T2 as seen in t with t < T1 < T2:
rate_default_in_coupondates = [1 - P(1); P(1:length(cf_schedule(:,1))-1)-P(2:length(cf_schedule(:,1)))];

% Computation of the (dirty) price of the risky bond:
PV = sum(discounts.*cf_schedule(:,2).*P(1:length(cf_schedule(:,1))))+R*100*sum(rate_default_in_coupondates.*discounts);

end % function PV_risky_bond_h
