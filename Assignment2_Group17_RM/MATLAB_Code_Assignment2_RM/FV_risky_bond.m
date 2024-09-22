function [ FV ] = FV_risky_bond(cf_schedule, Q, ZC_curve, R)
% Function to compute the dirty 1y forward prices for a given risky bond
% INPUTS:
% cf_schedule:  Table of cash flows of corporate bonds such that
%               -column #1: cash flow date(year frac)
%               -column #2: cash flow amount (US $)
% Q:            Rating transition matrix
% ZC_curve:     table of zero-coupon rates (continuous compounding) such that
%               -column #1: maturity (year frac)
%               -column #2: MID rate
% R:            Recovery rate
% OUTPUTS:
% FV:           Column vector of dirty 1y forward prices, each row
%               corresponding to a rating 
%               (row #1: IG; row #2: HY; row #3: Defaulted)

% Preallocation of the memory for the vector FV:
FV = zeros(3,1);

% Schedule of remaining payments:
% cf_schedule = cf_schedule(3:4,:);

% Zero coupon curve at the cash flows dates obtained via linear interpolation:
ZeroRates = interp1(ZC_curve(:,1),ZC_curve(:,2),cf_schedule(2:end,1),"linear");
% Discounts at the cash flows dates computed in t=1y
discounts = exp(-(cf_schedule(2:end,1)).*ZeroRates);
% Forward discounts at the cash flow dates:
fwd_discounts = discounts(2:end)./discounts(1);

% Alternative method to compute the default probabilities at 6 months
% by computing the transistion matrix at 6 months to have by using the
% time-homogeneity of the Markov chain process:
% Q_6m = sqrtm(Q);
 
% Hazard rate (assumed constant):
h_IG = -log(1-Q(1,3));
h_HY = -log(1-Q(2,3));

% Survival probabilities at 1y, 1.5y and 2y, computed at 1y in the IG case:
% survProb_IG = [1; 1-Q_6m(1,3); 1-Q(1,3)]; % by the alternative method
survProb_IG = [1; exp(-h_IG.*[0.5; 1])];
% Survival probabilities at 1y, 1.5y and 2y, computed at 1yin the HY case::
% survProb_HY = [1; 1-Q_6m(2,3); 1-Q(2,3)]; % by the alternative method
survProb_HY = [1; exp(-h_HY.*[0.5; 1])];

% Calculate forward values
% Case of remaining IG:
FV(1) = cf_schedule(3:end,2)'*(fwd_discounts.*survProb_IG(2:3))+100*R*(survProb_IG(1:2)-survProb_IG(2:3))'*fwd_discounts; % ?????
 
% Case of migrating to HY:
FV(2) = cf_schedule(3:end,2)'*(fwd_discounts.*survProb_HY(2:3))+100*R*(survProb_HY(1:2)-survProb_HY(2:3))'*fwd_discounts; % ?????

% Case of default:
FV(3) = R*100;


end % function FV_risky_bond