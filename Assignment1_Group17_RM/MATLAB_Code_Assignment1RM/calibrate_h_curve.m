function h_curve = calibrate_h_curve(ZC_curve,cf_schedule_1y, Bond_dirty_price_1y,cf_schedule_2y,Bond_dirty_price_2y,R)
% Calibration of the hazard curve through a numerical method given two risky
% bonds with maturities 1 year and 2 years respectively, the zero coupon
% curve and the recovery rate
%
% INPUT:
% ZC_curve:             Table of zero-coupon rates (continuous compounding)
%                       -column #1: maturity (year frac)
%                       -column #2: MID rate
% cf_schedule_1y:       Table of cash flows of corp. bonds with maturity 1 year
%                       -column #1: cash flow date(year frac)
%                       -column #2: cash flow amount (US $)
% Bond_dirty_price_1y:  Price of the risky bond (dirty) with maturity 1 year                  
% cf_schedule_2y:       Table of cash flows of corp. bonds with maturity 2 years
%                       -column #1: cash flow date (year frac)
%                       -column #2: cash flow amount (US $)
% Bond_dirty_price_2y:  Price of the risky bond (dirty) with maturity 2 years
% R:                    Recovery rate
%
% OUTPUT:
% h_curve:              Table of piece-wise hazard rates

% First year hazard rate calibration:
eqn = @(h) PV_risky_bond_h(cf_schedule_1y,[1 h],ZC_curve,R)-Bond_dirty_price_1y;
h1 = fzero(eqn,0);
% Second year hazard rate calibration:
eqn = @(h) PV_risky_bond_h(cf_schedule_2y,[1 h1; 2 h],ZC_curve,R)-Bond_dirty_price_2y;
h2 = fzero(eqn,0);

h_curve = [1.00 h1; 2.00 h2];

end % function calibrate_h_curve