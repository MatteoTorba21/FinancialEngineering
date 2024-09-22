function z = calibrate_z(ZC_curve,cf_schedule, Bond_dirty_price)
% Computation of the z-spread through a numerical method given a risky bond
%
% INPUT:
% ZC_curve:             Table of zero-coupon rates (continuous compounding)
%                       -column #1: maturity (year frac)
%                       -column #2: MID rate
% cf_schedule:          Table of cash flows of corp. bonds with
%                       -column #1: cash flow date(year frac)
%                       -column #2: cash flow amount (US $)
% Bond_dirty_price:     Price of the risky bond (dirty)                  
%
% OUTPUT:
% z:                    Z-spread of the bond

% Equation to be solved
eqn = @(z_var) PV_risky_bond_Z(z_var,cf_schedule,ZC_curve)-Bond_dirty_price;
% Numerical solution
z = fzero(eqn,0);

end % function calibrate_z