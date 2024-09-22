function h_curve = check_h(cf_schedule_1y,PV_1y,cf_schedule_2y,PV_2y,ZC_curve,R)
% Compute the hazard rate by reversing the formula to compute the price
%
% INPUT:
% cf_schedule_1y:      Table of cash flows of corp. bonds with expiry 1y with
%                      -column #1: cash flow date (year frac)
%                      -column #2: cash flow amount (US $)
% PV_1y:               Price of the risky bond (dirty) with expiry 1y
% cf_schedule_2y:      Table of cash flows of corp. bonds with expiry 2y with
%                      -column #1: cash flow date (year frac)
%                      -column #2: cash flow amount (US $)
% PV_2y:               Price of the risky bond (dirty) with expiry 2y
% ZC_curve:            Table of zero-coupon rates (continuous compounding)
%                      -column #1: maturity (year frac)
%                      -column #2: MID rate
% R:                   Recovery rate
%
% OUTPUT:
% h_curve:             Table of piece-wise hazard rates

% Spline interpolation of the zero coupon rates for all the cash flow dates:
rates = interp1(ZC_curve(:,1),ZC_curve(:,2),cf_schedule_2y(:,1),'spline');
% Discounts for pure interest rates:
discounts = exp(-cf_schedule_2y(:,1).*rates);

% Find the hazard rate between t=0 and t=1y
% Solve the reversed equation in function of k=exp(-h1*0.5)
a = cf_schedule_1y(2,2)*discounts(2)-R*100*discounts(2);
b = cf_schedule_1y(1,2)*discounts(1)+R*100*(discounts(2)-discounts(1));
c = 100*R*discounts(1)-PV_1y;

k = (-b+sqrt(b^2-4*a*c))/(2*a);
% Find h1:
h1 = -2*log(k);
% Update the result:
h_curve = [1.00 h1];

% Find the hazard rate between t=1 and t=2y
% Solve the reversed equation in function of k=exp(-h2*0.5)
a = cf_schedule_2y(4,2)*discounts(4)*exp(-h1)-100*R*exp(-h1)*discounts(4);
b = cf_schedule_2y(3,2)*discounts(3)*exp(-h1)+100*R*exp(-h1)*(discounts(4)-discounts(3));
c = cf_schedule_2y(1,2)*discounts(1)*exp(-h1*0.5)+cf_schedule_2y(2,2)*discounts(2)*exp(-h1)+...
    100*R*discounts(1)*(1-exp(-h1*0.5))+100*R*discounts(2)*(exp(-h1*0.5)-exp(-h1))+100*R*discounts(3)*exp(-h1)-PV_2y;

k = (-b+sqrt(b^2-4*a*c))/(2*a);
% Find h2:
h2 = -2*log(k);
% Update the result:
h_curve = [h_curve; 2.00 h2];

end % function check_h
