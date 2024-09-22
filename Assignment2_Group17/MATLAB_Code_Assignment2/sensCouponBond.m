function MacD = sensCouponBond(setDate, couponPaymentDates, fixedRate, dates, discounts)
% computation of Macauly duration for an InterBank coupon bond
%
%INPUT
% setDate:              settlement date
% couponPaymentDates:   vector of dates in which we have coupon payments
% fixedRate:            fixed Interest Rate
% dates:                vector of dates (output from bootstrap)
% discount:             vector of discount factors (output from bootstrap)
%
%OUTPUT
% MacD:                 Macaulay duration

Notional = 1e7;         % notional
European_30_360 = 6;    % yearfrac 30/360 European

% vector with settlement date and coupon payment dates
couponDates = [setDate; couponPaymentDates];

% find the Discounts Factor at coupon payment dates by passing to the rates and linear interpolating
discountCouponPayments = InterpDFviaRates(dates, discounts, couponPaymentDates);

% time between two dates (in order to compute the coupons)
delta = [yearfrac(couponDates(1:end-1),couponDates(2:end),European_30_360)];

% vector of the coupons (for the last date we have to add the face value)
coupon = [fixedRate.*delta(1:end-1); 1+(fixedRate.*delta(end))];

% denominator of the Macauly duration 
MacD_den = sum(coupon.*discountCouponPayments)*Notional; 
 
% numerator of the Macauly duration
MacD_num = sum(coupon.*yearfrac(couponDates(1)*ones(length(couponDates)-1,1),...
    couponDates(2:end),European_30_360).*discountCouponPayments)*Notional;

% Macauly duration
MacD = MacD_num/MacD_den; 

% BondPrice = MacD_den; % price of the bond

% at par DV01_z = MacaulyDuration * BondPrice * 1bp
% proxy_DV01_z = MacD * BondPrice * 1e-4;   % (=5.201670995688768e+03)

end % function sensCouponBond

