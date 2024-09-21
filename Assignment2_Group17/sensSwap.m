function [DV01, BPV, DV01_z] = sensSwap(setDate, fixedLegPaymentDates, fixedRate, dates, discounts,discounts_DV01)
% computation of DV01, BPV and DV01_z
%
%INPUT
% setDate:                settlement date
% fixedLegPaymentDates:   vector of dates for fixed leg payments 
% fixedRate:              fixed Interest Rate
% dates:                  vector of dates (output from bootstrap)
% discount:               vector of discount factors (output from bootstrap)
% discount_DV01:          vector of discount factors (output from bootstrap with shifted rates)
%
%OUTPUT
% DV01:                   parallel shift of market rates
% BPV:                    Basic Point Value
% DV01_z:                 parallel shift of zero rates

European_30_360 = 6;  % yearfrac 30/360 European

% vector with settlement date and fixed leg payment dates
fixedDates = [setDate; fixedLegPaymentDates];

% find the Discounts Factor at fixed leg payment dates by passing to the rates and linear interpolating
discountFixedLegPayments = InterpDFviaRates(dates, discounts, fixedLegPaymentDates);

% Basis Point Value
BPV = discountFixedLegPayments'*(yearfrac(fixedDates(1:end-1), fixedDates(2:end),European_30_360));

fixedLegNPV = fixedRate*BPV;

floatingLegNPV = 1-discountFixedLegPayments(end);

NPV = -fixedLegNPV+floatingLegNPV;

% find the Discounts Factor at fixed leg payment dates by passing to the shifted rates and linear interpolating
discountFixedLegPayments_shifted = InterpDFviaRates(dates,discounts_DV01,fixedLegPaymentDates);

BPV_shifted = discountFixedLegPayments_shifted'*(yearfrac(fixedDates(1:end-1),fixedDates(2:end),European_30_360));

fixedLegNPV_shifted = fixedRate*BPV_shifted;

floatingLegNPV_shifted = 1-discountFixedLegPayments_shifted(end);

NPV_shifted = -fixedLegNPV_shifted+floatingLegNPV_shifted;

% parallel shift of market rates
DV01 = abs(NPV-NPV_shifted);


% vector of zero Rates
zZeroRates = zeroRates(dates,discounts); 

% vector of zero Rates shifted by 1bp
zZeroRates_shifted = zZeroRates + 1e-4;  

% linear interpolation of the shifted zero Rates
zInterpZeroRates_shifted = interp1(dates(2:end),zZeroRates_shifted,fixedLegPaymentDates);

% find the Discounts Factor at fixed leg payment dates by passing to the shifted zero rates and linear interpolating
zDiscount_shifted = DiscountsFromRates([dates(1); fixedLegPaymentDates],zInterpZeroRates_shifted);

zBPV_shifted = zDiscount_shifted'*(yearfrac(fixedDates(1:end-1),fixedDates(2:end),European_30_360));

zFixedLegNPV_shifted = fixedRate * zBPV_shifted;

zFloatingLegNPV_shifted = 1 - zDiscount_shifted(end);

zNPV_shifted = - zFixedLegNPV_shifted + zFloatingLegNPV_shifted;

% parallel shift of zero rates
DV01_z = abs(NPV - zNPV_shifted);


end % function sensSwap