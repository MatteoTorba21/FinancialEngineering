function S_asw = AssetSwapSpread(dates,discounts,BondDirtyPrice,couponRate,paymentDates,settlementDate)
% Function to compute the asset swap spread Over Euribor3m
%
% INPUTS:
% dates:                Vector of dates obtained from the bootstrap
% discounts:            Vector of dicsounts obtained from the bootstrap
% BondDirtyPrice:       Price of the defaultable bond
% couponRate:           Annual coupon rate of the bond
% paymentDates:         Vector of the payment dates of the bond
%
% settlementDate:       Settlement date
% OUTPUTS:
% S_asw:                Asset swap spread
%

SwapDayCount = 6; % Yearfrac fixed leg convenction: 30/360 European
FLCount = 2;      % Yearfrac floating leg convenction: Act/360

% Fixed leg
% Discount factors at the payment dates:
DFAtPayment = InterpDFviaRates(dates,discounts,paymentDates);
% Deltas for the fixed leg:
deltas = yearfrac([settlementDate; paymentDates(1:end-1)],paymentDates,SwapDayCount);
% Price of the non defaultable Bond:
BondPrice = couponRate*deltas'*DFAtPayment+1*DFAtPayment(end);

% Floating leg
% Floating leg dates:
FloatingLegDates = findFloatingLegDates(settlementDate,3);
% Deltas for the floating leg:
deltas_floatingleg = yearfrac([settlementDate; FloatingLegDates(1:end-1)],FloatingLegDates,FLCount);
% Discount factors fot the floating leg cashflows:
DF_floatingleg = InterpDFviaRates(dates,discounts,FloatingLegDates);
% Computation of the BPV for the floating leg:
BPV_floatingleg = deltas_floatingleg'*DF_floatingleg;

% Asset swap spread:
S_asw = (BondPrice-BondDirtyPrice)/BPV_floatingleg;

end % function AssetSwapSpread