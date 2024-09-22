%% Assignment_2
% Group 17, AY2023-2024
% Alessandro Torazzi, Matteo Torba, Giovanni Urso, Chiara Zucchelli

clear all;
close all;
clc;

%% Settings
formatData='dd/mm/yyyy'; %Pay attention to your computer settings 

%% Read market data
% this fuction works on Windows OS. Pay attention on other OS.

[datesSet, ratesSet] = readExcelData('MktData_CurveBootstrap', formatData);

%% Bootstrap
% Computes Euribor 3m bootstrap with a single-curve model

% dates includes SettlementDate as first date
% determine the discounts via the bootstrap
[dates, discounts] = bootstrap(datesSet, ratesSet);

% display the results obtained in the bootstrap
disp('Dates obtained:')
disp(datestr(dates))
disp('Discount Factors obtained:')
disp(discounts)

%% Compute Zero Rates
% compute the Zero Rates from the Discount Factors
zRates = zeroRates(dates, discounts);

% display the Zero Rates obtained
disp('Zero Rates obtained:')
disp(zRates)

%% Plot Results
% Function to plot discount factors and zero rates
plotDFandRates(dates,discounts);

%% Exercise 2

shift_1bp = 1e-4;      % Shift of the curve
Notional = 1e7;        % Notional
fixedRate = 2.8173e-2; % Fixed rate

% Increase the rates used for the bootstrap of the original curve of 1 bp in parallel shift
ratesSetShifted = ratesSet;
ratesSetShifted.depos = ratesSet.depos + shift_1bp;
ratesSetShifted.futures = ratesSet.futures + shift_1bp;
ratesSetShifted.swaps = ratesSet.swaps + shift_1bp;

% Determine the Discount Factors for the shifted rates via the bootstrap
[dates_DV01, discounts_DV01] = bootstrap(datesSet, ratesSetShifted);

setDate = datesSet.settlement; % Settlement date
swapMaturity = 6; % Maturity of the swap which has to be priced
fixedLegPaymentDates = datesSet.swaps(1:swapMaturity); % Vector of the fixed leg payment dates

% Determine DV01, BVP and DV01_z
[DV01, BPV, DV01_z] = sensSwap(setDate, fixedLegPaymentDates, fixedRate, dates, discounts,discounts_DV01);
 
DV01 = DV01*Notional
BPV = BPV*Notional 
DV01_z = DV01_z*Notional
 
couponPaymentDates = fixedLegPaymentDates; % vector of coupon payment dates
% determine Macauly Duration
MacD = sensCouponBond(setDate, couponPaymentDates, fixedRate, dates, discounts)


%% Exercise 3

idx6y = find(dates == datesSet.swaps(6)); 
delta = yearfrac(dates(idx6y),dates(idx6y+1),6); 
BP = discounts(idx6y)+1-discounts(idx6y+1)*(1+delta*mean(ratesSet.swaps(7,:)));
disp('Bond Price:')
disp(BP*100) 


