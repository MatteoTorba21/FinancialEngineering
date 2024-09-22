%% Assignment_5 Group_20, AY2023-2024
% Margherita Bencini, Alessandro Torazzi, Matteo Torba, Giovanni Urso
clear; close all; clc;

%% YearFrac Convenctions
DepoDayCount = 2;  % yearfrac Act/360
IBDayCount = 3;    % yearfrac Act/365
SwapDayCount = 6;  % yearfrac 30/360 European

%% Bootstrap
% Set the format data:
formatData='dd/mm/yyyy';
% Read market data:
[datesSet, ratesSet] = readExcelData('MktData_CurveBootstrap', formatData);
% Computes Euribor 3m bootstrap with a single-curve model dates: includes settlementDate as first date,
% determine the discounts via the bootstrap:
[dates, discounts] = bootstrap(datesSet, ratesSet);

%%  1. Case study: Certificate Pricing
% Settlement date:
settlementDate = dates(1);
% Time to maturity of the contract:
yearsTTM = 4;
% Payment dates of party A
FLDates = findFloatingLegDates(settlementDate,yearsTTM,[datenum('21-Mar-2008') datenum('24-Mar-2008') datenum('15-Ago-2008') datenum('24-Dec-2008') datenum('25-Dec-2008') datenum('26-Dec-2008') datenum('31-Dec-2008') datenum('01-Jan-2009') datenum('10-Apr-2009') datenum('13-Apr-2009') datenum('01-May-2009') datenum('24-Dec-2009') datenum('25-Dec-2009') datenum('31-Dec-2009') datenum('01-Jan-2010') datenum('02-Apr-2010') datenum('05-Apr-2010') datenum('24-Dec-2010') datenum('31-Dec-2010') datenum('22-Apr-2011') datenum('25-Apr-2011') datenum('15-Ago-2011') datenum('26-Dec-2011')]);
FLDates = [settlementDate; FLDates]; % Add the settlement date
% Discount factors at party A payment dates
discountsFinal = InterpDFviaRates(dates,discounts,FLDates);
% Payment dates of party B
yearlydates = FLDates(4*[0:yearsTTM]+1);
% Discount factors at party B payment dates
yearlyDiscounts = discountsFinal(4*[0:yearsTTM]+1);
% Time fractions
timeFractions = yearfrac(yearlydates(1:end-1),yearlydates(2:end),IBDayCount);
% Forward zero rates
r = -log(yearlyDiscounts(2:end)./yearlyDiscounts(1:end-1))./timeFractions;
% Given data:
sigma1 =16.1/100;       % Volatility of the ENEL stock
sigma2=20/100;          % Volatility of the AXA stock
div1 = 2.5/100;         % Dividend of the ENEL stock
div2 = 2.7/100;         % Dividend of the AXA stock
correlation = 40/100;   % Correlation between the two stocks
S1 = 100;               % ENEL stock price
S2 = 200;               % AXA stock price
S_pol = 100e-4;         % Spread over libor
X = 2/100;              % Upfront
P = 95/100;             % Protection

% Simulate the two stocks dynamics:
[asset1,asset2,asset1_minus,asset2_minus] = simulateAsset(S1,S2,sigma1, sigma2,r,timeFractions,div1,div2,correlation);
% Compute the ratios needed for the coupon value
SimAsset1 = sum(asset1(:,2:end)./asset1(:,1:end-1),2);
SimAsset2 = sum(asset2(:,2:end)./asset2(:,1:end-1),2);
SimAsset1_minus = sum(asset1_minus(:,2:end)./asset1_minus(:,1:end-1),2);
SimAsset2_minus = sum(asset2_minus(:,2:end)./asset2_minus(:,1:end-1),2);
% Compute the Payoff of the coupon of party A (which will be multiplied by alpha)
Payoff = max((SimAsset1 + SimAsset2)/8 - P,0);
Payoff_minus = max((SimAsset1_minus + SimAsset2_minus)/8 - P,0);
% Payoff_AV = (Payoff_minus+Payoff)/2;
Payoff_AV = [Payoff;Payoff_minus];
% Discounted payoff
discPayoff = (Payoff_AV)*yearlyDiscounts(end);
% NPV of the payments of party A:
LegA = 1 - yearlyDiscounts(end) + sum(discountsFinal(2:end).*yearfrac(FLDates(1:end-1), FLDates(2:end), DepoDayCount))*S_pol + (1-P)*discountsFinal(end);
% alpha value:
alpha = (LegA-X)./mean(discPayoff);

% Confidence interval: method 1 via normfit
confidence = 0.95; % confidence level

% Confidence intervals
t = norminv(confidence/2)
meanP = mean(discPayoff);
SE = sqrt(var(discPayoff)/length(discPayoff));
LB = (1/meanP - t * SE) * (LegA - X);
UB = (1/meanP + t * SE) * (LegA - X);


disp('--- 1. CASE STUDY: CERTIFICATE PRICING --- ')
disp(' ')
fprintf('Participation coefficient: α = %.6f  \nConfidence interval: [ %.6f , %.6f ] \n\n', alpha, LB, UB)

%% 2. Exercise: Pricing Digital option
load("cSelect20230131_B.mat")
% Notional of digital
Notional = 10e6;
% Coefficient of digital payoff
digitalPayoff = 5/100;
% Time to maturity set to 1 y
TTM = yearfrac(settlementDate,yearlydates(2),IBDayCount);
% Difference between digital Black and correct digital price
difference = Digital(cSelect,r(1),TTM,Notional,digitalPayoff);

disp('--- 2. EXERCISE: PRICING DIGITAL OPTION --- ')
disp(' ')
fprintf('Difference: diff = %.4f  \n\n', difference)
%% 3. Exercise: Pricing
runPricingFourier
%% 4. Case study: Volatility surface calibration
runCalibration

