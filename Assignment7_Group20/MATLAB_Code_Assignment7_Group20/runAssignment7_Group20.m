%% Assignment_7 Group_20, AY2023-2024
% Margherita Bencini, Alessandro Torazzi, Matteo Torba, Giovanni Urso
clc; clear; close all;

%% Case study: Certificate Pricing
% Load dates and discounts from the bootstrap done in Assignment 2:
load("dates.mat");
load("discounts.mat");
% Load data from the struct given for Assignment 5:
load("cSelect20230131_B.mat")
% Convenctions for yearfrac function:
ACT360 = 2; ACT365 = 3; ACT30360 = 6;
% Parameter in the Laplace exponent (for NIG)
alpha = 0.5;
% Calibrate the NIG parameters:
[sigma, volvol, nu] = NMMCalibration(alpha);

disp('--- 1. CASE STUDY: CERTIFICATE PRICING --- ')
fprintf(['Calibrated parameters for NIG model: \n' ...
    'σ (average volatility):  %f    ' ...
    'η (volatility skew): %f    ' ...
    'k (vol of vol): %f\n\n'], sigma, nu, volvol);


%% a. Value the upfront X%
% Strike
K = 3200;
% Dividend
div = cSelect.dividend;
% Notional 
Principal = 1e8;
% Load the spot
S0 = cSelect.reference;
% Fix seed
rng(20)
% Find the floating leg Dates up to 3 years and the correspective discounts
% (useful for point d)
Maturity = 3;
FLDates_3y = findFloatingLegDates(dates(1),Maturity, eurCalendar);
FLDates_3y = [dates(1); FLDates_3y];
discounts3m_3years = InterpDFviaRates(dates,discounts,FLDates_3y);
% Save the values up to 2 years
FLDates_2y = FLDates_3y(1:9);
discounts3m_2years = discounts3m_3years(1:9);
deltasYearly = yearfrac(dates(1),FLDates_2y(5),ACT365);
% Compute the forward F(0,1)
F0 = S0*exp(-div*deltasYearly)/discounts3m_2years(5);
% Find the coupon reset date
CouponResetDate =  addtodate(dates(1), 364, 'day');
if ~isbusday(CouponResetDate,eurCalendar)
% If the date obtained is not a business day, the next business day is taken,
% according to the 'modified follow' convenction
    CouponResetDate = busdate(CouponResetDate, 'modifiedfollow', eurCalendar);
end

% Compute the zero rates
zRates = zeroRates(dates,discounts);
% Compute time to maturity of the reset date
TTM = yearfrac(dates(1),CouponResetDate,ACT365); 
% Compute the zero rate for 1Y
zRate1year = interp1(dates(2:end),zRates,FLDates_2y(5));
% Compute the zero rate and discount for reset date
zRateResetDate1y = interp1(dates(2:end),zRates,CouponResetDate);
discountReset = InterpDFviaRates(dates,discounts,CouponResetDate);
% Compute the forward at reset date
F0ResetDate = F0*exp(zRateResetDate1y*TTM - zRate1year*deltasYearly-div*(TTM-deltasYearly));
% Compute the log moneyness 
k = log(F0ResetDate/K);

% COMPUTATION OF THE UPFRONT VIA THE CLOSED FORMULA (USING LEWIS FORMULA)
% Compute the probability as the price of a 1-year digital option
% Logarithm of Laplace transform function
lnL = @(w) TTM./volvol*(1-alpha)/alpha*(1-(1+(w.*volvol*sigma^2)/(1-alpha)).^alpha);
% Characteristic function
Phi = @(x) exp(-1i*x.*lnL(nu) + lnL((x.^2 + 1i*(1+2*nu).*x)./2));
% Function to be integrated
f = @(u) real(exp(1i*u*k).*Phi(u)./(1i*u));
% Integration
I = integral(@(u) f(u), 0, Inf, "AbsTol",1e-16,"RelTol",1e-16); 
% Compute the probability
ProbNOTER = 1/2 + 1/pi*I;

% Compute the year fractions
deltas = yearfrac(FLDates_2y(1:end-1),FLDates_2y(2:end),ACT360);

% Compute the NPV of party A in case of Early Redemption
spread = 1.3/100;
legA_1y = 1-discounts3m_2years(5)+(discounts3m_2years(2:5)'*deltas(1:4))*spread;
% Compute the NPV of party A when it arrives to maturity
legA_2y = 1-discounts3m_2years(end)+(discounts3m_2years(2:end)'*deltas(1:end))*spread;
% Compute the Coupon paid by B in case of Early Redemption
Coupon1y = 6/100*yearfrac(dates(1),FLDates_2y(5),ACT30360)*discounts3m_2years(5);
% Compute  the Coupon paid by B when it arrives to maturity
Coupon2y = 2/100*yearfrac(FLDates_2y(5),FLDates_2y(end),ACT30360)*discounts3m_2years(end);
% Compute the NPV without the Upfront
NPV = (legA_1y- Coupon1y)*(1-ProbNOTER) + ProbNOTER*(legA_2y-Coupon2y);
% Compute the Upfront

X = NPV*Principal;
disp('1a. Value the upfront X%')
fprintf('Upfront via Lewis closed formula : X = %.2f € \n', X)

%% MonteCarlo check
% Compute F0 at Reset Date
F0ResetDate = F0*exp(zRateResetDate1y*TTM - zRate1year*deltasYearly-div*(TTM-deltasYearly));
% Simulate Asset
rng(6)
% Number of simulations 
Nsim = 2e5;
% Asset simulation using NIG dynamics
Asset = SimulateNIG(deltasYearly,sigma,volvol,nu,Nsim/2,F0ResetDate,[1; discountReset],div,lnL);
% NPV of the contract calculation for wìeach simulation
Payoff_simulation = ((Asset<K)*(legA_1y- Coupon1y) + (Asset>=K)*(legA_2y-Coupon2y))*Principal;
% Mean and confidence interval for X%
confidence = 0.05;
[NPVMC,~,CIMC] = normfit(Payoff_simulation,confidence);

fprintf('Upfront via Monte-Carlo: X = %.2f € \nConfidence interval at %.0f%% : [ %.2f € , %.2f € ] \n\n', NPVMC, 100*(1 - confidence), CIMC(1), CIMC(2))
%% b.+e. Closed Black formula for the upfront value
% Obtain Black impied volatility
blksigma = interp1(cSelect.strikes,cSelect.surface,K);
% Struct initialization 
data = cSelect;
data.K = K;
% Black N(d_2) computation
d2 = log(F0ResetDate/K)/(sqrt(TTM)*blksigma) - 1/2*sqrt(TTM)*blksigma;
ProbERBlack = normcdf(d2);
% Upfront using Black model
NPVBlack =(legA_1y- Coupon1y)*(1-ProbERBlack) + ProbERBlack*(legA_2y-Coupon2y);
% Error between NIG and Black models
errBlack = NPV-NPVBlack;
disp('1b. & 1e. Closed Black formula for the upfront value')
fprintf('Upfront using Black model: X = %.2f € \nDifference between the NIG and Black models: %.2f € \n', NPVBlack*Principal, errBlack*Principal)

% Black digital risk correction
% Difference in the digital price considering the digital risk and corrected N(d_2)
difference = Digital(data,zRateResetDate1y,TTM,1,1)*exp(zRateResetDate1y*TTM);
ProbDigitalCorrected = difference + ProbERBlack;
% NPV with corrected probabilities
NPVBlackCorrected =(legA_1y- Coupon1y)*(1-ProbDigitalCorrected) + ProbDigitalCorrected*(legA_2y-Coupon2y);
fprintf('Upfront using Black model considering Digital risk: X = %.2f € \n\n', NPVBlackCorrected*Principal)

%% d. Structured bond with a three-year expiry
% Reset date for second year computation
ResetDate2y = addtodate(dates(1), 366+363, 'day');
if ~isbusday(ResetDate2y,eurCalendar)
% If the date obtained is not a business day, the next business day is taken,
% according to the 'modified follow' convenction
    ResetDate2y = busdate(ResetDate2y, 'modifiedfollow', eurCalendar);
end
% Interpolate the zero rate
zRateResetDate2y = interp1(dates(2:end),zRates,ResetDate2y);
% yearfrac for second reset date
TTM2y = yearfrac(dates(1),ResetDate2y,ACT365);
% Forward at 2Y Reset date (for MC NIG simulation)
F0ResetDate2y = S0*exp(zRateResetDate2y*TTM2y-div*(TTM2y));
discountReset2y = InterpDFviaRates(dates,discounts,ResetDate2y);
% Simulate Asset
Nsim = 1e5;
Asset = SimulateNIG([TTM;TTM2y-TTM],sigma,volvol,nu,Nsim/2,F0ResetDate2y,[1; discountReset; discountReset2y],div,lnL);

% Compute the year fractions up to 3 years
deltas_3y = yearfrac(FLDates_3y(1:end-1),FLDates_3y(2:end),ACT360);
% Compute the NPV of party A when it arrives to maturity
legA_3y = 1-discounts3m_3years(end)+(discounts3m_3years(2:end)'*deltas_3y(1:end))*spread;
% Compute the Coupon paid by B in case of Early Redemption a 1 year
Coupon1y = 6/100*yearfrac(dates(1),FLDates_3y(5),ACT30360)*discounts3m_3years(5);
% Compute the Coupon paid by B in case of Early Redemption at 2 years
Coupon2y = 6/100*yearfrac(FLDates_3y(5),FLDates_3y(9),ACT30360)*discounts3m_3years(9);
% Compute  the Coupon paid by B when it arrives to maturity (3 years)
Coupon3y = 2/100*yearfrac(FLDates_3y(9),FLDates_3y(end),ACT30360)*discounts3m_3years(end);

% Compute the NPV for each simulation
Payoff_simulation_3Y = (- (Coupon1y)*(Asset(:,1)<K)  - (Coupon2y)*(Asset(:,2)<K)+(legA_3y-Coupon3y)*(Asset(:,2)>K & Asset(:,1)>K) + ...
+ (Asset(:,2)<K | Asset(:,1)<K) * legA_2y) * Principal;
% Mean and confidence interval for the upfront in 3Y case
[NPVMC_3Y,~,CIMC_3Y] = normfit(Payoff_simulation_3Y,confidence);
disp('1d. Structured bond with a three-year expiry')
fprintf('Upfront for 3Y contract via Monte-Carlo: X = %.2f € \nConfidence interval at %.0f%% : [ %.2f € , %.2f € ] \n\n', NPVMC_3Y, 100*(1 - confidence), CIMC_3Y(1), CIMC_3Y(2))

%% Bermudian Swaption Pricing via Hull-White
% Run the script to price the Bermudan
TrinomialTree;




