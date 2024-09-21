%% Assignment_6 Group_20, AY2023-2024
% Margherita Bencini, Alessandro Torazzi, Matteo Torba, Giovanni Urso
clc; clear; close all;
format long
% Read the market data about the volatilities
Caps_vol = xlsread("Caps_vol_20-2-24.xlsx","D2:R17");
Caps_vol = Caps_vol([1,3:size(Caps_vol,1)],:);
Caps_vol = Caps_vol(:,3:end);   % Volatilities matrix
Strikes = xlsread("Caps_vol_20-2-24.xlsx","F1:R1")'/100; % Vector of strikes
% Maturities in years:
YTM = [1:10,12,15]';
Caps_vol = Caps_vol(1:length(YTM),:);
% Notional value of the contract of the exercise:
Notional = 50000000;
%% a. Bootstrap the market discounts for the 16-feb-24
% Set the format data:
formatData='dd/mm/yyyy';
% Read market data for the bootstrap:
[dates, rates] = readExcelData("MktData_CurveBootstrap_20-2-24.xls", formatData);

% Quoted swap expiries:
SwapExpiries = [1:12 15 20 25 30 40 50]';
% Yearly swap expiries up to 50 years:
FullSwapExpiries = [1:50]';
% Bootstrap discounts up to 50 years:
[dates_bootstrap, discounts_bootstrap] = interpolateAndLaunchBootsrap(dates,rates,SwapExpiries,FullSwapExpiries);

disp('--- a. Bootstrap the market discounts for the 16-feb-24 --- ')
disp('Dates               Discounts')
for ii=1:length(dates_bootstrap)
    fprintf('%s         %.15f\n', datestr(dates_bootstrap(ii), 'dd-mmm-yyyy'),  discounts_bootstrap(ii));
end
disp(' ')
%% b. Determine the upfront X
% Start by determing the upfront via the flat vols:
X_flat = PriceXFlatVol(dates, SwapExpiries, rates, Strikes, Caps_vol, YTM);

% Calibrate the spot volatilities
spotvol = launchCalibration(dates,dates_bootstrap, discounts_bootstrap,Strikes,Caps_vol,YTM,1);
% Determine the upfront with via the spot volatilities
X = priceX(dates, SwapExpiries, rates, Strikes, spotvol);

disp('--- b. Determine the upfront X --- ')
fprintf('Via flat volatilities: X = %.6f%%\n', 100*X_flat)
fprintf('Via spot volatilities: X = %.6f%%\n\n', 100*X)
%% c. Delta bucket sensitivities
disp('--- c. Delta bucket sensitivities --- ')
disp(' ')
% Depos Delta bucket inizialization
XDeposBucket = zeros(4,1);
parfor ii = 1:4
    % Rates struct inizialization
    newrates = rates;
    % Increase each bucket iteratively
    newrates.depos(ii,:) = newrates.depos(ii,:) + 1e-4;
    fprintf('Depo rate %d of 4 shifting \n', ii)
    % Store the X value after rate shift
    XDeposBucket(ii) = priceX(dates, SwapExpiries, newrates, Strikes, spotvol);
end
% Futures Delta bucket inizialization
XFuturesBucket = zeros(7,1);
parfor ii = 1:7
    % Rates struct inizialization
    newrates = rates;
    % Increase each bucket iteratively
    newrates.futures(ii,:) = newrates.futures(ii,:) + 1e-4;
    fprintf('Future rate %d of 7 shifting \n', ii)
    % Store the X value after rate shift
    XFuturesBucket(ii) = priceX(dates, SwapExpiries, newrates, Strikes, spotvol);
end
% Swaps Delta bucket inizialization
XSwapsBucket = zeros(size(rates.swaps,1)-1, 1);
parfor ii = 2:size(rates.swaps,1)
    % Rates struct inizialization
    newrates = rates;
    % Increase each bucket iteratively
    newrates.swaps(ii,:) = newrates.swaps(ii,:) + 1e-4;
    fprintf('Swap rate %d of 17 shifting \n', ii-1)
    % Store the X value after rate shift
    XSwapsBucket(ii-1) = priceX(dates, SwapExpiries, newrates, Strikes, spotvol);

end
% Delta-Bucket sensitivities
% The delta is the difference of the upfront after and before the shift
DeposDelta = - X + XDeposBucket;
FuturesDelta = - X + XFuturesBucket;
SwapsDelta = - X + XSwapsBucket;

fprintf('Δ-Depos = %.10f%%\n', DeposDelta*100)
fprintf('Δ-Futures = %.10f%%\n', FuturesDelta*100)
fprintf('Δ-Swaps = %.10f%%\n',SwapsDelta*100)


% Compute the difference between sum of all delta-bucket sensitivities and
% total delta, doing also a plot
% Dates:
SwDates = findSwapDates(dates.settlement, SwapExpiries, eurCalendar)';
DeltaDates = [dates.depos(1:4); dates.futures(1:7,2); SwDates(2:end)];
% Sensitives:
DeltaSensitivities = [DeposDelta;FuturesDelta;SwapsDelta];
% Compute and plot the difference
diffDelta = DeltaDiffPlot(rates,dates, SwapExpiries, Strikes, X, DeltaSensitivities, DeltaDates, spotvol);
%% d. Total Vega
% Recalibrate the spot volatilities:
spotvol_shift = launchCalibration(dates,dates_bootstrap, discounts_bootstrap,Strikes,Caps_vol+1,YTM,0);
% Compute the upfront with 1 bp shifted flat vols
XTotalVega = priceX(dates, SwapExpiries, rates, Strikes, spotvol_shift);
% The vega is the difference of the upfront after and before the shift 
TotalVega = - X + XTotalVega;
disp(' ')
disp('--- d. Total Vega sensitivity --- ')
fprintf('Total Vega = %.2f € \n\n',Notional * TotalVega)

%% e. Vega-Bucket
disp('--- e. Vega bucket sensitivities --- ')
% Vega bucket inizialization
XVegaBucket = zeros(size(Caps_vol,1), 1);
%iterate on the years
parfor ii = 1:size(Caps_vol,1)
    % Flat vols inizialization
    NewCapsVol = Caps_vol;
    % Shift vols for each bucket iteratively
    NewCapsVol(ii, :) = NewCapsVol(ii, :) + 1;
    fprintf('Volatility bucket %d of %d shifting \n', ii, size(Caps_vol,1));
    % Recalibrate the spot volatilities:
    spotvol_new = launchCalibration(dates,dates_bootstrap, discounts_bootstrap,Strikes,NewCapsVol,YTM,0);
    % Compute the upfront with shifted vols by bucket
    XVegaBucket(ii) = priceX(dates, SwapExpiries, rates, Strikes, spotvol_new);
end
% Compute vega as the difference of upfront prices
VegaBucket = - X + XVegaBucket;
fprintf('Vega = %.2f €\n',VegaBucket*Notional)

% difference between sum of all vega-bucket sensitivities and total vega
diffVega = Notional*VegaDiffPlot(VegaBucket,TotalVega,YTM);
fprintf('Total Vega - sum of Vega buckets: %.10f%\n',diffVega);
disp(' ')
%% f. Coarse grained Delta sensitivities and hedging 
disp(' ')
disp('--- f. Coarse grained Delta sensitivities and hedging ---')
% Define the coarse grained peaks
CGBuckets = [2;5;10;15];
Deltas = [DeposDelta;FuturesDelta;SwapsDelta];
% Retrive the par swap rates for the swaps to be used in hedging
parSwapRates = mean(rates.swaps([2,5,10,13],:),2);
% Compute deltas for each CG bucket for upfront and swaps
SwDates = findSwapDates(dates.settlement, SwapExpiries, eurCalendar)';

BucketDates = [dates.depos(1:4);dates.futures(1:7,2);SwDates(2:end)];
DeltaBucketSw = DeltaBucketSwap(dates, SwapExpiries, rates, CGBuckets, parSwapRates);
[DeltaCGX, DeltaCGSwaps] = CGDelta(Deltas,dates,SwapExpiries,CGBuckets, BucketDates,DeltaBucketSw);
% Solve the linear system for the notionals matching each CG bucket delta
SwapsNotional = - DeltaCGSwaps\DeltaCGX * Notional;
for ii = 1 : length(CGBuckets)
    fprintf('Coarse grained Delta with peak in %d years = %f%% \n',CGBuckets(ii), DeltaCGX(ii)*100)
end
disp('--- Hedging with swaps ---')
for ii = 1 : length(CGBuckets)
    fprintf('Notional of swap with %d years maturity: %.2f € \n',CGBuckets(ii), SwapsNotional(ii))
end
figure()
hold on
plot([0;CGBuckets],[1;1;0;0;0],'LineWidth',2);
hold on
plot([0;CGBuckets],[0;0;1;0;0], 'LineWidth',2);
plot([0;CGBuckets],[0;0;0;1;0], 'LineWidth',2);
plot([0;CGBuckets],[0;0;0;0;1], 'LineWidth',2);
legend('2Y coarse grained shift','5Y coarse grained shift','10Y coarse grained shift','15Y coarse grained shift','Location', 'northwest',FontSize=18)
xlabel('Years', FontSize=18)
ylabel('Weights',FontSize=18)
title('Coarse grained weights',FontSize=22)


%% g. Vega hedging
disp(' ')
disp('--- g. Vega and Delta hedging ---')
% ATM strike for 5Y cap to be used in hedging
ATM5Y = mean(rates.swaps(5,:));

% Cap price computation
Cap5Y = priceCapViaSpotVol(dates, SwapExpiries, rates, Strikes, ATM5Y, 5, spotvol);
% Cap price with a full flat vol increase
Cap5Y_vega = priceCapViaSpotVol(dates, SwapExpiries, rates, Strikes, ATM5Y, 5, spotvol_shift);
% Vega of the cap as the difference of the prices
VegaCap = Cap5Y_vega - Cap5Y;
% Cap notional to hedge the total vega of the structured product
Cap5YNotional = -Notional*TotalVega/VegaCap ;

% Rate curve shift by 1 bp
allRatesShifted = rates;
allRatesShifted.depos = allRatesShifted.depos + 1e-4;
allRatesShifted.futures = allRatesShifted.futures + 1e-4;
allRatesShifted.swaps = allRatesShifted.swaps + 1e-4;

% Total delta as the difference of the upfront prices
X_totalDelta = priceX(dates, SwapExpiries, allRatesShifted, Strikes, spotvol) - X;
% Total delta of the 5Y cap
Cap5Y_totalDelta = priceCapViaSpotVol(dates, SwapExpiries, allRatesShifted, Strikes, ATM5Y,5, spotvol) - Cap5Y;
% Swap DV01 (since at par it is its BPV)
Swap_totalDelta = payerSwapNPV(dates, SwapExpiries, allRatesShifted, 5, mean(rates.swaps(5,:)));
% Swap notional computation , matching its delta with the portfolio's one
Notional5YSwap = - (Cap5YNotional*Cap5Y_totalDelta + Notional * X_totalDelta)/Swap_totalDelta;

fprintf('Vega for an ATM 5Y Cap = %f €\n',Cap5Y_vega);
fprintf('Notional of a ATM 5Y Cap to hedge the Vega for the structured product : %.2f €\n',Cap5YNotional);
fprintf('Delta for a 5Y Swap = %f €\n',Swap_totalDelta);
fprintf('Notional of a 5Y Swap to hedge the Delta for the portfolio: %.2f €\n',Notional5YSwap);


%% h. hedging of the coarse-grained vega
disp(' ')
disp('--- h. hedging of the coarse-grained Vega ---')
%Define the two peaks
VegaCGPeaks = [5;15];
%Define maturities vector
YTM = [1:10,12,15]';

%Define the Strike as the ATM15y Swap rate
ATM15Y = mean(rates.swaps(13,:));
%Define vectors of Strikes,Maturities and Coarse-Grained Buckets
CapStrike = [ATM5Y;ATM15Y];
CapMaturity = [5,15];
CGBuckets = [5;15];

%Compute the vega of the two Caps
VegaCaps = BucketVegaCap(dates, SwapExpiries, rates, Strikes, Caps_vol, YTM,CapStrike,CapMaturity,spotvol);
%Compute the Coarse-Grained Vega
[VegaCGX,VegaCGCaps] = CGVega(VegaBucket,dates,SwapExpiries,CGBuckets,YTM,VegaCaps);                       

%Hedging the Vega with a Cap, solving to find the notional
Notionals5_15 = - Notional * (VegaCGCaps\VegaCGX);
fprintf('Notional of 5Y Cap to Vega-hedge = %f €\n',Notionals5_15(1));
fprintf('Notional of 15Y Cap to Vega-hedge = %f €\n',Notionals5_15(2));


%Hedging also the Coarse Grained delta with a swap
%Compute mid-market swap rates
parSwapRates = mean(rates.swaps([5,13],:),2);

%Compute Delta Bucket of the Swaps
DeltaBucketSw = DeltaBucketSwap(dates, SwapExpiries, rates, CGBuckets, parSwapRates);
%Compute Delta Bucket of the Caps
DeltaCap = DeltaBucketCap(dates, SwapExpiries, rates, Strikes, CapMaturity, CapStrike, spotvol);
% Compute Coarse-Grained Delta of the Notional, the Swaps and the Caps
[DeltaCGX, DeltaCGSwaps, DeltaCGCaps] = CGDeltaCaps(Deltas, dates, SwapExpiries, CGBuckets, BucketDates, DeltaBucketSw, DeltaCap);

% Find the notionals of the Swaps that allow us to Delta-hedge
NotionalSwap5_15Y = - DeltaCGSwaps\(DeltaCGCaps* Notionals5_15 + Notional * (DeltaCGX));
fprintf('Notional of 5Y Swap to Δ-hedge = %f €\n',NotionalSwap5_15Y(1));
fprintf('Notional of 15Y Swap to Δ-hedge = %f €\n',NotionalSwap5_15Y(2));

