%% Assignment_3 Group_17, AY2023-2024
% Alessandro Torazzi, Matteo Torba, Giovanni Urso, Chiara Zucchelli
clear all; close all; clc;

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

%% Exercise 1
disp('--- EXERCISE 1: ASSET SWAP SPREAD --- ')
disp(' ')

bondDirtyPrice = 1.01;                      % Price of the defaultable bond
couponRate = 3.9/100;                       % Annual coupon rate
Maturity = 3;                               % Maturity of the bond
paymentDates = datesSet.swaps(1:Maturity);  % Vector of the coupon payments

% Computation of the Asset Swap Spread:
S_asw = AssetSwapSpread(dates,discounts,bondDirtyPrice,couponRate,paymentDates,datesSet.settlement);
fprintf('The asset swap rate is equal to: s^{asw} = %.4f bps\n\n',S_asw*1e4)

%% Exercise 2: CDS Bootstrap
disp('--- EXERCISE 2: CDS BOOTSTRAP ---')
disp(' ')

ISP_R = 0.4;                                      % Recovery of ISP
datesCDS = datesSet.swaps(1:7);                   % CDS dates for known spreads
ISP_CDSSpreads = 1e-4*[29; 32; 35; 39; 40; 41];   % CDS Spreads given of ISP

% Spline interpolation of the missing CDS
MissingCDS = interp1([datesCDS(1:5); datesCDS(7)],ISP_CDSSpreads,datesCDS(6),'spline');
% Update of the CDS spreads vector:
idx = find(datesCDS>=datesCDS(6),1);
ISP_CDSSpreads = [ISP_CDSSpreads(1:idx-1); MissingCDS; ISP_CDSSpreads(idx:end)];
fprintf('ISP CDS spread foy t=6years obtained via spline interpolation: %.1f bps \n\n',MissingCDS*10000 )

% Intensities and survival probabilities approximation by neglecting the accrual:
flag = 1;
tic; % Start timer
[ISP_datesCDS_approx,ISP_survProbs_approx,ISP_intensities_approx] = bootstrapCDS(dates,discounts,datesCDS,ISP_CDSSpreads,flag,ISP_R);
elapsed_time = toc; % Stop timer
disp('1) Intensities approximation by neglecting the accrual:')
fprintf('Intensities (bp): λ(0,1)=%.1f; λ(1,2)=%.1f, λ(2,3)=%.1f, λ(3,4)=%.1f, λ(4,5)=%.1f, λ(5,6)=%.1f, λ(6,7)=%.1f \n',ISP_intensities_approx*10000 )
fprintf('Execution time of approximation: %.6f seconds\n \n', elapsed_time); % Print the elapsed time

% Exact intensities and survival probabilities (considering the accrual):
flag = 2;
tic; % Start timer
[IPS_datesCDS_exact,ISP_survProbs_exact,ISP_intensities_exact] = bootstrapCDS(dates,discounts,datesCDS,ISP_CDSSpreads,flag,ISP_R);
elapsed_time = toc; % Stop timer
disp('2) Exact intensities considering the accrual:')
fprintf('Intensities (bp): λ(0,1)=%.1f; λ(1,2)=%.1f, λ(2,3)=%.1f, λ(3,4)=%.1f, λ(4,5)=%.1f, λ(5,6)=%.1f, λ(6,7)=%.1f \n',ISP_intensities_exact*10000 )
fprintf('Execution time of exact computation: %.6f seconds\n \n', elapsed_time); % Print the elapsed time

% Intensities and survival probabilities via Jarrow-Turnbell model:
flag = 3;
tic; % Start timer
[ISP_datesCDS_JT,ISP_survProbs_JT,ISP_intensities_JT] = bootstrapCDS(dates,discounts,datesCDS,ISP_CDSSpreads,flag,ISP_R);
elapsed_time = toc; % Stop timer
disp('3) Intensities via Jarrow-Turnbull approximation:')
fprintf('Intensities (bp): %.1f (constant over the 7 years period by assumption) \n',ISP_intensities_JT(1)*10000 )
fprintf('Execution time of Jarrow-Turnbull approximation: %.6f seconds\n \n', elapsed_time); % Print the elapsed time

% Error by neglecting the accrual:
diff_accrual = max(abs(ISP_intensities_approx-ISP_intensities_exact));
disp('Error by neglecting the accrual')
fprintf('Error = max |λ_{approx}-λ_{exact}| =  %.10f bps \n \n', diff_accrual*1e4)

%% Plot of the intensities
% Interpretation of the Jarrow-Turnbell model:as can be seen graphically, 
% it is approximatly a mean of the previously computed exact/approx intensities.
% The displaied numerical computations confirm this supposition 
plotBootstrapCDSResults(datesSet.settlement,dates,discounts,datesCDS,ISP_CDSSpreads,ISP_R,ISP_intensities_approx,ISP_intensities_exact);

%% Exercise 3 
disp('--- EXERCISE 3: PRICE FIRST TO DEFAULT ---')
disp(' ')

UCG_R = 0.45;                                       % Recovery of UCG
datesCDS = datesSet.swaps(1:7);                     % CDS dates for known spreads
UCG_CDSSpreads = 1e-4*[34; 39; 45; 46; 47; 47];     % CDS Spreads of UCG:
MaturityDate = '20-feb-2012';                       % Maturity of the First to Default (FtD)
TTM = find(datesSet.swaps==datenum(MaturityDate));  % Time to maturity of the FtD
rho = 0.2;                                          % Correlation
N = 10000;                                          % Number of simulations

% Spline interpolation of the missing CDS of UCG:
MissingCDS = interp1([datesCDS(1:5); datesCDS(7)],UCG_CDSSpreads,datesCDS(6),'spline');
% Update of the CDS spreads of UCG vector:
idx = find(datesCDS>=datesCDS(6),1);
UCG_CDSSpreads = [UCG_CDSSpreads(1:idx-1); MissingCDS; UCG_CDSSpreads(idx:end)];
fprintf('UCG CDS spread foy t=6years obtained via spline interpolation: %.1f \n\n',MissingCDS*10000 )

% Exact intensities and survival probabilities of UCG:
flag = 1;
[UCG_datesCDS_approx,UCG_survProbs_approx,UCG_intensities_approx] = bootstrapCDS(dates,discounts,datesCDS,UCG_CDSSpreads,flag,UCG_R);
% Exact intensities and survival probabilities of UCG:
flag = 2;
[UCG_datesCDS_exact,UCG_survProbs_exact,UCG_intensities_exact] = bootstrapCDS(dates,discounts,datesCDS,UCG_CDSSpreads,flag,UCG_R);

% Spread of a First to Default, confidence interval and percentages of default
[s,CI,P] = firstToDefault2Names(dates,discounts,TTM,ISP_R,UCG_R,datesCDS,ISP_intensities_exact,UCG_intensities_exact,rho,N);
fprintf('Spread of the First To Default from MC: %.4f bps \nwith Confidence Interval with confidence level α=95%%: [%.4f %.4f] bps\n \n', s*1e4, CI(1)*1e4, CI(2)*1e4)
fprintf(['Percentage of no defaults: %.1f%% \n' ...
    'Percentage of only ISP default: %.1f%% \n' ...
    'Percentage of only UCG default: %.1f%% \n'...
    'Percentage of both ISP and UCG defaults: %.1f%% \nwith correlation: ρ=%.2f \n \n'], P(1,1), P(1,2), P(2,1), P(2,2), rho);

% Approx spread of a First to Default, confidence interval and percentages of default (by neglecting the accrual):
[sNoAccrual,CINoAccrual,PNoAccrual] = firstToDefault2NamesNoAccrual(dates,discounts,TTM,ISP_R,UCG_R,datesCDS,ISP_intensities_approx,UCG_intensities_approx,rho,N);
fprintf('Approximated Spread of the First To Default from MC BY NEGLECTING THE ACCRUAL: %.4f bps \nwith Confidence Interval with confidence level α=95%%: [%.4f %.4f] bps\n \n', sNoAccrual*1e4, CINoAccrual(1)*1e4, CINoAccrual(2)*1e4)

% Plot the default probabilities for different correlations
plotCorrelationRho(dates,discounts,TTM,ISP_R,UCG_R,datesCDS,ISP_intensities_exact,UCG_intensities_exact);

%% Check LI Model MC 
% Check the correct implementation of MonteCarlo Li model by comparing the 
% quoted and computed CDS spreads

[s_CDS_ISP, CI_CDS_ISP, P_CDS_ISP] = CDSspread(dates,discounts,TTM,ISP_R,datesCDS,ISP_intensities_exact,N);
fprintf('Quoted spread of the ISP CDS: %.4f bps\n',ISP_CDSSpreads(TTM)*1e4);
fprintf('Spread of the ISP CDS from MC: %.4f bps \nwith Confidence Interval with confidence level α=95%%: [%.4f %.4f] bps\n \n', s_CDS_ISP*1e4, CI_CDS_ISP(1)*1e4, CI_CDS_ISP(2)*1e4);

[s_CDS_UCG, CI_CDS_UCG, P_CDS_UCG] = CDSspread(dates,discounts,TTM,UCG_R,datesCDS,UCG_intensities_exact,N);
fprintf('Quoted spread of the UCG CDS: %.4f bps\n',UCG_CDSSpreads(TTM)*1e4);
fprintf('Spread of the UCG CDS from MC: %.4f bps \nwith Confidence Interval with confidence level α=95%%: [%.4f %.4f] bps\n \n', s_CDS_UCG*1e4, CI_CDS_UCG(1)*1e4, CI_CDS_UCG(2)*1e4);

%% Spread vs Correlation

N = 5000;                         % Number of Monte_Carlo simulations  
rho_vector = [-0.99:0.06:0.99]';  % Vector of correlations

% Plot the Spreads computed for different correlations:
plotSpreadForDifferentCorrelations(dates,discounts,TTM,ISP_R,UCG_R,datesCDS,ISP_intensities_exact,...
                                   UCG_intensities_exact,rho_vector,N,ISP_CDSSpreads,UCG_CDSSpreads)





