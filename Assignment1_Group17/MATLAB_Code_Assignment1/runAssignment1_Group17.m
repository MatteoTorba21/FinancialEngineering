%% Assignment_1
%  Group 17, AA2023-2024
%  Alessandro Torazzi, Matteo Torba, Giovanni Urso, Chiara Zucchelli
clc, clear, close

%% Pricing parameters
S0 = 1;                   % underlying price
K = 1;                    % strike price
r = 0.03;                 % ttm-zero-rate
TTM = 1/4;                % time-to-maturity
sigma = 0.22;             % volatility
flag = 1;                 % flag:  1 call, -1 put
d = 0.06;                 % dividend Yield
KI = 1.3;                 % barrier
threshold = 1e-4;         % threshold required for the numerical methods' error
S0min = 0.7; S0max = 1.5; % minimum and maximum underlying prices (for Vega calculation)
valueDate = '15/02/2008'; % value date
EEfreq = 1/12;            % early exercise frequency
dmin = 0; dmax = 0.06;    % maximum and minimum dividend Yields
Noptions = 1e6;           % numeber of contracts

%% Quantity of interest
B = exp(-r*TTM); % Discount

%% Pricing 
F0 = S0*exp(-d*TTM)/B;     % Forward in G&C Model
OptionPrice = zeros(1,3);
for pricingMode=1:3 % 1 ClosedFormula, 2 CRR, 3 MonteCarlo
    switch (pricingMode)
        case 1  % ClosedFormula
            M = 100;
        case 2  % CRR
            M = 100;
        case 3  % Monte-Carlo
            M = 10^6;
        otherwise
    end
    OptionPrice(pricingMode) = EuropeanOptionPrice(F0,K,B,TTM,sigma,pricingMode,M,flag);
end
% Price for 1 contract 
% printOptionPrices(OptionPrice,0)

NotionalOptionPrice = OptionPrice*Noptions; % option price for 1 Mln of contracts
printOptionPrices(NotionalOptionPrice,1)

%% Errors Rescaling 
% plot Errors for CRR varing number of steps
% Note: both functions plot also the Errors of interest as side-effect 
[nCRR,errCRR] = PlotErrorCRR(F0,K,B,TTM,sigma);

% plot Errors for MC varing number of simulations N 
[nMC,stdEstim] = PlotErrorMC(F0,K,B,TTM,sigma);

% We find the number of steps for CRR which reaches the threshold required:
M_opt_CRR = findM(errCRR,threshold,nCRR)

% We find the number of Monte-Carlo simulations which reaches the threshold required:
M_opt_MC = findM(stdEstim,threshold,nMC)

%% KI Option
% European barrier price via CRR
M = 1000;
euroBarrierCRR = EuropeanOptionKICRR(F0,K, KI,B,TTM,sigma,M);
Notional_euroBarrierCRR = euroBarrierCRR*Noptions % option price for 1 Mln of contracts

% European barrier price via MC
M = 1e6;
euroBarrierMC = EuropeanOptionKIMC(F0,K, KI,B,TTM,sigma,M);
Notional_euroBarrierMC = euroBarrierMC*Noptions % option price for 1 Mln of contracts

% Closed formula (digital + call with K = barrier)
% the amount of digital is (barrier - strike)
PlotEuroBarrierKIPayoff(K,KI,1);
euroBarrierClosed = EuropeanOptionKIClosed(F0,K, KI,B,TTM,sigma);
Notional_euroBarrierClosed = euroBarrierClosed*Noptions % option price for 1 Mln of contracts

%% KI Option Vega
% In order to plot the vega in function of the underlying, we pass a vector of forwards to the function VegaKI
S0vector = S0min:0.01:S0max;       % row vector for the underlying prices
F0vector = S0vector*exp(-d*TTM)/B; % row vector of the forwards 

% We compute the Vega of the option using two numerical methods (CRR and MC) and a closed formula
% We printed also the execution time of the three computations to determine which one is more efficient
figure();
hold on;
title('Vega')
for flagNum=1:3 % 1 CRR, 2 MonteCarlo, 3 closed formula 
    switch (flagNum)
        case 1 % CRR
            N = 5e3;
            tic; % Start timer
            vegaCRR = VegaKI(F0vector,K,KI,B,TTM,sigma,N,flagNum);
            elapsed_time = toc; % Stop timer
            % Print the elapsed time
            fprintf('Execution time to compute the Vega via the CRR: %.6f seconds\n', elapsed_time);
            % We multiply the vector of the computed vegas by the number of
            % otpions to have the vega expressed in euros                                                                                  
            plot(S0vector,vegaCRR*Noptions,'LineWidth',2); 
            
        case 2  % Monte-Carlo
            N = 1e5;
            tic;% Start timer
            vegaMC = VegaKI(F0vector,K,KI,B,TTM,sigma,N,flagNum);
            elapsed_time = toc; % Stop timer
            % Print the elapsed time
            fprintf('Execution time to compute the Vega via the MC: %.6f seconds\n', elapsed_time);
            % We multiply the vector of the computed vegas by the number of
            % options to have the vega expressed in euros
            plot(S0vector,vegaMC*Noptions,'LineWidth',2);

        case 3  % closed formula
            tic; % Start timer
            vegaClosedFormula = VegaKI(F0vector,K,KI,B,TTM,sigma,0,flagNum);
            elapsed_time = toc; % Stop timer
            % Print the elapsed time
            fprintf('Execution time to compute the Vega via the closed formula: %.6f seconds\n', elapsed_time);
            % We multiply the vector of the computed vegas by the number of
            % options to have the vega expressed in euros
            plot(S0vector, vegaClosedFormula*Noptions,'LineWidth',2);
    end
end
plot(S0vector, zeros(length(S0vector),1), 'k--', 'LineWidth',1)
legend('CRR', 'MC', 'Closed formula','FontSize', 24)
xlabel('S0', 'FontSize', 18)
ylabel('Vega (expressed in EUR)', 'FontSize', 18)
hold off

%% Antithetic Variables
% We compute through the antithetic variables method the prices of the
% options and we determine the error of this method
MCAVTerr = ErrorMCAVT(F0,K,B,TTM,sigma,flag);
% We verify graphically that the error is reduced with respect to the
% error of the 'classic' Monte-Carlo method
PlotMCAVTerror(nMC,stdEstim,MCAVTerr);

%% Bermudan Options
% We compute the price of the Bermudan option via the CRR approach
N = 1e3; % CRR steps
BermudanOptionPriceCRR = BermudanOptionCRR(F0,K,r,TTM,EEfreq,sigma,N);
NotionalBermudanOptionPrice = BermudanOptionPriceCRR*Noptions

%% Bermudan Option w.r.to Dividend Yields
% We plot the European Price and the Bermudan price for different dividend yields
n = 20; % discretization of the dividend yield
compareEuropeanBermudan(S0,K,r,TTM,sigma,dmin,dmax,n,EEfreq,Noptions)








