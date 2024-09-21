%% Trinomial tree to price a swaption
% Load dates and discounts from the bootstrap done in Assignment 2:
load("dates.mat");
load("discounts.mat");
% Convenctions for yearfrac function:
ACT360 = 2; ACT365 = 3; ACT30360 = 6;

T_start = 2;    % Start time to enter in the swaption
K = .05;        % Swaption's strike
Maturity = 10;  % Maturity of the underlying swap
N = 1000;       % Number of time steps in a year

% HW parameters
a = 11/100;
sigma = 0.8/100;

% Find the yearly dates up to 10 years (using the modified follow convenction)
% and the correspective discounts
FLDates_10y = findLegDatesYearly(dates(1),Maturity, eurCalendar);
FLDates_10y = [dates(1); FLDates_10y];
discountsyearly = InterpDFviaRates(dates,discounts,FLDates_10y);
% Keep only the discounts from 2 years on
discountsyearly = discountsyearly(T_start+1:end) ;
% Compute the yearly year fractions starting from 2 years:
deltas = yearfrac(FLDates_10y(1+T_start:end-1),FLDates_10y(2 + T_start:end),ACT30360);

% From the bootsrap results compute the forward discounts matrix: B(t0=0,ti,tj) for ti=2,...,9 and tj=ti,...,10 
fwdD_matrix = zeros(length(discountsyearly),Maturity-T_start);
for ii=1:Maturity-T_start
    fwdD_matrix(ii:(Maturity-T_start +1),ii) = discountsyearly(ii:end)./discountsyearly(ii);
end

disp('--- 2. EXERCISE: BERMUDIAN SWAPTION PRICING VIA HULL-WHITE ---')
disp('2a. Price the Bermudan Swaption')

tic % Start timer
% Compute the Bermudan option's price and the European Swaption prices:
[BermudanPrice,EuropeanSwaptionPrices] = buildTrinomialTree(T_start,K,Maturity,N,a,sigma,dates,discounts,fwdD_matrix,deltas,1);
elapsed_time = toc; % Stop timer
% Display the result:
fprintf('\nFor N=%d steps each year\nComputational time: %d\nBermudan Swaption price: %.8f%% \n\n',N,elapsed_time,100*BermudanPrice)

%% b. Check the implementation
% Check that the implementation of the tree is correct by evaluating
% European swaptions prices and confronting the results with the ones
% obtained via the Jamshidian formula.
disp('2b. Check the implementation')
% Compute the prices of the European swaptions via the Jamshidian formula:
Put_Prices = JamshidianPrices(a,sigma,Maturity,0,T_start,fwdD_matrix,K,deltas,discountsyearly,1);
% Compute the error
Error = (abs(EuropeanSwaptionPrices-Put_Prices'))';
% Display the results:
disp('----------------------------------------------------------------------------------------------------')
fprintf('Expiry |  European Swaption Value via tree |  Put option price via Jamshidian formula |  Error\n')
disp('----------------------------------------------------------------------------------------------------')
for ii=1:length(EuropeanSwaptionPrices)
    fprintf('%d      |  %.8f%%                      |  %.8f%%                             |  %d\n',ii+1, 100*EuropeanSwaptionPrices(ii), 100*Put_Prices(ii), Error(ii))
end
disp(' ')
%% c. Define upper and lower bound
% As an upper bound, the sum of the European options is taken, as a lower
% bound the max among the European options is taken:
disp('2c. Define upper and lower bound')
fprintf('Upper bound for the Bermudan: %.8f%%\n', sum(100*Put_Prices))
fprintf('Lower bound for the Bermudan: %.8f%%\n\n', max(100*Put_Prices))

%% Evaluation of the error for different values of N
disp('2 EXTRA: Evaluation of the error for different values of N')
% Initialization of the vectors:
Error = [];
BermudanPrice = zeros(10,1);
for N=2.^(1:10)
    tic % Start timer
    % Compute the Bermudan Swaption price and the European options prices:
    [BermudanPrice(log2(N)),EuropeanSwaptionPrices] = buildTrinomialTree(T_start,K,Maturity,N,a,sigma,dates,discounts,fwdD_matrix,deltas,0);
    elapsed_time = toc; % Stop timer
    % Update the errors matrix
    Error = [Error abs(EuropeanSwaptionPrices-Put_Prices')'];
    % Display the results:
    fprintf('\nFor N=%d steps each year\nComputational time: %d\nBermudan Swaption price: %.8f%% \n', N, elapsed_time, 100*BermudanPrice(log2(N)))
    fprintf('Bounds for the Bermudan: [ %.8f%% , %.8f%% ]\n', max(100*Put_Prices), sum(100*Put_Prices))
    disp('-----------------------------------------------------------------------------------------------------')
    fprintf('Expiry |  European Swaption Value via tree |  Put option price via Jamshidian formula |  Error\n')
    disp('-----------------------------------------------------------------------------------------------------')
    for ii=1:length(EuropeanSwaptionPrices)
        fprintf('%d      |  %.8f%%                      |  %.8f%%                             |  %d\n', ii+1, 100*EuropeanSwaptionPrices(ii), 100*Put_Prices(ii), Error(ii))
    end
    disp('-----------------------------------------------------------------------------------------------------')
end

% Plot of the error between the prices via the tree and the prices obtained
% with the Jamshidian formula for different N values:
figure()
loglog(2.^(1:10),max(Error),'--o',LineWidth=2)    % Plot the error
hold on
loglog(2.^(1:10),2.^(-1:-1:-10),'--',LineWidth=2) % Plot 1/N
grid on
xlabel('Number of time steps each year (N)',FontSize=24)
ylabel('Error',FontSize=24)
title('Check for the tree precision',FontSize=28)
legend('Error','1/N',FontSize=28)
hold off

% Plot the Bermudan price computed for different values of N:
figure()
plot(2.^(1:10),100*BermudanPrice,'--o',LineWidth=2) % Plot the Bermudan price
grid on
xlabel('Number of time steps each year (N)',FontSize=24)
ylabel('Bermudan Price (%)',FontSize=24)
title('Bermudan price for different number of time steps',FontSize=28)
