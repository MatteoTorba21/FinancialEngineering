clear all
load('cSelect20230131_B.mat')
TTM = 1.0027397;    % Time to maturity
discount_TTM = 0.961345775982266;   % Discount 1 y
K = cSelect.strikes;    % Strikes
MktVol = cSelect.surface;   % Quoted volatilities
Spot = cSelect.reference;   %Spot price
dividend = cSelect.dividend;    % dividend yield
F0 = Spot*exp(-dividend*TTM)/discount_TTM;  % Fwd price at time 0

% FFT parameters and definition of price function w.r.t. moneynesses and
% model parameters
M = 15;
N = 2^M;
dx = 0.0025;
x1 = -dx*(N-1)/2;
x = [x1:dx:-x1];
mkt_moneyness = log(F0./K);
MktPrice = blsprice(Spot,K,-log(discount_TTM)/TTM,TTM,MktVol,dividend);
MktPrice = blkprice(F0, K, -log(discount_TTM)/TTM, TTM,MktVol);
alpha = 1/3;
lnL = @(w,params) TTM./params(2)*(1-alpha)/alpha*(1-(1+(w.*params(2)*params(1)^2)/(1-alpha)).^alpha);
Phi = @(x,params) exp(-1i*x.*lnL(params(3),params) + lnL((x.^2 + 1i*(1+2*params(3)).*x)./2,params));
f = @(x,params) 1/(2*pi).*Phi(-x-1i/2,params).*1./(x.^2 + 1/4);
I = @(parameters) integralViaFFT(@(x) f(x,parameters),M,dx);
% Call prices computed in quoted moneynesses in function of parameters
C = @(parameters) Spot.*exp(-dividend*TTM).*(1-exp(-x/2).*I(parameters));
modelPricer =@(parameters) interp1(x,C(parameters),mkt_moneyness);
% lsqnonlin options
options = optimoptions('lsqnonlin','MaxFunctionEvaluations',30000','MaxIterations',10000,'OptimalityTolerance', 1e-8,'Display', 'off');
% Least squares optimization to find the optimal model parameters (calibration peformed on call prices using L2 as metric)
calibratedParams = lsqnonlin(@(params) ((modelPricer(params)-MktPrice)),[0.2,1,0.2],[0,0,-10],[2,10,10],options);
% Prices with calibrated parameters
modelPrices = modelPricer(calibratedParams);
% Plot of call prices w.r.t. moneyness
figure()
plot(mkt_moneyness,modelPrices,'.', 'LineWidth',2,'MarkerSize',23)
grid on
hold on
plot(mkt_moneyness,MktPrice,"s",'MarkerSize',14,'LineWidth',2)
xlabel('Moneyness',Fontsize=28)
ylabel('Call price',Fontsize=28)
legend('Model prices', 'Quoted prices', 'Location', 'northwest', 'FontSize', 24)
title(['Quoted and model call prices w.r.t. moneyness'],FontSize=28);
% Plot of call prices w.r.t. strikes
figure()
plot(K,modelPrices,'.', 'LineWidth',2,'MarkerSize',23)
grid on
hold on
plot(K,MktPrice,"s",'MarkerSize',14,'LineWidth',2)
xlabel('Strike',Fontsize=28)
ylabel('Call price',Fontsize=28)
legend('Model prices','Quoted prices',FontSize=24)
title(['Quoted and model call prices w.r.t. strike price'],FontSize=28);
calibrationError = sum((modelPrices-MktPrice).^2)/length(K);
% Black implied volatility for model prices
modelVol = blkimpv(F0, K, -log(discount_TTM)/TTM, TTM, modelPrices);
% Plot of volatility smiles
figure()
plot(K,MktVol,'.', 'MarkerSize',23,'LineWidth',4)
grid on
hold on
plot(K,modelVol, 's','MarkerSize',15,'LineWidth',2)
xlabel('Strike',Fontsize=28)
ylabel('Volatility',Fontsize=28)
legend('Model implied volatility','Quoted volatility',FontSize=24)
title(['Quoted and model volatility smiles'],FontSize=28);
figure()
[ax, h1, h2] = plotyy(K, MktPrice-modelPrices, K, MktVol-modelVol);
xlabel('Strike')
ylabel('error')
xlabel(ax(1), 'Strike',Fontsize=28);
ylabel(ax(1), 'Price error',Fontsize=28);
ylabel(ax(2), 'Volatility error',Fontsize=28);
title(['Model price and volatility error'],FontSize=28);
set(h1, 'LineStyle', 'none');
set(h2, 'LineStyle', 'none');
set(h1,'Marker','.','MarkerSize',20);
set(h2,'Marker', '*','MarkerSize',7,'LineWidth',2);
legend('Prices difference','Volatility difference',FontSize=24)
grid(ax(1), 'on');
grid(ax(2), 'on');

disp('--- 4. CASE STUDY: VOLATILITY SURFACE CALIBRATION --- ')
fprintf(['Calibrated parameters with α = 1/3: \n' ...
    'σ (average volatility):  %f\n' ...
    'η (volatility skew): %f\n' ...
    'k (vol of vol): %f\n'], calibratedParams(1), calibratedParams(3), calibratedParams(2));
fprintf('MSE of prices calibration: %f\n', calibrationError);
fprintf('MSE of volatility: %f\n', sum(modelVol-MktVol).^2)/length(K);
