function [sigma, k, nu] = NMMCalibration(alpha)
% Function to run the calibration of the NIG pameters
%
% INPUT:
% alpha:    parameter in the Laplace exponent
%
% OUTPUT:
% sigma:    average volatility
% k:        volvol
% nu:       volatility skew

load('cSelect20230131_B.mat')
TTM = 1.0027397;            % Time to maturity
discount_TTM = 0.961345775982266;   % Discount 1 y
K = cSelect.strikes;        % Strikes
MktVol = cSelect.surface;   % Quoted volatilities
Spot = cSelect.reference;   %Spot price
dividend = cSelect.dividend;% dividend yield
F0 = Spot*exp(-dividend*TTM)/discount_TTM;  % Fwd price at time 0

% FFT parameters and definition of price function w.r.t. moneynesses and
% Model parameters definition
M = 15;
N = 2^M;
dx = 0.0025;
% Computation of the remaining model parameters
x1 = -dx*(N-1)/2;
x = [x1:dx:-x1];
mkt_moneyness = log(F0./K);
% Market prices computation via Black formula 
MktPrice = blkprice(F0, K, -log(discount_TTM)/TTM, TTM,MktVol);
% Definition of the integrand function
lnL = @(w,params) TTM./params(2)*(1-alpha)/alpha*(1-(1+(w.*params(2)*params(1)^2)/(1-alpha)).^alpha);
Phi = @(x,params) exp(-1i*x.*lnL(params(3),params) + lnL((x.^2 + 1i*(1+2*params(3)).*x)./2,params));
f = @(x,params) 1/(2*pi).*Phi(-x-1i/2,params).*1./(x.^2 + 1/4);
I = @(parameters) integralViaFFT(@(x) f(x,parameters),M,dx);
% Call prices computed in quoted moneynesses in function of parameters
C = @(parameters) Spot.*exp(-dividend*TTM).*(1-exp(-x/2).*I(parameters));
modelPricer = @(parameters) interp1(x,C(parameters),mkt_moneyness);
% lsqnonlin options
options = optimoptions('lsqnonlin','MaxFunctionEvaluations',300000','MaxIterations',100000,'OptimalityTolerance', 1e-15,'Display', 'off');
% Least squares optimization to find the optimal model parameters (calibration peformed on call prices using L2 as metric)
calibratedParams = lsqnonlin(@(params) ((modelPricer(params)-MktPrice)),[0.2,1,0.2],[0,0,-10],[2,20,20],options);

% Define the outputs:
sigma = calibratedParams(1);
k = calibratedParams(2);
nu = calibratedParams(3);
end
