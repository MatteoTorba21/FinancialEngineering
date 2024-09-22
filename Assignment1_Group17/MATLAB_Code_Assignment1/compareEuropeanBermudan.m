function [] = compareEuropeanBermudan(S0,K,r,TTM,sigma,dmin,dmax,N,EEfreq,Noptions)
%
%INPUT
% S0:        underlying price
% K:         strike
% r:         zero-rate assumed constant
% TTM:       time-to-maturity
% sigma:     volatility
% dmin:      minimum dividend Yield
% dmax:      maximum dividend Yield
% N:         lenght of dividend discretization
% EEFreq:    early execise frequency 

BermudianStepsCRR = 1e3; % number of intervals of each period in CRR tree for Bermudian option

dividends = linspace(dmin,dmax,N); % row vector of dividend Yields
B = exp(-r*TTM);                   % discount factor
F0 = S0*exp(-dividends*TTM)/B;     % row vector of forwards computend with dividend Yields
EuroPrices = zeros(size(F0));
BermudanPrices = zeros(size(F0));
% for each value of the forward we compute the European option price and
% the Bermudian option price 
for i = 1 : length(F0)
     EuroPrices(i) = EuropeanOptionPrice(F0(i),K,B,TTM,sigma,1,0,1);
     BermudanPrices(i) = BermudanOptionCRR(F0(i),K,r,TTM,EEfreq,sigma,BermudianStepsCRR);
end

NotionalEuroPrices = Noptions.*EuroPrices;
NotionalBermudanPrices = Noptions.*BermudanPrices;

% we compare the trend of the European option price with the trend of the
% Bermudian option price (both with dividend Yields)
figure()
grid on
plot(dividends,NotionalEuroPrices,'r','LineWidth',2);
hold on
plot(dividends,NotionalBermudanPrices,'b','LineWidth',2);
xlabel('Dividend yield', 'FontSize',18);
legend('European','Bermudan','FontSize',24);
hold off

end % function compareEuropeanBermudan
