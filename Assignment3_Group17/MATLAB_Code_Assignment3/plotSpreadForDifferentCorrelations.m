function [] = plotSpreadForDifferentCorrelations(dates,discounts,TTM,R1,R2,datesCDS,intensities_exact_1,...
    intensities_exact_2,rho_vector,N,CDSSpreads_1,CDSSpreads_2)
% Function to plot the spreads of the First to Default with different correlations
%
% INPUTS:
% dates:                   Vector of dates obtained from the bootstrap
% discounts:               Vector of discounts obtained from the bootstrap
% TTM:                     Time to maturity expressed in years
% R1:                      Recovery rate of the first obligor
% R2:                      Recovery rate of the second obligor
% datesCDS:                Vector of dates of the Credit Default Swap
% intensities_exact_1:     Vector of the exact intensities of the first obligor
% intensities_exact_2:     Vector of the exact intensities of the second obligor
% rho:                     Vector of Default correlation of the two obligors
% N:                       Number of Monte-Carlo simulations
% CDSSpreads_1:            Vector of spreads of the Credit Default Swap of the first obligor
% CDSSpreads_2:            Vector of spreads of the Credit Default Swap of the first obligor
%

% Initialization of vectors of spreads, bid and ask:
s = zeros(size(rho_vector));
bid = zeros(size(rho_vector))';
ask = zeros(size(rho_vector))';

for ii = 1:length(rho_vector)
    [s(ii),CI] = firstToDefault2Names(dates,discounts,TTM,R1,R2,datesCDS,intensities_exact_1,intensities_exact_2,rho_vector(ii),N);
    bid(ii) = CI(1);
    ask(ii) = CI(2);
end
figure()
plot(rho_vector,s*1e4,'-*','LineWidth',4)
hold on
plot(rho_vector,bid*1e4,'LineWidth',2,'LineStyle','--')
plot(rho_vector,ask*1e4,'LineWidth',2,'LineStyle','--')
plot(rho_vector, ones(size(rho_vector)) * CDSSpreads_1(TTM)*1e4,'LineWidth',2,'LineStyle','--')
plot(rho_vector, ones(size(rho_vector)) * CDSSpreads_2(TTM)*1e4,'LineWidth',2,'LineStyle','--')
xfill = [rho_vector; flipud(rho_vector)];
yfill = [ask'; flipud(bid')];
fill(xfill, yfill , 'b', 'FaceAlpha', 0.3);

legend('Spread of FTD','Spread of FTD bid','Spread of FTD ask', 'Spread of ISP CDS' , 'Spread of UCG CDS', 'Location','southeast','FontSize',18);
title('FTD Spread vs Correlation','FontSize', 24)
xlabel('Correlation', 'FontSize', 18);
ylabel('Spread (in bps)', 'FontSize', 18);

end % function plotSpreadForDifferentCorrelations