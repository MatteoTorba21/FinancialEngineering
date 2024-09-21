function [] = plotCorrelationRho(dates,discounts,TTM,R1,R2,datesCDS,intensities_exact_1,intensities_exact_2)
% Function to plot the default probabilities (probability only 1 name default, probability only 2 name 
% default, proability of both defaults) for different correlations
%
% INPUTS:
% dates:                 Vector of dates obtained from the bootstrap
% discounts:             Vector of discounts obtained from the bootstrap
% TTM:                   Time to maturity expressed in years
% R1:                    Recovery rate of the first obligor
% R2:                    Recovery rate of the second obligor
% datesCDS:              Vector of dates of the Credit Default Swap
% intensities_exact_1:   Vector of the exact intensities of the first obligor
% intensities_exact_2:   Vector of the exact intensities of the second obligor
%

N = 1000;                           % Number of simulations
rho_vector = [-.99 -.5 0 .5  .99];  % Vector of correlations
for ii = 1:length(rho_vector)
    [s,CI,P] = firstToDefault2Names(dates,discounts,TTM,R1,R2,datesCDS,intensities_exact_1,intensities_exact_2,rho_vector(ii),N);
    P1(ii) = P(1,2);
    P2(ii) = P(2,1);
    P3(ii) = P(2,2);
end
figure()
plot(rho_vector,P1,'r*','LineWidth', 8)
hold on
plot(rho_vector,P2,'bo','LineWidth', 8)
plot(rho_vector,P3,'gs','LineWidth', 8)
plot(rho_vector,P1,'r','LineWidth', 2)
plot(rho_vector,P2,'b','LineWidth', 2)
plot(rho_vector,P3,'g','LineWidth', 2)

legend('Probability of only ISP default', 'Probability of only UCG default', 'Probability of both defaults','FontSize',18)
title('Probabilities vs Rho', 'FontSize', 24)
xlabel('Rho', 'FontSize', 18);
ylabel('Probabilities (%)','FontSize', 18);

end % function plotCorrelationRho