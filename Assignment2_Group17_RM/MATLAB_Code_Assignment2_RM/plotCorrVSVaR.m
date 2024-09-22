function plotCorrVSVaR(rho,VaR_Def_200,VaR_DM_200,VaR_Def_20,VaR_DM_20)
% Function to plot the obtained results for different values of correlation

figure(1)
plot(rho,VaR_Def_200,'--*',LineWidth=1.5)
hold on
plot(rho,VaR_DM_200,'--o',LineWidth=1.5)
plot(rho,VaR_Def_20,'--*',LineWidth=1.5)
plot(rho,VaR_DM_20,'--o',LineWidth=1.5)
legend('VaR Default only for 200 names','VaR Default and migration for 200 names','VaR Default only for 20 names','VaR Default and migration for 20 names','FontSize',18,'Location','northwest')
xlabel('ρ','FontSize',18)
ylabel('VaR','FontSize',18)
title('Plot the computed VaR for different levels of correlation ρ','FontSize',24)
hold off


end % function plotCorrVSVaR