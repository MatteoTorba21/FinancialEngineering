function [] = plotBootstrapCDSResults(settlement,dates,discounts,datesCDS,CDSSpreads,R,intensities_approx,intensities_exact)
% function to plot of the JT intensities VS exact and approx intensities for different maturities
%
% INPUT:
% settlement:           Settlement day
% dates:                Dates obtained from the bootstrap
% discounts:            Discounts obtained from the bootstrap
% datesCDS:             Vector of dates of the Credit Default Swap
% CDSSpreads:           Vector of spreads of the Credit Default Swap
% R:                    Recovery value
% intensities_approx:   Approximated Intensities
% intensities_exact:    Exact Intensisties (with the accrual)
%

% Initialization of the vector of intensities of JT for different maturities:
intensities_JT = zeros(length(datesCDS),1);

% Plot of the JT intensities VS exact and approx intensities for different maturities:
figure()
for ii = 1:length(datesCDS)-1
    [~,~,intensities_JT_iter] = bootstrapCDS(dates,discounts,datesCDS(1:ii),CDSSpreads(1:ii),3,R);
    intensities_JT(ii) = intensities_JT_iter(1);
    xx = linspace(settlement,datesCDS(ii),1000);
    subplot(2,3,ii)
    if ii == 1
        grid on
        plot(xx,ones(size(xx))*intensities_approx(1)*1e4,LineWidth=2)
        hold on
        plot(xx,ones(size(xx))*intensities_exact(1)*1e4,LineWidth=2);
        plot(xx,ones(size(xx))*intensities_JT_iter*1e4,LineWidth=2);
    else
        grid on
        plot(xx,interp1(datesCDS(1:ii), intensities_approx(1:ii)*1e4, xx,'next','extrap'),LineWidth=2)
        hold on
        plot(xx,interp1(datesCDS(1:ii), intensities_exact(1:ii)*1e4, xx,'next','extrap'),LineWidth=2)
        plot(xx,interp1(datesCDS(1:ii), intensities_JT_iter*1e4, xx,'next','extrap'),LineWidth=2)
    end
    legend('Approx','Exact','JT','Location','northwest','FontSize',14)
    ylabel('Intensities (in bps)','FontSize',16);
    xlabel('Date','FontSize',16);
    title(['Intensities Over Time up to ',num2str(ii),'y'],'FontSize',18);

    set(gca, 'XTick', [])
    xlim([settlement, datesCDS(ii)]);
    if ii == 1
        ylim([intensities_JT_iter*1e4-10, intensities_JT_iter*1e4+10]);
    end

    hold off
end

% Plot of the JT intensities VS exact and approx intensities for the whole CDS curve:
figure()
[~,~,intensities_JT_iter] = bootstrapCDS(dates,discounts,datesCDS,CDSSpreads,3,R);
intensities_JT(7) = intensities_JT_iter(1);
grid on
xx = linspace(settlement,datesCDS(end),1000);
plot(xx,interp1(datesCDS, intensities_approx*1e4, xx,'next','extrap'),LineWidth=2)
hold on
plot(xx,interp1(datesCDS, intensities_exact*1e4, xx,'next','extrap'),LineWidth=2)
plot(xx,interp1(datesCDS, intensities_JT_iter*1e4, xx,'next','extrap'),LineWidth=2)

legend('Approx','Exact','JT','Location','northwest','FontSize',18)
ylabel('Intensities (in bps)','FontSize',18);
xlabel('Date','FontSize',18);
title('Intensities Over Time','FontSize',24);

% rotate x-axis tick labels for better readability (optional)
xtickangle(45);
ax = gca;  % get the current axes
ax.XTick = datesCDS(1:10:end);  % show every 10th date on the x-axis
datetick('x', 'dd-mm-yyyy'); % format the x-axis as dates
xlim([settlement, datesCDS(end)]);

hold off

% Numerical computation of the mean lambdas:
mean_values_exact = zeros(length(intensities_exact),1);
for ii=1:length(intensities_exact)
    mean_values_exact(ii) = mean(intensities_exact(1:ii));
end
% Difference between the mean exact lambdas and JT approximation:
differences_JTVSmean = mean_values_exact-intensities_JT;

disp('Mean values of the exact intenities VS JT approximation up to different maturities (results in bps):')
fprintf('for maturity in 1 year:  λ_{mean}=%.1f    λ_{JT}=%.1f    difference=%.1f\n',mean_values_exact(1)*1e4,intensities_JT(1)*1e4,differences_JTVSmean(1)*1e4)
fprintf('for maturity in 2 years: λ_{mean}=%.1f    λ_{JT}=%.1f    difference=%.1f\n',mean_values_exact(2)*1e4,intensities_JT(2)*1e4,differences_JTVSmean(2)*1e4)
fprintf('for maturity in 3 years: λ_{mean}=%.1f    λ_{JT}=%.1f    difference=%.1f\n',mean_values_exact(3)*1e4,intensities_JT(3)*1e4,differences_JTVSmean(3)*1e4)
fprintf('for maturity in 4 years: λ_{mean}=%.1f    λ_{JT}=%.1f    difference=%.1f\n',mean_values_exact(4)*1e4,intensities_JT(4)*1e4,differences_JTVSmean(4)*1e4)
fprintf('for maturity in 5 years: λ_{mean}=%.1f    λ_{JT}=%.1f    difference=%.1f\n',mean_values_exact(5)*1e4,intensities_JT(5)*1e4,differences_JTVSmean(5)*1e4)
fprintf('for maturity in 6 years: λ_{mean}=%.1f    λ_{JT}=%.1f    difference=%.1f\n',mean_values_exact(6)*1e4,intensities_JT(6)*1e4,differences_JTVSmean(6)*1e4)
fprintf('for maturity in 7 years: λ_{mean}=%.1f    λ_{JT}=%.1f    difference=%.1f\n\n',mean_values_exact(7)*1e4,intensities_JT(7)*1e4,differences_JTVSmean(7)*1e4)

end % function plotBootstrapCDSResults