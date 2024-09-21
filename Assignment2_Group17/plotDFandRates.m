function plotDFandRates(dates, discounts)
%Input:
% dates:      vector of dates (output from bootstrap)  
% discounts:  vector of discount factors (output from bootstrap)

%discount curve
figure()
grid on
yyaxis left;
plot(dates,discounts,'-o')
ylabel('Discount Factor','FontSize',18);

% rotate x-axis tick labels for better readability (optional)
xtickangle(45);
%zero-rates
yyaxis right;
plot(dates(2:end),zeroRates(dates,discounts)*100,'-x')
ylabel('Zero Rate','FontSize',18)
xlabel('Date','FontSize',18);
legend('Discount Factors','Zero Rates', 'FontSize',24);
title('Discount Factors and Rates Over Time','FontSize',24);
% rotate x-axis tick labels for better readability (optional)
xtickangle(45);
ax = gca;  % get the current axes
ax.XTick = dates(1:10:end);  % show every 10th date on the x-axis
datetick('x', 'dd-mm-yyyy'); % format the x-axis as dates
xlim([dates(1), dates(end)]);

end % function plotDFandRates