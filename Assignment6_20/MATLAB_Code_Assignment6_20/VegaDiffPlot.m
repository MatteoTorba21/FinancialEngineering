function vegaDiff = VegaDiffPlot(VegaBucket,TotalVega,YTM)
%VEGADIFFPLOT firstly computes the difference between the total vega and the
%sum of all vega-bucket sensitivities; then it plots the cumulative sum of the bucket vegas against the total  
% INPUT:
% VegaBucket:       Vega-Bucket sensitivities
% TotalVega:        Total-Vega value
% YTM:              quoted maturities
%
% OUTPUT:  
% diff:             difference between sum of all vega-bucket sensitivities and total vega


% computing the difference

vegaDiff = TotalVega - sum(VegaBucket);


% Plot 
figure
hold on
plot(YTM,cumsum(VegaBucket), '-o', 'LineWidth',2, 'DisplayName', 'Cumulative Vega-Bucket sum');
plot(YTM, TotalVega*ones(length(YTM),1), '-','LineWidth',2, 'DisplayName','Vega-Total');
legend('show',Location='west',FontSize = 20)
title ('Vega-Bucket VS Total-Vega',FontSize=20)
xlabel('Dates',FontSize=20);
ylabel('Vegas',FontSize=20);
hold off
end

