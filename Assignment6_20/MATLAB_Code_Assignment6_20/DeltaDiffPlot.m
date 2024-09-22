function diff = DeltaDiffPlot(rates, dates, SwapExpiries, Strikes, X, DeltaSensitivities, DeltaDates, spotvol)
% SENSPLOT plot cumulated sum of Delta-Bucket sensitivities with respect to
% Total Delta, and compute the difference between Total Delta and total sum
% of the Delta Buckets
% INPUTS:
% dates:        Struct of dates of the financial instruments quoted in the market
% SwapExpiries: Vector of the quoted swap expiries
% rates:        Struct of dates of the financial instruments quoted in the market
% Strikes:      Vector of strikes of quoted volatilities
% X:            Upfront value
% DeltaSensitivities: Delta-Bucket sensitivities
% DeltaDates:   Dates of quoted rates
% spotvol:      Matrix of spot volatilities
%
% OUTPUT: 
% diff: Total Delta - cumulated sum of Delta Buckets

% Shifted Rates
allRatesShifted = rates;
allRatesShifted.depos = allRatesShifted.depos + 1e-4;
allRatesShifted.futures = allRatesShifted.futures + 1e-4;
allRatesShifted.swaps = allRatesShifted.swaps + 1e-4;

% Total delta as the difference of the upfront prices
X_totalDelta = priceX(dates, SwapExpiries, allRatesShifted, Strikes, spotvol) - X;
disp(' ')
fprintf('Δ-Total: %.10f%\n',X_totalDelta);
diff = X_totalDelta-sum(DeltaSensitivities);
disp(' ')
fprintf('Δ-Total - sum of Δ buckets: %.10f%\n',diff);
disp(' ')

% Plot 
figure
hold on
plot(DeltaDates,cumsum(DeltaSensitivities), '-o', 'LineWidth',2, 'DisplayName', 'Cumulative Delta-Bucket sum');
plot(DeltaDates, X_totalDelta*ones(length(DeltaDates),1), '-','LineWidth',2, 'DisplayName','Delta-Total');
legend('show',Location='SouthEast',FontSize = 20)
title ('Delta-Bucket VS Total-Delta',FontSize = 25)
xlabel('Dates',FontSize=20);
ylabel('Deltas',FontSize=20);
hold off
end

