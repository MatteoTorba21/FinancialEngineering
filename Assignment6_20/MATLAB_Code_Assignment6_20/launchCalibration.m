function spotvol = launchCalibration(dates,dates_bootstrap,discounts_bootstrap,Strikes,Caps_vol,YTM,PlotVolSmile)
% Function to calibrate the spot volatilities from the flat volatilities,
% by imposing that they give the same Caps prices.
% 
% INPUTS:
% dates                 Struct of dates of the financial instruments quoted in the market
% dates_bootstrap       Vector of dates obtained with the bootstrap
% discounts_bootstrap   Vector of discounts obtained with the bootstrap
% Strikes               Vector of strikes of quoted volatilities
% Caps_vol              Matrix of volatilities
% YTM                   Maturities of the market volatilities
% PlotVolSmile:         Flag variable
%                         -If PlotVolSmile=1 plot the spot volatilities smile
%                         -If PlotVolSmile=0 do not plot the volatilities smile
%
% OUTPUTS:
% spotvol:      Matrix of spot volatilities

% Set plot input to false if not specified
if nargin < 7
        PlotVolSmile = 0;
end

% Determine the upfront via the spot vols
% Start by computing the Caps price for the volatility grid
% Find the payment dates of the Caps (every 3 months)
FLDates = findFloatingLegDates(datenum(dates.settlement),50, eurCalendar);
FLDates = [dates.settlement; FLDates];  % add the settlement date
% Interpolate the discounts at payment dates:
discounts3m = InterpDFviaRates(dates_bootstrap,discounts_bootstrap,FLDates);
% yearfrac Act/360:
ACT360 = 2;
% yearfrac Act/365:
ACT365 = 3;
% Compute the year fractions:
deltas = yearfrac(FLDates(1:end-1),FLDates(2:end),ACT360);

% Fill the matrix of the Cap prices for the different strikes and
% volatilities:
for ii = 1:length(YTM)
    capsPrice(ii,:) = capPrice(discounts3m,FLDates,deltas,Strikes,Caps_vol(ii,:)',YTM(ii));
end

% Calibrate the spot volatilites:
spotvol = [];
for ii = 1:length(Strikes)
    spotvol = [spotvol calibratevol(Strikes(ii),capsPrice(:,ii),Caps_vol(1,ii),YTM,FLDates,discounts3m,deltas,ACT365)];
end

% Plot the spot volatilities smile
if PlotVolSmile==1
    for ii=1:length(YTM)
        figure(ii)
        plot(Strikes,Caps_vol(ii,:),'--o', LineWidth=2)
        hold on
        plot(Strikes,spotvol(ii,:),'--o', LineWidth=2)
        legend('Flat volatilities','Spot volatilities')
        xlabel('Strikes')
        ylabel('Volatilities')
        title(['Flat VS Spot volatilities for YTM = ', num2str(YTM(ii))])
        hold off
    end
    figure();
    [X, Y] = meshgrid(Strikes, FLDates(1:YTM(end)*4-1));  % Create a grid of strikes and dates
    surf(X, Y, spotvol);  % Plot the surface
    xlabel('Strike');
    ylabel('Date');
    zlabel('Spot Volatility');
    title('Volatility Surface');
    hold on;
    [X, Y] = meshgrid(Strikes, FLDates(YTM*4+1));  % Create a grid of strikes and dates
    plot3(X, Y, Caps_vol, 'o', 'MarkerFaceColor', 'r');
    legend('Spot vol surface','Quoted flat vol')
end

end