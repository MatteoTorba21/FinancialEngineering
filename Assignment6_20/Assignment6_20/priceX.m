function x = priceX(dates, SwapExpiries, rates, Strikes, spotvol)
% Function to determine the percentage value (wrt the notional) of the
% upfront in the contract given in the text via the spot volatilities: the
% upfront is determine by passing from the falt volatilities to the spot
% volatilities.
%
% INPUTS:
% dates:        Struct of dates of the financial instruments quoted in the market
% SwapExpiries: Vector of the quoted swap expiries
% rates:        Struct of rates of the financial instruments quoted in the market
% Strikes:      Vector of strikes of quoted volatilities
% spotvol:      Matrix of spot volatilities
%
% OUTPUT:
% x:            Upfront's value of the contract

% set plot input to false if not specified
% if nargin < 7
%         PlotVolSmile = 0;
% end
% Yearly swap expiries up to 50 years:
FullSwapExpiries = [1:50]';
% Bootstrap:
[dates_bootstrap, discounts_bootstrap] = interpolateAndLaunchBootsrap(dates,rates,SwapExpiries,FullSwapExpiries);

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
% for ii = 1:length(YTM)
%     capsPrice(ii,:) = capPrice(discounts3m,FLDates,deltas,Strikes,Caps_vol(ii,:)',YTM(ii));
% end

% Calibrate the spot volatilites:
% spotvol = [];
% for ii = 1:length(Strikes)
%     spotvol = [spotvol calibratevol(Strikes(ii),capsPrice(:,ii),Caps_vol(1,ii),YTM,FLDates,discounts3m,deltas,ACT365)];
% end
% 
% % Plot the spot volatilities smile
% if PlotVolSmile==1
%     for ii=1:length(YTM)
%         figure(ii)
%         plot(Strikes,Caps_vol(ii,:),'--o', LineWidth=2)
%         hold on
%         plot(Strikes,spotvol(ii,:),'--o', LineWidth=2)
%         legend('Flat volatilities','Spot volatilities')
%         xlabel('Strikes')
%         ylabel('Volatilities')
%         title(['Flat VS Spot volatilities for YTM = ', num2str(YTM(ii))])
%         hold off
%     end
%     figure();
%     [X, Y] = meshgrid(Strikes, FLDates(1:YTM(end)*4-1));  % Create a grid of strikes and dates
%     surf(X, Y, spotvol);  % Plot the surface
%     xlabel('Strike');
%     ylabel('Date');
%     zlabel('Spot Volatility');
%     title('Volatility Surface');
%     hold on;
%     [X, Y] = meshgrid(Strikes, FLDates(YTM*4+1));  % Create a grid of strikes and dates
%     plot3(X, Y, Caps_vol, 'o', 'MarkerFaceColor', 'r');
%     legend('Spot vol surface','Quoted flat vol')
% 
% end

% Compute the forward discounts:
fwd_discounts = discounts3m(2:end)./discounts3m(1:end-1);
% fwd_discounts = fwd_discounts(2:ncapl+1); % select only the needed forward discounts
% Compute the forward Libor:
L = (1./fwd_discounts -1)./deltas;
% Strikes of the structured cap:
K = [0.032, 0.035, 0.04];
% Number of caplets to consider for the structured cap:
ncaplets = [19, 20, 20];

% Price of the structured cap computed via the spot volatlities:
Cap = PriceStructuredCaps(K, Strikes, L, deltas, FLDates, spotvol, ncaplets, discounts3m);
% BPV:
BPV = discounts3m(3:61)'*deltas(2:60);
% Net present value of the contract:
NPV_Spot = 0.009*BPV + 1 - discounts3m(2) - 0.01*discounts3m(2)*deltas(1) + Cap;
% Value of the upfront:
x = NPV_Spot;

end