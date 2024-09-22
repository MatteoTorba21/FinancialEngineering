function Vega = BucketVegaCap(dates, SwapExpiries, rates, Strikes, Caps_vol, YTM,CapStrikes,CapMaturity,spotvol)
% Function to compute the Vega-bucket of Caps by shifting the market
% volatility of 1 bp
% 
% INPUTS: 
% dates:        Struct of dates of the financial instruments quoted in the market     
% SwapExpiries:     Vector of the quoted swap expiries
% rates:        Struct of rates of the financial instruments quoted in the market
% Strikes:      Vector of strikes of quoted volatilities
% Caps_vol:     Matrix of quoted volatilities
% YTM:          Maturities of the market volatilities
% CapStrikes:   Vector of cap's strikes
% CapMaturity:  Vector of cap's maturities
% spotvol:      Matrix of spot volatilities
%
% OUTPUTS:
% Vega:         Matrix of Vega-Bucket of the Caps (rows: different bucket dates, columns: different strikes and maturities of caps)

% Initialization of the output:
Vega = zeros(size(Caps_vol,1), length(CapMaturity));
% Iterate on the years
for ii = 1:size(Caps_vol,1)
    % Flat vols inizialization
    NewCapsVol = Caps_vol;
    % Shift vols for each bucket iteratively
    NewCapsVol(ii, :) = NewCapsVol(ii, :) + 1;
    % Yearly swap expiries up to 50 years:
    FullSwapExpiries = [1:50]';
    % Bootstrap:
    [dates_bootstrap, discounts_bootstrap] = interpolateAndLaunchBootsrap(dates,rates,SwapExpiries,FullSwapExpiries);
    % Calibration of the spot volatility:
    NewCapsSpotVol = launchCalibration(dates,dates_bootstrap,discounts_bootstrap,Strikes,NewCapsVol,YTM,0);
    for k = 1:length(CapMaturity)
    % Compute the upfront with shifted vols by bucket
    Vega(ii,k) = priceCapViaSpotVol(dates, SwapExpiries, rates, Strikes, CapStrikes(k),CapMaturity(k),NewCapsSpotVol)- priceCapViaSpotVol(dates, SwapExpiries, rates, Strikes, CapStrikes(k),CapMaturity(k), spotvol);
    end

end
end