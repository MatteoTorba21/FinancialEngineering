function  DeltaCap = DeltaBucketCap(dates, quotedSwaps, rates, Strikes, CapMaturity, CapStrike, spotvol)
% Function computes the Delta Bucket of Caps, shifting rates one at a time
% for the given strikes and maturities
% 
% INPUTS:
% dates:        Struct of dates of the financial instruments quoted in the market
% quotedSwaps:  Vector of the quoted swap expiries
% rates:        Struct of rates of the financial instruments quoted in the market
% Strikes:      Vector of strikes of quoted volatilities
% CapMaturity:  Vector of cap's maturities
% CapStrike:    Vector of cap's strikes
% spotvol:      Matrix of spot volatilities
%
%OUTPUT:
%DeltaCap:      Matrix of Delta Bucket in which the rows are the different
%maturities and the coloumn the different Cap evaluated, in terms of
%maturity and strike

% Depos Delta bucket initialization
Delta_depos = zeros(4,length(CapStrike));
for ii = 1:4
    % Rates struct inizialization
    newrates = rates;
    % Increase each bucket iteratively
    newrates.depos(ii,:) = newrates.depos(ii,:) + 1e-4;
    % Store the X value after rate shift
    for jj=1:length(CapStrike)
        Delta_depos(ii,jj) = priceCapViaSpotVol(dates, quotedSwaps, newrates, Strikes, CapStrike(jj), CapMaturity(jj), spotvol)-priceCapViaSpotVol(dates, quotedSwaps, rates, Strikes, CapStrike(jj), CapMaturity(jj), spotvol);
    end
end

% Futures Delta bucket inizialization
Delta_futures = zeros(7,length(CapStrike));
for ii = 1:7
    % Rates struct inizialization
    newrates = rates;
    % Increase each bucket iteratively
    newrates.futures(ii,:) = newrates.futures(ii,:) + 1e-4;
    % Store the X value after rate shift
    for jj=1:length(CapStrike)
        Delta_futures(ii,jj) = priceCapViaSpotVol(dates, quotedSwaps, newrates, Strikes, CapStrike(jj), CapMaturity(jj), spotvol)-priceCapViaSpotVol(dates, quotedSwaps, rates, Strikes, CapStrike(jj), CapMaturity(jj), spotvol);
    end
end

% Swaps Delta bucket inizialization
Delta_swaps = zeros(size(rates.swaps,1)-1, length(CapStrike));
for ii = 2:size(rates.swaps,1)
    % Rates struct inizialization
    newrates = rates;
    % Increase each bucket iteratively
    newrates.swaps(ii,:) = newrates.swaps(ii,:) + 1e-4;
    % Store the X value after rate shift
    for jj=1:length(CapStrike)
        Delta_swaps(ii-1,jj) = priceCapViaSpotVol(dates, quotedSwaps, newrates, Strikes, CapStrike(jj), CapMaturity(jj), spotvol)-priceCapViaSpotVol(dates, quotedSwaps, rates, Strikes, CapStrike(jj), CapMaturity(jj), spotvol);
    end
end
%Aggregating results
DeltaCap = [Delta_depos;Delta_futures;Delta_swaps];
end