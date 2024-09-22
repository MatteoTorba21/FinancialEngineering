function spotvol = calibratevol(Strike,Capprices,sigma1,YTM,dates,discounts,deltas,DayCount)
% Function to calibrate a vector of spot volatilities for a given strike
% and a given vector of maturities (in years) by using the prices of the
% Caps.
%
% INPUTS:
% Strike:       Strike value
% Capprices:    Matrix of the Cap prices for the different strikes and
%               volatilities
% sigma1:       Vector of volatilities given a strike
% YTM:          Maturities of the market volatilities
% dates:        Vector of the payment dates of the Caps (every 3 months)
% discounts:    Vector of the discounts at payment dates
% deltas:       Vector of the year fractions between payment dates
% DayCount:     Day count Act365 convenction for yearfrac
%
% OUTPUTS:
% spotvol:      Vector of the spot volatilities for different YTM and the
%               given strike     

% spotvol: all the spot volatilities for each caplet (3m) 
% deltas

% Initialization of the vector of the spots volatilities: for the first
% year the spot volatilities coincide with the flat volatilities
spotvol =[sigma1*ones(3,1); zeros((YTM(end)-1)*4,1)];
% Vector of the forward discounts:
fwd_discounts = discounts(2:end)./discounts(1:end-1);
for ii = 2:length(YTM)
    % Number of caps in the considered period:
    ncap = (YTM(ii)-YTM(ii-1))*4;
    % Difference between the caps prices:
    DC = Capprices(ii)-Capprices(ii-1);
    % Libor rates:
    L = (1./fwd_discounts((YTM(ii-1)*4+1):(ncap+YTM(ii-1)*4)) -1)./deltas((YTM(ii-1)*4+1):(ncap+YTM(ii-1)*4));
    % Inizialization of previous spot vol in the bootstrap procedure
    sigmaS = spotvol(YTM(ii-1)*4-1);
    % Function handle that computes the spot vols up to the next quoted flat
    % volatility. The function linear interpolates the spot vols in
    % function of the last spot vol of the period
    vecSigma = @(sigmaf) sigmaS + cumsum(deltas((YTM(ii-1)*4+1):(ncap+YTM(ii-1)*4)))/sum(deltas((YTM(ii-1)*4+1):(ncap+YTM(ii-1)*4)))*(sigmaf-sigmaS);
    % Year fraction computation for each caplet maturity
    yf = yearfrac(dates(1),dates(YTM(ii-1)*4+1:YTM(ii-1)*4 + ncap),DayCount);
    % d computation in function of the last spot volatility of the bootstrap period    
    d = @(sigmaf) (L-Strike)./(vecSigma(sigmaf)*1e-4.*sqrt(yf));

    % sumcapl =@(sigmaf)  sum(discounts(YTM(ii-1)*4+2:YTM(ii-1)*4+ncap+1).*deltas((YTM(ii-1)*4+1):(ncap+YTM(ii-1)*4)).*((L-Strike).*normcdf(d(sigmaf)) + 1e-4*vecSigma(sigmaf).*sqrt(yearfrac(dates(1),dates(YTM(ii-1)*4+1:YTM(ii-1)*4 + ncap),DayCount)).*normpdf(d(sigmaf))));
    % Function that computes the sum of the caplets of the bootstrap period
    % in function of the last spot volatility
    sumcapl =@(sigmaf)  sum(discounts(YTM(ii-1)*4+2:YTM(ii-1)*4+ncap+1).*deltas((YTM(ii-1)*4+1):(ncap+YTM(ii-1)*4)).*((L-Strike).*normcdf(d(sigmaf)) + 1e-4*vecSigma(sigmaf).*sqrt(yf).*normpdf(d(sigmaf))));
    % caplets = @(sigmaf) discounts(YTM(ii-1)*4+1:YTM(ii-1)*4+ncap).*deltas((YTM(ii-1)*4+1):(ncap+YTM(ii-1)*4)).*((L-Strike).*normcdf(d(sigmaf)) + 1e-4*vecSigma(sigmaf).*sqrt(yearfrac(dates(1),dates(YTM(ii-1)*4+1:YTM(ii-1)*4 + ncap),DayCount)).*normpdf(d(sigmaf)))
    % Difference of sum of quoted caplets and the caplets found plugging
    % the spot vols to be bootstrapped
    diffeqtn = @(sigmaf) abs(sumcapl(sigmaf) - DC);
    % Optiond of the solver
    options = optimoptions('fsolve', 'TolFun', 1e-16, 'Display', 'off');
    % Calibrated spot volatilities via the fsolve function:
    sigCalibrated = fsolve(diffeqtn,90,options);
    % Calibrated spot volatilities via the fsolve function
    spotvol(YTM(ii-1)*4:(YTM(ii-1)*4+ncap-1)) = vecSigma(sigCalibrated);
    % check = sumcapl((sigCalibrated))-DC
end


end