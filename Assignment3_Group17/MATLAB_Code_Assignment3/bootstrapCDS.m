function [datesCDS,survProbs,intensities] = bootstrapCDS(datesDF,discounts,datesCDS,spreadsCDS,flag,recovery)
% Bootstrap for a Credit Default Swap
%
% INPUTS:
% datesDF:      Dates obtained from the bootstrap
% discounts:    Discounts obtained from the bootstrap
% datesCDS:     Vector of dates of the Credit Default Swap
% spreadsCDS:   Vector of spreads of the Credit Default Swap
% flag:         If flag is equal to:
%               1 -> approximation neglecting the “accrual” term
%               2 -> exact computation considering the “accrual” term
%               3 -> Jarrow-Turnbull approximation
% recovery:     Recovery value
%
% OUTPUTS:
% datesCDS:     Complete vector of dates of the Credit Default Swap
% survProbs:    Survival probabilities
% intensities:  Intensities
%

% Yearfrac convenction: 30/360 European
SwapDayCount = 6;
% Yearfrac convenction: ACT/365
Act365 = 3;

%survival probability initilization
survProbs = zeros(length(datesCDS)+1,1);

%survival probability in 0
survProbs(1) = 1; 

% intensity initialization
intensities = zeros(length(datesCDS),1);

% Interpolation from bootstapped discounts via rates:
DFCDS = InterpDFviaRates(datesDF,discounts,datesCDS);
% Year fractions to compute the survival probabilities:
deltas = yearfrac([datesDF(1); datesCDS(1:end-1)],datesCDS,SwapDayCount);

switch(flag)

    case 1 % approximation neglecting the accrual
        % Survival probability at 1 year:
        survProbs(2) = (1-recovery)*DFCDS(1)/(spreadsCDS(1)*deltas(1)*DFCDS(1)+(1-recovery)*DFCDS(1));
        % Survival probabilities for years from 2 to length(datesCDS)
        for ii = 2:length(datesCDS)
            survProbs(ii+1) = (-spreadsCDS(ii) * (deltas(1:ii-1).*DFCDS(1:ii-1))' * survProbs(2:ii) ...
                            + (1-recovery) * DFCDS(1:ii-1)' * (survProbs(1:ii-1)-survProbs(2:ii)) ...
                            + (1-recovery) * DFCDS(ii) * survProbs(ii))...
                            /(spreadsCDS(ii)*deltas(ii)*DFCDS(ii) + (1-recovery) * DFCDS(ii));
        end
        intensities = - log(survProbs(2:end)./survProbs(1:end-1))./...
                        yearfrac([datesDF(1); datesCDS(1:end-1)],datesCDS,Act365);
                
    case 2 % exact
        % Survival probability at 1 year:
        survProbs(2) = ((1-recovery)*DFCDS(1)-spreadsCDS(1)*deltas(1)/2*DFCDS(1))/(spreadsCDS(1)*deltas(1)/2*DFCDS(1)+(1-recovery)*DFCDS(1));
        % Survival probabilities for years from 2 to length(datesCDS)
        for ii = 2:length(datesCDS)
            survProbs(ii+1) = (-spreadsCDS(ii) * (survProbs(2:ii)+survProbs(1:ii-1))' * (deltas(1:ii-1).*DFCDS(1:ii-1))/2 ... 
                            + (1-recovery) * DFCDS(1:ii-1)' * (survProbs(1:ii-1)-survProbs(2:ii))...
                            - spreadsCDS(ii) * survProbs(ii) * deltas(ii) * DFCDS(ii)/2 ...
                            + (1-recovery) * DFCDS(ii) * survProbs(ii))...
                            /(spreadsCDS(ii)*deltas(ii)*DFCDS(ii)/2 + (1-recovery)*DFCDS(ii));
        end    
        intensities = - log(survProbs(2:end)./survProbs(1:end-1))./...
                        yearfrac([datesDF(1); datesCDS(1:end-1)],datesCDS,Act365);
                    
    case 3 % Jarrow-Turnbull
        intensities = spreadsCDS(end)/(1 - recovery)*ones(length(datesCDS),1);
        survProbs(2:end) = exp(-intensities.*yearfrac(datesDF(1)*ones(length(intensities),1), datesCDS, Act365));
    otherwise

end

% Remove the probability of default at the settlement date:
survProbs = survProbs(2:end);

end % function bootstrapCDS