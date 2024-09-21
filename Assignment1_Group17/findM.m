function M_hat=findM(errors,bps,M)
% This functions gives back the value of M which gives an error
% lower than the required threshold
%
%INPUT:
% errors:  vector of the errors
% bps:     threshold required
% M:       vector of the time steps or Monte-Carlo simulations
%
%OUTPUT:
% M_hat:   number of time steps or Monte-Carlo simulations needed to reach the treshold

M_hat = 0;

for ii=1:length(errors)
    if errors(ii)<=bps
        M_hat = M(ii);
    end
    if M_hat ~= 0
        break
    end
end

if M_hat == 0
    M_hat = M(end);
    error('The threshold is not reached')
end

end %function findM