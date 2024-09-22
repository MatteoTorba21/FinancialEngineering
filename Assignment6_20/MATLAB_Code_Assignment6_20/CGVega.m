function [VegaCGX,VegaCGCaps] = CGVega(VegaBucket, dates,SwapExpiries, CGBuckets, YTM, VegaCap)
%UNTITLED3 
SwDates = findSwapDates(dates.settlement, SwapExpiries, eurCalendar)';
% CG peaks dates
CGPeaks = findSwapDates(dates.settlement, CGBuckets, eurCalendar)';
MaturityDates = findSwapDates(dates.settlement,YTM, eurCalendar)';
% Vector of all quoted dates
% BucketDates = [dates.depos;dates.futures(:,2);SwDates];
% Vector of each bucket start
startBucket = [dates.settlement;CGPeaks(1:end-1)];
% Vector of each bucket end
endBucket = [CGPeaks(2:end);SwDates(find(CGPeaks(end)==SwDates)+1)];
% Upfront price vector for each bucket inizialization
VegaCGX = zeros(size(CGBuckets));
% Swaps matrix inizialization
VegaCGCaps = zeros(length(CGBuckets));
% Bucket cycle
for ii = 1: length(CGPeaks)
    % Check if it is first bucket
    if ii == 1
        % CG weights computation
        weights = interp1([startBucket(ii) CGPeaks(ii) endBucket(ii)], [1,1,0], MaturityDates,'linear');
    else
        % CG weights computation
        weights = interp1([startBucket(ii) CGPeaks(ii) endBucket(ii)], [0,1,0], MaturityDates,'linear');
    end
    weights(isnan(weights)) = 0;
    VegaCGX(ii) = weights'*VegaBucket;
    VegaCGCaps(ii,:) = weights'*VegaCap;
end