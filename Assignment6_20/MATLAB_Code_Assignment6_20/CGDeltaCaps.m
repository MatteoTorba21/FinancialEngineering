function [DeltaCGX, DeltaCGSwaps, DeltaCGCaps] = CGDeltaCaps(DeltaBucket, dates, SwapExpiries, CGBuckets, BucketDates, DeltaBucketSw, DeltaBucketCap)
% Function to determine the upfront prices after each coarse grained bucket
% shift and each delta for both swaps and caps having as maturity the coarse grained
% peaks
%
% INPUTS:
% DeltaBucket   
% dates:        Struct of dates of the financial instruments quoted in the market
% SwapExpiries: Vector of the quoted swap expiries
% CGBuckets:    Course grained buckets' peaks
% BucketDates   Vector of dates included in buckets
% DeltaBucketSw:      Matrix of Delta Buckets of the Swaps
% DeltaBucketCap:     Matrix of Delta Buckets of the Caps

% OUTPUT:
% DeltaCGX:           Vector of Upfront's values of the contract after each bucket shift
% DeltaCGSwaps:       Matrix of swaps DV01 (bucket peaks on rows, swap maturities on columns)
% DeltaCGCaps:        Matrix of Caps DV01 (bucket peaks on rows, Cap maturities on columns)


% Quoted swaps dates computation
SwDates = findSwapDates(dates.settlement, SwapExpiries, eurCalendar)';
% CG peaks dates
CGPeaks = findSwapDates(dates.settlement, CGBuckets, eurCalendar)';
% Vector of all quoted dates
% BucketDates = [dates.depos;dates.futures(:,2);SwDates];
% Vector of each bucket start
startBucket = [dates.settlement;CGPeaks(1:end-1)];
% Vector of each bucket end
endBucket = [CGPeaks(2:end);SwDates(find(CGPeaks(end)==SwDates)+1)];
% Upfront price vector for each bucket inizialization
DeltaCGX = zeros(size(CGBuckets));
% Swaps matrix inizialization
DeltaCGSwaps = zeros(length(CGBuckets));
% Caps matrix inizialization
DeltaCGCaps = zeros(length(CGBuckets));
% Bucket cycle
for ii = 1: length(CGPeaks)
    % Check if it is first bucket
    if ii == 1
        % CG weights computation
        weights = interp1([startBucket(ii) CGPeaks(ii) endBucket(ii)], [1,1,0], BucketDates,'linear');
    else
        % CG weights computation
        weights = interp1([startBucket(ii) CGPeaks(ii) endBucket(ii)], [0,1,0], BucketDates,'linear');
    end
    weights(isnan(weights)) = 0;
    
    %Compute the three coarse-grained deltas
    DeltaCGX(ii) = weights'*DeltaBucket;
    DeltaCGSwaps(ii,:) = weights'*DeltaBucketSw;
    DeltaCGCaps(ii,:) = weights'*DeltaBucketCap;
end
end

