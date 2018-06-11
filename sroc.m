function [auc, TP, FP] = sroc(A,B)
%   ROC analysis to distinguish two conditions.
%
%   Usage: 
%    [auc, TP, FP] = sroc(A,B)
%
%   Input:
%    A: vector for condition 1
%    B: vector for condition 2
%       Preferrably mean(A) > mean(B), then AUC values above 0.5 will 
%       indicate better classification.
%
%   Output:
%    auc: Area under ROC curve (i.e. the ROC index), calculated with
%         numerical trapezoidal integration
%    TP:  True positive rate
%    FP:  False positive rate


% sorting input vectors in ascending order
A = sort(A);
B = sort(B);

for n=1:length(B)-1
    %specifying all levels of discrimination threshold based on amps in condition B
    cuts(n+1) = (B(n)+B(n+1))/2;
end
%ends with slightly above max variance across conditions
cuts=[cuts max([A; B])+(B(2)-B(1))];

%for each level of discrimination threshold, calculate true positive and false positive rates
for c=1:length(cuts)
    TP(c)=length(find(A>=cuts(c)))/length(A);
    FP(c)=length(find(B>=cuts(c)))/length(B);
end

auc =trapz(sort(FP),sort(TP));

end
