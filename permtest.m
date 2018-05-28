function [h, p, diff, diff_null] = permtest(a,varargin)
%   Non-parametric permutation test to compare means using
%   within-participant shuffling. Use as a paired-sample permutation test. 
%
%   H0: the two populations have equal means
%   HA: the two populations have unequal means
%
%   Usage:
%     [h p diff diff_null] = permtest(a)                          (compare to zero)
%     [h p diff diff_null] = permtest(a,b)                        (compare to a particular value, or compare the mean of two distributions)
%     [h p diff diff_null] = permtest(a,b,npermutes,pthresh,tail)
%
%   Input:
%     a:          distribution 1
%     b:          (optional) distribution 2, or a single value. default: 0
%     npermutes:  (optional) the number of permutations. default: 1000
%     pthresh:    (optional) the p-threshold for significance. default: 0.05
%     tail:       (optional) test one-tailed (specify 'right' or 'left') or two-tailed (specify 'both'). default: two-tailed
%                   if there is no a priori expectation, then use 'tail' = 'both'
%                   if the a priori expectation is that a > b, then use 'tail' = 'right'
%                   if the a priori expectation is that a < b, then use 'tail' = 'left'
%
%   Output:
%     h:          significance (1 or 0)
%     pval:       p-value of permutation test. Discard H0 if pval is small.
%     diff:       mean(a)-mean(b)
%     diff_null:  the permuted null distribution
%
% RL van den Brink, 2017
% r.l.van.den.brink@fsw.leidenuniv.nl

%% check the input

if nargin == 0
    error('not enough input arguments')
elseif nargin == 1
    b = zeros(size(a));
    npermutes = 1000;
    pthresh   = 0.05;
    tail = 'both';    
elseif nargin == 2
    npermutes = 1000;
    pthresh   = 0.05;
    tail = 'both';
    b = varargin{1};
elseif nargin == 3
    pthresh   = 0.05;
    b = varargin{1};
    npermutes = varargin{2};
    tail = 'both';
elseif nargin == 4
    b = varargin{1};
    npermutes = varargin{2};
    pthresh = varargin{3};
    tail = 'both';
elseif nargin > 5
    error('too many input arguments')
else
    b = varargin{1};
    npermutes = varargin{2};
    pthresh = varargin{3};
    tail = varargin{4};
end

if isempty(b); b                 = zeros(size(a)); end
if isempty(npermutes); npermutes = 1000;   end
if isempty(pthresh); pthresh     = 0.05;   end
if isempty(tail); tail           = 'both'; end

if (length(a) ~= length(b)) && length(b) == 1
    b = zeros(size(a)) + b;
elseif (length(a) ~= length(b)) && length(b) ~= 1
    error('The data in a paired sample test must be the same size.')
end

%make sure the input distributions are row vectors
if size(a,2) > size(a,1); a = a'; end
if size(b,2) > size(b,1); b = b'; end

%% test

%compute difference in mean
diff = mean(a)-mean(b);

%compute permuted null distribution of mean differences
diff_null = zeros(npermutes,1);
for permi = 1:npermutes
    bnull = [a b];
    idx   = rand(size(a)) < 0.5;
    idx   = logical([idx 1-idx]);
    anull = bnull(idx);
    bnull(idx) = [];
    
    diff_null(permi) = mean(anull)-mean(bnull);
end

%calculate p-value
if strcmpi(tail,'both')    
    p = 1-sum(abs(diff)>abs(diff_null))/npermutes;
elseif strcmpi(tail,'right')
    p = 1-sum(diff>diff_null)/npermutes;
elseif strcmpi(tail,'left')
    p = 1-sum(diff<diff_null)/npermutes;
end

%compare p to the alpha level
h = p < pthresh;

end