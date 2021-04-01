function X = getPotency(num,ops)
%getPotency identifies the potency of a number or a vector
%
% [base potency] = getPotency(x,options)
%
% option = 1 : base is rounded down
% option = 2 : base is not rounded down
%
% Examples:
% x = 3.4e7;
% getpotency(x,1) is [3 7]
% getpotency(x,2) is [3.4 7]
%
% y = [0.01 2000];
% getpotency(y,1) is [1 -2]
%                    [2  3]

%% Check for pos. or neg. potency

ResultSize = max(size(num));
num_koeff = zeros(ResultSize,1);
num_round = zeros(ResultSize,1);
num_pot   = zeros(ResultSize,1);

for Iter = 1 : ResultSize
    tmp = num2str(num(Iter),'%10.8e');              % transform number into #.###e+00 format
    
    num_koeff(Iter) = str2double(tmp(1:10));        % define coefficient
    num_round(Iter) = floor(num_koeff(Iter));       % round coefficient down
    num_pot(Iter)   = str2double(tmp(end-2:end));   % get potency from e+00
end

if     ( ops == 1 );  X = [num_round num_pot];
elseif ( ops == 2 );  X = [num_koeff num_pot];
else   error('getPotency must be called like this: X = getPotency(num,ops)');
end

