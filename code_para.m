%{
test code:

a = 0.01;b=99.99;c=1.01;
para = code_para(a,b,c);
de_para = decode_para(para);

%}

function [ out_para ] = code_para( varargin )
%CODE_PARA Code parameters into one, the output parameter can be decoded
%with decode_para to get the original parameters. 
% CAUTION: round off error may occur when more than 3 parameters are coded.

% to ensure parameters to be correctly decoded, each parameter should be in
% range [0.01:0.01:999.99]. ALERT: input parameter will be rounded off if it
% has higher precision than two decimal point.

%Niki 2015/6/19.


n_input = length(varargin);
if nargin > 4
    error('round off error may occur when more than 3 parameters are coded')
end
out_para = 0;

precision = 0.01;
scale_order = 5;

for in = 1:n_input
    
    in_para = varargin{in};
    % check input:
    if any(~ismember(round(in_para/precision),1:10^scale_order-1)) % multiply 10 to avoid floating point error.
        fprintf('code_para: input %d,',in)
        error('input parameter should be in range: [0.01:0.01:99.99]')
    end
    
    out_para = out_para + in_para*10^((in-1)*scale_order);
    
end

