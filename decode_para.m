function [ out ] = decode_para( para )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% the outputs are in range: [0.01:0.01:999.99].
precision = 0.01;
scale_order = 5;

para = para(:)'/precision;

col_out = floor(log10(max(para))/scale_order)+1;
out = zeros(length(para),col_out);

for i = 1:col_out
    order = 10^(scale_order*(i-1));
    out(:,i) = mod(floor(para/order), 10^scale_order)*precision;
end

end

