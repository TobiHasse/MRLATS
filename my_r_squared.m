function [r2] = my_r_squared(observed, predicted, log_transform,rm_zeros)
% this function should compute R2 for liniar regression
% written by Tobias Hasse tobiack@udel.edu July 25, 2020
if rm_zeros
    predicted(observed == 0) = [];
    observed( observed == 0) = [];
end
if log_transform
    observed = log(observed);
    predicted = log(predicted);
end

avg = mean(observed);
sst = sum((observed-avg).^2);
ssr = sum((observed-predicted).^2);

r2 = 1 - ssr/sst;

end
