function output = dict_correlation(a, error, w)
% w = window size
% error = N x T (called r)
% a = N x T x M

[N, T, M] = size(a);
output = zeros(N, w, M);

padded_a = padarray(a, [0 floor(w/2) 0]);

% convolve coefficients over all basis elements
for basis_num = 1:M;
    for t = 1:w
        for j = 1:N
            output(j, t, basis_num) = (padded_a(j, t:T + t - 1, basis_num) * error(j, :)')/T;
        end
    end
end

end