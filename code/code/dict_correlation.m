function output = dict_correlation(a, error, w)
% w = window size
% error = N x T (called r)
% a = N x T x M

[N, T, M] = size(a);
output = zeros(N, w, M);

padded_a = zeros(N, T + w, M);
w_half = floor(w/2);

% Pad a
for k = 1:M
   for j = 1:T
        padded_a(:, j + w_half, k) = a(:, j, k);
   end
end

% convolve coefficients over all basis elements
for basis_num = 1:M;
    for t = 1:w
        for j = 1:N
            %size(output(j, t, basis_num))
            %size(padded_a(j, t:T + t - 1, basis_num))
            %size(error(j, :)')
            output(j, t, basis_num) = padded_a(j, t:T + t - 1, basis_num) * error(j, :)';
        end
    end
end

end