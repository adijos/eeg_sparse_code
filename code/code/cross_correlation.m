function output = cross_correlation(Phi, r)
% Phi = N x w x M
% Residual = N x T (called r)
% Output = N x T x M

[N, w, M] = size(Phi);
T = size(r, 2);
output = zeros(N, T, M);

padded_r = zeros(N, T + w);
w_half = floor(w/2);

% Pad r
for j = 1:T
    padded_r(:, j + w_half) = r(:, j);
end

% convolve coefficients over all basis elements
for basis_num = 1:M;
%     for t = 1:T;
%         temp = reshape(Phi(:, :, basis_num), N, w) .* padded_r(:, t:t + w - 1);
%         output(:, t, basis_num) = sum(temp, 2);
%     end
    
    flipped_basis = fliplr(reshape(Phi(:, :, basis_num), N, w));
    for i = 1:N;
        result = conv(flipped_basis(i, :), r(i, :));
        output(i, :, basis_num) = result(w_half + 1:end - w_half);
    end
end

end