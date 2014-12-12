function output = cross_correlation(Phi, r)
% Phi = N x w x M
% Residual = N x T (called r)
% Output = N x T x M

[N, w, M] = size(Phi);
T = size(r, 2);
output = zeros(N, T, M);

padded_r = padarray(r, [0 (w - 1)], 'post');

% padded_r = zeros(N, T + w);
w_half = floor(w/2);

% Pad r
% for j = 1:T
%     padded_r(:, j + w_half) = r(:, j);
% end

% convolve coefficients over all basis elements
for basis_num = 1:M;
%     for t = 1:T;
%         temp = reshape(Phi(:, :, basis_num), N, w) .* padded_r(:, t:t + w - 1);
%         output(:, t, basis_num) = sum(temp, 2);
%     end
    
    flipped_basis = fliplr(Phi(:, :, basis_num));
    padded_flip = padarray(flipped_basis, [0 (T - 1)], 'post');
    for i = 1:N;
        convy = ifft(fft(padded_r(i, :)) .* fft(padded_flip(i, :)));
        output(i, :, basis_num) = convy(1 + floor(w/2):end - floor(w/2));
%         result = conv(flipped_basis(i, :), r(i, :));
%         output1(i, :, basis_num) = result(w_half + 1:end - w_half);
%         norm(output(i, :, basis_num) - output1(i, :, basis_num))
    end
end

end