function [recon] = reconstruct(Phi, a)
% N = # of dimensions (electrodes)
% w = size of basis time window
% M = # of dictionary elements
% T = the size of the data
% a = N x T x M
% Phi = N x w x M
[N, w, M] = size(Phi);
T = size(a, 2);
recon = zeros(N, T);

padded_a = zeros(N, T + w, M);
w_half = floor(w/2);

% pad a with w_half zeros in front and half
for i = 1:M
    for j = 1:T
        padded_a(:, j + w_half, i) = a(:, j, i);
    end
end

% convolve coefficients over all basis elements
for basis_num = 1:M;
%     flipped_basis = fliplr(reshape(Phi(:, :, basis_num), N, w));
%     for t = 1:T;
%         temp = flipped_basis .* reshape(padded_a(:, t:t + w - 1, basis_num), N, w);
%         recon(:, t) = recon(:, t) + sum(temp, 2);
%     end
    for i = 1:N;
        result = conv(reshape(Phi(i, :, M), 1, w), reshape(a(i, :, M), 1, T));
        recon(i, :) = recon(i, :) + result(w_half + 1:end - w_half);
    end
end

end