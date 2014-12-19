function [recon] = reconstruct(Phi, a)
% N = # of dimensions (electrodes)
% w = size of basis time window
% M = # of dictionary elements
% T = the size of the data
% a = N x T x M
% Phi = N x w x M
[N, M, w] = size(Phi);
T = size(a, 2);
recon = zeros(N, T);
recon1 = zeros(N, T);

padded_a = padarray(a, [0 (w - 1)], 'post');
padded_phi = padarray(Phi, [0 0 (T-1)], 'post');
size(padded_a)

% padded_a = zeros(N, T + w - 1, M);
% padded_phi = zeros(N, T + w - 1, M);


% pad a with w_half zeros in front and half
% padded_a(:, 1 + floor(w/2):end - floor(w/2), :) = a;
% padded_phi(:, 1 + floor(T/2): end - floor(T/2), :) = Phi;

% convolve coefficients over all basis elements
for basis_num = 1:M;
%     flipped_basis = fliplr(reshape(Phi(:, :, basis_num), N, w));
%     for t = 1:T;
%         temp = flipped_basis .* reshape(padded_a(:, t:t + w - 1, basis_num), N, w);
%         recon(:, t) = recon(:, t) + sum(temp, 2);
%     end
    
    for i = 1:N;
        convy = ifft(fft(reshape(padded_a(i, M, :), 1, T + w - 1)) .* fft(reshape(padded_phi(i, M, :), 1, T + w - 1)));
        recon(i, :) = recon(i, :) + convy(1 + floor(w/2):end - floor(w/2));
%         recon(i, :) = recon(i, :) + conv(a(i, :, M), Phi(i, :, M), 'same');
%         norm(recon(i, :) - recon1(i, :))
    end
end

end