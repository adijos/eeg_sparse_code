function a = time_sparsify(I, Phi, lambda)
% online learning
% N = # of dimensions (electrodes)
% w = size of basis time window
% M = # of dictionary elements
% T = the size of the data
% a = N x T x M
% Phi = N x w x M
% I = 16 x T

% fixed numbers
beta = 1000;
sigma = .1;
eta = 5e-4;
num_iter = 100;


[N, w, M] = size(Phi);
T = size(I, 2);
error = zeros(N, num_iter);
a = randn(N, T, M);

for i = 1:num_iter;
    recon = reconstruct(Phi, a);
    e = I - recon;
    error(:, i) = sum(e.^2, 2);
    a_prime = cross_correlation(Phi, e) - lambda * S_prime(a, beta, sigma);
    a = a + eta * a_prime;
end

figure(2)
for i=1:16;
    subplot(4,4,i)
    plot(error(i,:))
end

figure(3)
for j = 1:16;
    subplot(4, 4, j)
    plot(a(j, :, 1));
end

end