function a = time_sparsify(I, Phi, trial)
% Function to Learn the Sparse Coefficients
% N = # of dimensions (electrodes)
% w = size of basis time window
% M = # of dictionary elements
% T = the size of the input data (I)
% a = N x T x M
% Phi = N x w x M
% I = N x T
% beta is
% sigma is
% eta is the learning rate for the coefficients
% lambda is the relative contribution of the sparse penalty

% fixed numbers
beta = 100;
sigma = .1;
eta = 2e-4;
num_iter = 100;

% Get relative sizes of things
[N, w, M] = size(Phi);
T = size(I, 2);

% Reconstruction error array
error = zeros(N, num_iter);

% Initalize a randomly
a = randn(N, T, M);

for i = 1:num_iter;
    i
    recon = reconstruct(Phi, a);
    e = I - recon;
    error(:, i) = sum(e.^2, 2);
    a_prime = cross_correlation(Phi, e) - S_prime(a, beta, sigma);
    sprintf('Correlation Completed');
    a = a + eta * a_prime;
end

figure(200*sigma)
for i=1:16;
    subplot(4,4,i)
    plot(error(i,:))
end
suptitle(strcat(num2str(sigma),'Coefficient Reconstruction Error Trial #', int2str(trial)));

figure(300*sigma)
for j = 1:16;
    subplot(4, 4, j)
    plot(a(j, :, 1));
end
suptitle(strcat(num2str(sigma),'Basis 1 Sparse Coefficients Trial #', int2str(trial)));

end