function plot_timePhi(Phi)
% assumes Phi is (# of bases x dimensions x time points)
[N, w, M] = size(Phi);
figure(1)
for figgy = 1:M;
    subplot(4, 4, figgy);
    concat_dims = zeros(16, w);
    for i = 1:N;
        concat_dims(i, :) = reshape(Phi(i, :, figgy), 1, w);
    end
    plot(concat_dims)
end
    
end