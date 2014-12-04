function [] = plot_recon(original, recon)
    figure(101)

    [N, T] = size(original);
    N_sqrt = ceil(sqrt(N));
    for i = 1:N
        subplot(N_sqrt, N_sqrt, i);
        hold on
        plot(original(i, :), 'b');
        plot(recon(i, :), 'g');
        hold off
    end
end