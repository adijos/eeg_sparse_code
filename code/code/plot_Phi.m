Phi = load('pre_phi.mat');
Phi = Phi.Phi;
[N M w] = size(Phi);

for i = 1:N
    figure(i)
    for j=1:M;
        subplot(4,4,j)
        plot(reshape(Phi(i, j, :), 1, w));
        xlim([0 500]);
    end
    suptitle(strcat('Electrode #', int2str(i)));
end