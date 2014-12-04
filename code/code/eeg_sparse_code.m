% BAYJAM - Dog 2%

num_trials = 100;
num_pre_data_files = 1;
num_inter_data_files = 500;

base_pre_file = './processed_dog_2/pre_';
base_inter_file = './processed_dog_2/inter_';

% # of time samples (each file spans 10 minutes, 600 seconds)
T = 3997;

% number of dictionary elements
M = 16;

% number of dimensions (electrodes)
N = 16;

lambda = 0.1;

% Basis Time window
w = 125; % Bruno said 125, we should try to figure out why (1/4 of 500 Hz)

% Phi - Dictionary
pre_Phi = randn(N, w, M);
dPhi = zeros(N, w, M);

% normalize dictionary
for i = 1:M
    for t = 1:w
        pre_Phi(:, t, i) = pre_Phi(:, t, i) * diag(1./sqrt(sum(pre_Phi(:, t, i) .* pre_Phi(:, t, i))));
    end
end

% Dictionary Learning Rate
eta = 1e-5;

error = [];
for t = 1:num_trials
    sprintf(strcat('Trial #', int2str(t)))
    
    for i = 1:num_pre_data_files
        sprintf(strcat('Data File #', int2str(i)))
        
        %Load data file
        data = load(strcat(base_pre_file, int2str(i), '.mat'));
        data = data.data; % 16 x time sampling amount
        
        % calculate coefficients for these data

        pre_a = time_sparsify(data, pre_Phi, lambda);
        
        recon = reconstruct(pre_Phi, pre_a);
        error = [error, sum((data - recon).^2,2)];
        % update bases
        
        e = data - recon;
        c = cross_correlation(pre_a, e);
        dPhi = mean(c, 2);
        dPhi = repmat(dPhi, 1, w, 1);
        pre_Phi = pre_Phi + eta * dPhi;
        
        % normalize dictionary
        for i = 1:M
            for t = 1:w
                pre_Phi(:, t, i) = pre_Phi(:, t, i) * diag(1./sqrt(sum(pre_Phi(:, t, i) .* pre_Phi(:, t, i))));
            end
        end
        
    end
    
    figure(5)
    for i=1:16;
        subplot(4,4,i)
        plot(pre_Phi(i, :, 1));
    end
end

pre_a = time_sparsify(data, pre_Phi, lambda);
recon = reconstruct(pre_Phi, pre_a);
plot_recon(data, reconstruct(pre_Phi, pre_a));

 figure(100)
 for i=1:16;
     subplot(4,4,i)
     plot(error(i,:))
 end
