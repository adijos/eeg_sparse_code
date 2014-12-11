% BAYJAM - Dog 2%

num_trials = 10;

num_pre_data_files = 1;
num_inter_data_files = 500;

base_file_pre = './Dog_2/Dog_2_preictal_segment_00';
base_file_inter = './Dog_2/Dog_2_interictal_segment_0';

% # of time samples (each file spans 10 minutes, 600 seconds)
T = 239766;

% number of dictionary elements
M = 16;

% number of dimensions (electrodes)
N = 16;

% Basis Time window
w = 501; % Bruno said 125, we should try to figure out why (1/4 of 500 Hz)

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
deta = 5e-4;

error = [];
for t = 1:num_trials
    sprintf(strcat('Trial #', int2str(t)))
    
    for i = 1:num_pre_data_files
        sprintf(strcat('Data File #', int2str(i)))
        
        if i < 10
            data = load(strcat(base_file_pre,'0',int2str(i)));
        else
            data = load(strcat(base_file_pre, int2str(i)));
        end
        data = getfield(data, strcat('preictal_segment_', int2str(i)));
        data = data.data;
        
        % calculate coefficients for these data

        pre_a = time_sparsify(data, pre_Phi, t);
        
        recon = reconstruct(pre_Phi, pre_a);
        error = [error, sum((data - recon).^2,2)];
        % update bases
        
        e = data - recon;
        dPhi = dict_correlation(pre_a, e, w);
        pre_Phi = pre_Phi + deta * dPhi;
        
        % normalize dictionary
        for j = 1:M
            for k = 1:w
                pre_Phi(:, k, j) = pre_Phi(:, k, j) * diag(1./sqrt(sum(pre_Phi(:, k, j) .* pre_Phi(:, k, j))));
            end
        end
        
        figure(5)
        for j=1:16;
            subplot(4,4,j)
            plot(pre_Phi(j, :, 1));
        end
        
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
