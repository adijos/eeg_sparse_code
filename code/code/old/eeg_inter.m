% BAYJAM - Dog 2%
% Each data file consists of 10 min of data of 500 Hz EEG from 16
% electrodes on a canine. The number of samples per file is 239766.

num_trials = 10;

num_pre_data_files = 25;
num_inter_data_files = 25;

base_file_pre = './Dog_2/Dog_2_preictal_segment_00';
base_file_inter = './Dog_2/Dog_2_interictal_segment_0';

% Number of dimensions (# of electrodes)
N = 16;

% number of dictionary elements
M = 16;

% Basis Time window
w = 501; % 500 = order of 1 second

% Phi - Dictionary
Phi = randn(N, w, M);
dPhi = zeros(N, w, M);

% normalize dictionary
for i = 1:M
    for t = 1:w
        Phi(:, t, i) = Phi(:, t, i) * diag(1./sqrt(sum(Phi(:, t, i) .* Phi(:, t, i))));
    end
end

% Dictionary Learning Rate
deta = 1e-3;

% Error storage to track reconstruction error across trials of dictionary
% learning
error = [];

for t = 1:num_trials
    sprintf(strcat('Trial #', int2str(t)))
    
    for i = 1:num_pre_data_files
        sprintf(strcat('Data File #', int2str(i)))
        
        if i < 10
            data = load(strcat(base_file_inter,'00',int2str(i)));
        else
            data = load(strcat(base_file_inter, '0', int2str(i)));
        end
        data = getfield(data, strcat('interictal_segment_', int2str(i)));
        full_data = data.data;
        
        sprintf('Data Loaded')
        
        % calculate coefficients for these data
        for repeat = 1:50
            % random sample of 500 * 5 (5 seconds)
            index = randi(size(full_data, 2) - 500 * 5)
            data = full_data(:, index:index + 500 * 5 - 1);
            pre_a = time_sparsify(data, Phi, t);
        
            sprintf('Coefficients Learned')
        
            recon = reconstruct(Phi, pre_a);
            error = [error, sum((data - recon).^2,2)];
            % update bases
        
            e = data - recon;
             dPhi = dict_correlation(pre_a, e, w);
             Phi = Phi + deta * dPhi;
            sprintf('Dictionary Learned');

            % normalize dictionary
            for j = 1:M
                for k = 1:w
                    Phi(:, k, j) = Phi(:, k, j) * diag(1./sqrt(sum(Phi(:, k, j) .* Phi(:, k, j))));
                end
            end

            figure(5)
            for j=1:16;
                subplot(4,4,j)
                plot(Phi(1, :, j));
                xlim([0 500]);
                ylim([-1.5 1.5]);
            end
            suptitle(strcat('Electrode 1 Trial #', int2str(t),'; subsample #', int2str(repeat)));
            figure(6)
            for j=1:16;
                subplot(4,4,j)
                plot(Phi(2, :, j));
                xlim([0 500]);
                ylim([-1.5 1.5]);
            end
            suptitle(strcat('Electrode 2 Trial #', int2str(t),'; subsample #', int2str(repeat)));
            figure(7)
            for j=1:16;
                subplot(4,4,j)
                plot(Phi(3, :, j));
                xlim([0 500]);
                ylim([-1.5 1.5]);
            end
            suptitle(strcat('Electrode 3 Trial #', int2str(t),'; subsample #', int2str(repeat)));
        end
    end
    
end

%pre_a = time_sparsify(data, Phi, 0);
%recon = reconstruct(Phi, pre_a);
%plot_recon(data, reconstruct(Phi, pre_a));

 figure(100)
 for i=1:16;
     subplot(4,4,i)
     plot(error(i,:))
 end
