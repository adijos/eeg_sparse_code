% BAYJAM - Dog 2%
% Each data file consists of 10 min of data of 500 Hz EEG from 16
% electrodes on a canine. The number of samples per file is 239766.

batch_size = 10;
num_trials = 5;

num_pre_data_files = 20;
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
Phi = randn(N, M, w);
dPhi = zeros(N, M, w);

% normalize dictionary
for i = 1:M
    for t = 1:w
        Phi(:, i, t) = Phi(:, i, t) * diag(1./sqrt(sum(Phi(:, i, t) .* Phi(:, i, t))));
    end
end

% Dictionary Learning Rate
eta = 1e-6;

% Error storage to track reconstruction error across trials of dictionary
% learning
error = [];

for k = 1:num_trials
    
    for i = 1:num_pre_data_files
        sprintf(strcat('Data File #', int2str(i)))
        
        if i < 10
            data = load(strcat(base_file_inter, '00', int2str(i)));
        else
            data = load(strcat(base_file_inter, '0', int2str(i)));
        end
        data = getfield(data, strcat('interictal_segment_', int2str(i)));
        full_data = data.data;
        
        sprintf('Data Loaded')
        
        % calculate coefficients for these data
        for repeat = 1:batch_size
            sprintf(strcat('Sample #', num2str(repeat)));
            % random sample of 500 * 5 (5 seconds)
            index = randi(size(full_data, 2) - 500 * 5);
            data = full_data(:, index:index + 500 * 5 - 1);
            T = size(data, 2);
            
            [a, Ihat, error] = sparsify_bruno(Phi, data);
        
            sprintf('Coefficients Learned')
            
            % Learn bases
            for t = 1:w
                tt = t:T - (w - t);
                dPhi(:, :, t) = dPhi(:, :, t) + error(:, tt) * a';
            end
        end
        
        Phi = Phi + eta * dPhi;
        sprintf('Dictionary Learned');

        %normalize dictionary
        for i = 1:M
            for t = 1:w
                Phi(:, i, t) = Phi(:, i, t) * diag(1./sqrt(sum(Phi(:, i, t) .* Phi(:, i, t))));
            end
         end

        figure(5)
        for j=1:16;
            subplot(4,4,j)
            plot(reshape(Phi(1, j, :),1,w));
            xlim([0 500]);
            ylim([-1.5 1.5]);
        end
        suptitle('Electrode 1');
        figure(6)
        for j=1:16;
            subplot(4,4,j)
            plot(reshape(Phi(2, j, :),1,w));
            xlim([0 500]);
            ylim([-1.5 1.5]);
        end
        suptitle('Electrode 2');
        figure(7)
        for j=1:16;
            subplot(4,4,j)
            plot(reshape(Phi(3, j, :),1,w));
            xlim([0 500]);
            ylim([-1.5 1.5]);
        end
        suptitle('Electrode 3');
    end   
end