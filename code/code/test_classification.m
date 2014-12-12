
load pre_phi.mat
load inter_phi.mat

num_pre_data_files = 20;
num_inter_data_files = 500;

base_file_pre = './Dog_2/Dog_2_preictal_segment_00';
base_file_inter = './Dog_2/Dog_2_interictal_segment_0';

% noise and sparse prior params
noise_var = 0.005;
beta = 2.5;
sigma = 0.316;
eta_a = 0.001;

% learning parameters
eta = 1.0;
VAR_GOAL = 0.1;
a_var = VAR_GOAL * ones(M,1);
var_eta = .01;
alpha = .02;
gain = sqrt(sum(sum(Phi .* Phi), 3))';

gain = gain .* ((a_var/VAR_GOAL).^alpha);

pre_correct = 0;
for i = 1:num_pre_data_files
        sprintf(strcat('Data File #', int2str(i)))
        
        if i < 10
            data = load(strcat(base_file_pre,'0',int2str(i)));
        else
            data = load(strcat(base_file_pre, int2str(i)));
        end
        data = getfield(data, strcat('preictal_segment_', int2str(i)));
        full_data = data.data;
        
        a=[];
        while isempty(a);

          % Extract random 5 second eeg clip

          index = randi(size(full_data, 2) - 500 * 5);
          data = full_data(:, index:index + 500 * 5 - 1);
          T = size(data, 2);

          % calculate coefficients

          sf = 1;
          while isempty(a) & sf <= 4
              [a, Ihat, e_pre] = newest_sparsify(pre_Phi, gain, full_data, noise_var, beta, sigma, eta_a/sf);
              sf = sf * 2;
          end
          
          sf = 1;
          while isempty(a) & sf <= 4
              [a, Ihat, e_inter] = newest_sparsify(inter_Phi, gain, full_data, noise_var, beta, sigma, eta_a/sf);
              sf = sf * 2;
          end
          
          if e_pre < e_inter
              pre_correct = pre_correct + 1
          end
          

        end
end
fprintf(strcat('Preictal Accuracy Rate:', num2str(pre_correct/num_pre_data_files)));

inter_correct = 0;
for i = 1:num_inter_data_files
        sprintf(strcat('Data File #', int2str(i)))
        
        if i < 10
            data = load(strcat(base_file_inter,'00',int2str(i)));
        else
            data = load(strcat(base_file_inter, '0', int2str(i)));
        end
        data = getfield(data, strcat('interictal_segment_', int2str(i)));
        full_data = data.data;
        
        a=[];
        while isempty(a);

          % Extract random 5 second eeg clip

          index = randi(size(full_data, 2) - 500 * 5);
          data = full_data(:, index:index + 500 * 5 - 1);
          T = size(data, 2);

          % calculate coefficients

          sf = 1;
          while isempty(a) & sf <= 4
              [a, Ihat, e_pre] = newest_sparsify(pre_Phi, gain, full_data, noise_var, beta, sigma, eta_a/sf);
              sf = sf * 2;
          end
          
          sf = 1;
          while isempty(a) & sf <= 4
              [a, Ihat, e_inter] = newest_sparsify(inter_Phi, gain, full_data, noise_var, beta, sigma, eta_a/sf);
              sf = sf * 2;
          end
          
          if e_pre > e_inter
              inter_correct = inter_correct + 1
          end
          

        end
end
fprintf(strcat('Interictal Accuracy Rate:', num2str(inter_correct/num_inter_data_files)));