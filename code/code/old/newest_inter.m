% Script for learning preictal dictionaries for Dog 2
% Data files are 10 min of 500 Hz sampled data. 239766 samples per file

% Data Variables
batch_size = 10; % # of samples per data file to train on
base_file_inter = './Dog_2/Dog_2_interictal_segment_0'; % data's home
num_inter_data_files = 20; % # of data files to train on

% Variable Sizes
N = 16; % # of electrodes
M = 16; % # of dictionary elements
w = 501; % length of basis windows
T = 2500; % length of data samples

% buffer margin
marg = floor(w/2);
num_valid = T - 2 * marg;

% noise and sparse prior params
noise_var = 0.005;
beta = 2500;
sigma = 0.116;
eta_a = 0.001;

% Phi - Dictionary
Phi = randn(N, M, w);

% learning parameters
eta = 1;
VAR_GOAL = 0.1;
a_var = VAR_GOAL * ones(M,1);
var_eta = .01;
alpha = .02;
gain = sqrt(sum(sum(Phi .* Phi), 3))';

while 1

  % iterate through all the data files
  for i = 1:num_inter_data_files
      sprintf(strcat('Data File #', int2str(i)))
        
      if i < 10
          data = load(strcat(base_file_inter,'00',int2str(i)));
      else
          data = load(strcat(base_file_inter, '0', int2str(i)));
      end
      data = getfield(data, strcat('interictal_segment_', int2str(i)));
      full_data = data.data;
      
      % Subtract out the mean and make unit variance
      full_data = full_data - mean(mean(full_data));
      full_data = full_data /sqrt(mean(mean(full_data.^2)));
      
      dPhi = zeros(N, M, w);
      

      for i = 1:batch_size

        fprintf('trial %d\n',i);

        a=[];
        while isempty(a);

          % Extract random 5 second eeg clip

          index = randi(size(full_data, 2) - 500 * 5);
          data = full_data(:, index:index + 500 * 5 - 1);
          T = size(data, 2);

          % calculate coefficients

          sf = 1;
          while isempty(a) & sf <= 4
              [a, Ihat, e] = newest_sparsify(Phi, gain, data, noise_var, beta, sigma, eta_a/sf);
              sf = sf * 2;
          end

        end
    
        % Update change in bases

        for t = 1:w
          tt = t:T - (w - t);
          dPhi(:, :, t) = dPhi(:, :, t) + e(:, tt) * a';
        end

        a_var = (1 - var_eta) * a_var + var_eta * sum(a .* a,2)/num_valid;
     

      end
  
  % update bases

  fprintf('updating bases\n');

  dPhi = dPhi/(batch_size * num_valid);
  Phi = Phi + eta * dPhi;
  
          figure(5)
        for j=1:16;
            subplot(4,4,j)
            plot(reshape(Phi(1, j, :),1,w));
            xlim([0 500]);
            ylim([-2.5 2.5]);
        end
        suptitle('Electrode 1');
        figure(6)
        for j=1:16;
            subplot(4,4,j)
            plot(reshape(Phi(2, j, :),1,w));
            xlim([0 500]);
            ylim([-2.5 2.5]);
        end
        suptitle('Electrode 2');
        figure(7)
        for j=1:16;
            subplot(4,4,j)
            plot(reshape(Phi(3, j, :),1,w));
            xlim([0 500]);
            ylim([-2.5 2.5]);
        end
        suptitle('Electrode 3');

  % normalize bases to match desired output variance

  gain = gain .* ((a_var/VAR_GOAL).^alpha);
  normPhi=sqrt(sum(sum(Phi.*Phi),3));
  for i=1:M
    Phi(:,i,:)=gain(i)*Phi(:,i,:)/normPhi(i);
  end
  end

  % display network
  
end
