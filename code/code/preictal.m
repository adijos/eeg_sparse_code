% Script for learning preictal dictionaries for Dog 2
% Data files are 10 min of 400 Hz sampled data. 239766 samples per file

% Data Variables
base_file_pre = './Dog_2/Dog_2_preictal_segment_00'; % data's home
batch_size = 100; % # of random samples per batch 
num_files = 20; % # of data files to train over
T = 1200; % size of data samples (with replacement), 3 seconds

% Dimensions
N = 16; % # of electrodes
M = 16; % # of dictionary elements
w = 401; % length of basis windows (1 second)

% buffer margin to pad
marg = floor(w/2);
num_valid = T - 2 * marg;

% noise and sparse prior params
noise_var = 0.005;
beta = 1000;
sigma = .01;
eta_a = 1e-9;

% Phi - Dictionary initialized randomly
Phi = randn(N, M, w);

% learning parameters
eta = 50; % dictionary learning rate
VAR_GOAL = 0.1; % ??
a_var = VAR_GOAL * ones(M,1); % ??
var_eta = 1e-20; % coefficients learning rate
alpha = .02; % ??
gain = sqrt(sum(sum(Phi .* Phi), 3))'; % ??

% Objective Function Record
E_record = [];

while 1
    
  % iterate through all the data files
  for i = 1:num_files
      sprintf(strcat('Data File #', int2str(i)))
        
      if i < 10
          data = load(strcat(base_file_pre,'0',int2str(i)));
      else
          data = load(strcat(base_file_pre, int2str(i)));
      end
      
      data = getfield(data, strcat('preictal_segment_', int2str(i)));
      full_data = data.data;
      
      % Subtract out the mean and make unit variance
      full_data = full_data - mean(mean(full_data));
      full_data = full_data /sqrt(mean(mean(full_data.^2)));
      
      dPhi = zeros(N, M, w);
      
      for j = 1:batch_size

          fprintf('trial %d\n',i);

          % Extract random 3 second eeg clip
          index = randi(size(full_data, 2) - T);
          data = full_data(:, index:index + T - 1);

          a=[];
          while isempty(a);

              % calculate coefficients
              sf = 1;
              while isempty(a) & sf <= 4
                  [a, Ihat, e] = sparsify(Phi, gain, data, noise_var, beta, sigma, eta_a/sf);
                  sf = sf * 2;
              end

          end
          
          % Update change in bases proportional to correlation between
          % reconstruction error and inferred coefficients
          for t = 1:w
              tt = t:T - (w - t);
              dPhi(:, :, t) = dPhi(:, :, t) + e(:, tt) * a';
          end

          % ??
          a_var = (1 - var_eta) * a_var + var_eta * sum(a .* a,2)/num_valid;
          
      end
      
      % update bases
      dPhi = dPhi/(batch_size * num_valid);
      Phi_old = Phi;
      Phi = Phi + eta * dPhi;
      
      % normalize bases to match desired output variance ??
      gain = gain .* ((a_var/VAR_GOAL).^alpha);
      normPhi=sqrt(sum(sum(Phi.*Phi),3));
      for i=1:M
          Phi(:,i,:)=gain(i)*Phi(:,i,:)/normPhi(i);
      end
      
      % Compute objective function
      E = (0.5 * 1/noise_var * sum(e(:).^2) + beta * sum(S(a(:)/sigma)))/num_valid;
      E_record = [E_record E];
      
      % Plot the objective function across batches
      figure(1);
      plot(E_record);
      title('Objective Function');
      xlabel('Number of Batches');
      
      figure(12);
      for channel=1:16;
        subplot(4,4,channel);
        imagesc(reshape(Phi(channel, :, :), M, w))
        az = 0;
        el = 90;
        view(az, el);
        colormap('gray');
       end
      % 
      
      
%       figure(2);
%       for j = 1:16
%           subplot(4, 4, j);
%           plot(reshape(Phi(1, j, :) - Phi_old(1, j, :), 1, w));
%           xlim([0 500])
%       end
%   
%       figure(2);
%       for j=1:16;
%           subplot(4, 4, j)
%           plot(reshape(Phi(1, j, :), 1, w));
%           xlim([0 w]);
%           ylim([-2.5 2.5]);
%       end
%       suptitle('Electrode 1');
%       
%       figure(6)
%       for j=1:16;
%           subplot(4, 4, j)
%           plot(reshape(Phi(2, j, :),1,w));
%           xlim([0 w]);
%           ylim([-2.5 2.5]);
%       end
%       suptitle('Electrode 2');
%       
%       figure(7)
%       for j=1:16;
%           subplot(4, 4, j)
%           plot(reshape(Phi(3, j, :), 1, w));
%           xlim([0 w]);
%           ylim([-2.5 2.5]);
%       end
%       suptitle('Electrode 3');
  end
end
          