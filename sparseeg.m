function [Phi] = sparseeg(datatype, pars)

% Script for learning preictal dictionaries for Dog 2
% Data files are 10 min of 400 Hz sampled data. 239766 samples per file
% datatype == 0 => preictal
% datatype == 1 => interictal

pars
% Data Variables
if datatype == 0,
    base_file = './data/Dog_2/Dog_2_preictal_segment_00';
else
    base_file = './data/Dog_2/Dog_2_interictal_segment_0';
end
batch_size = pars.batch_size; % # of random samples per batch 
num_files = pars.num_files; % # of data files to train over
T = pars.T; % size of data samples (with replacement), 3 seconds

% visualize a fixed set of data after some number of iterations?
visualize_fixed = 10; %pars.visualize_fixed; % iterations till visualize, default 0
if visualize_fixed > 0,
    fixed_data = load(strcat(base_file,'13'));
    fixed_data = getfield(fixed_data, strcat('preictal_segment_13'));
    full_data = fixed_data.data;
    
    % Subtract out the mean and make unit variance
    full_data = full_data - mean(mean(full_data));
    full_data = full_data /sqrt(mean(mean(full_data.^2)));
    idx = randi(size(full_data, 2) - T);
    fixed_data = full_data(:, idx:idx + T - 1);
end


% Dimensions
N = pars.N; % # of electrodes
M = pars.M; % # of dictionary elements
w = pars.w; % length of basis windows (1 second)

% buffer margin to pad
marg = floor(w/2);
num_valid = T - 2 * marg;

% noise and sparse prior params
noise_var = pars.noise_var; % see lines 9 and 59 in sparsify
% beta determines the sparse penalty. larger beta == greater sparse penalty
beta = pars.beta; 
% sigma determines the threshold at which coefficient values
% either get pushed down to 0 or boosted up. The derivative pushes
% coefficients greater if above this value and down towards zero if below
sigma = pars.sigma;
eta_a = pars.eta_a; % coefficient learning rate

% Phi - Dictionary initialized randomly
Phi = randn(N, M, w);

% learning parameters
eta = pars.eta; % dictionary learning rate
VAR_GOAL = pars.VAR_GOAL; % see
a_var = VAR_GOAL * ones(M,1); % see line 133
var_eta = pars.var_eta; % see line 133
alpha = pars.alpha; % ??
gain = sqrt(sum(sum(Phi .* Phi), 3))'; % ??

% Objective Function Record
E_record = [];
iter = 0; % count training iteration for recording
% Variance Difference Record (a_var - VAR_GOAL)
var_diff = [];

while 1
    
  % iterate through all the data files
  for i = 1:num_files
      sprintf(strcat('Data File #', int2str(i)))
      
      if datatype == 0,
              if i < 10
                  data = load(strcat(base_file,'0',int2str(i)));
              else
                  data = load(strcat(base_file, int2str(i)));
              end

              data = getfield(data, strcat('preictal_segment_', int2str(i)));
              full_data = data.data;
      elseif datatype == 1,
              if i < 10
                  data = load(strcat(base_file,'00',int2str(i)));
              elseif i > 9 && i < 100,
                  data = load(strcat(base_file,'0',int2str(i)));
              else,
                  data = load(strcat(base_file,int2str(i)));
              end
              data = getfield(data, strcat('interictal_segment_',int2str(i)));
              full_data = data.data;
      end
      
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

          % Update the variance of the coefficients based on the variance
          % of the last coefficient sample learned. Var_eta determines to 
          % what degree the coefficient variance is updated according to
          % the new sample. If var_eta is low, the variance will change only
          % slightly on each sample and vice versa.
          a_var = (1 - var_eta) * a_var + var_eta * sum(a .* a,2)/num_valid;
          
      end
      
      % update bases
      dPhi = dPhi/(batch_size * num_valid);
      Phi = Phi + eta * dPhi;
      
      % The gain determines how the dictionaries should be adjusted.
      % VAR_GOAL is the desired variance of the coefficients. In order to
      % push the coefficient variance a_var towards VAR_GOAL, the
      % dictionaries get pushed accordingly in order to compensate. This is
      % what the gain is serving. If a_var > VAR_GOAL, then the dictionary
      % values will increase by the gain so that the coefficient values
      % become smaller and thus reduce their variance. Vice versa if a_var
      % < VAR_GOAL. Alpha is a tiny value use to modulate the ratio. It is
      % unclear why it was selected as it is. Plot this!
      gain = gain .* ((a_var/VAR_GOAL).^alpha);
      normPhi = sqrt(sum(sum(Phi.*Phi),3));
      for i=1:M
          Phi(:,i,:) = gain(i)*Phi(:,i,:)/normPhi(i);
      end
      
      % Plot the difference in coefficient and goal variance
      var_diff = [var_diff (a_var - VAR_GOAL)];
      size(var_diff)
      figure(109);
      hold all;
      for l = 1:M
          plot(var_diff(l, :));
      end
      title('Coefficient Variance - Goal Variance');
      
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
       
      if visualize_fixed > 0 && mod(iter,visualize_fixed) == 0 && iter > 0,
          [a, Ihat, e] = sparsify(Phi, gain, data, noise_var, beta, sigma, eta_a/sf);
          
          figure(323);
          for channel=1:16;
              subplot(4,4,channel);
              plot(a(channel,:))
          end
          suptitle('Sparse Coefficients of Fixed Data')
          
          figure(342);
          for channel=1:16;
              subplot(4,4,channel);
              hold on
              plot(fixed_data(channel,:),'k')
              plot(Ihat(channel,:), 'b')
          end
          suptitle('Reconstructions of Fixed Data')
              
      end
      
      % count iteration of training
      iter = iter + 1
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
          