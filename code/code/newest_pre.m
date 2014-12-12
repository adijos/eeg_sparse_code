% sparsenet.m - sparse coding algorithm for video
% 
% Before running you must first define Phi (L x M x szt)

%szt = window size (our w)
%M = # of bases
%L = # of dimensions
%e = reconstruction error
%a = coefficients

% number of trials per movie chunk
batch_size = 10;

base_file_pre = './Dog_2/Dog_2_preictal_segment_00';

num_pre_data_files = 20;

N = 16;
M = 16;
w = 501;
T = 2500;

% margin of bad data
BUFF=4;
topmargin=15;  % for van hateren movies

% margin in time (assumes szt odd)
marg = floor(w/2);
num_valid = T - 2 * marg;

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

% Phi - Dictionary
Phi = randn(N, M, w);

% normalize bases to match desired output variance

gain = gain .* ((a_var/VAR_GOAL).^alpha);
normPhi=sqrt(sum(sum(Phi.*Phi),3));
for i=1:M
    Phi(:,i,:)=gain(i)*Phi(:,i,:)/normPhi(i);
end

while 1

  % choose a movie for this batch
  for i = 1:num_pre_data_files
      sprintf(strcat('Data File #', int2str(i)))
        
      if i < 10
          data = load(strcat(base_file_pre,'0',int2str(i)));
      else
          data = load(strcat(base_file_pre, int2str(i)));
      end
      data = getfield(data, strcat('preictal_segment_', int2str(i)));
      full_data = data.data;
      
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
    
        % collect update stats

        for t = 1:w
          tt = t:T - (w - t);
          dPhi(:, :, t) = dPhi(:, :, t) + e(:, tt) * a';
        end

        a_var = (1 - var_eta) * a_var + var_eta * sum(a .* a,2)/num_valid;
        
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
  % update bases

  fprintf('updating bases\n');

  dPhi = dPhi/(batch_size * num_valid);
  Phi = Phi + eta * dPhi;

  % normalize bases to match desired output variance

  gain = gain .* ((a_var/VAR_GOAL).^alpha);
  normPhi=sqrt(sum(sum(Phi.*Phi),3));
  for i=1:M
    Phi(:,i,:)=gain(i)*Phi(:,i,:)/normPhi(i);
  end

  % display network
  
end
