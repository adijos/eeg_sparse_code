% get default parameters
pars = default_pars(struct);
pars.basis_M = 400; % time window of basis (400hz -> 1s)
pars.basis_N = 1; % time of bases 1 dimensional
pars.num_bases = 20; % 20 basis elements to learn (made 20 so diff from 16)
pars.patch_M = 1000; % time window of patches (4000hz -> 10s)
pars.patch_N = 1; % time patches 1 dimensional
num_pre_data_files = 20; % ONE FOR TESTING!
coef_pars = default_coef_pars(struct);
pars.beta = 130; % 200 default

patch_M = pars.patch_M;
patch_N = pars.patch_N;
num_channels = 16; % 16 "channels" or electrodes
num_patches = 1000; % number of patches
num_test_patches = 1; % number of patches to test on
batch_size = 50;
num_trials = 10;

base_file_pre = './data/Dog_2/Dog_2_preictal_segment_00';
% 
% A = load('bases.mat');
% A = A.A;

%Generate white noise bases, normalize
rand('seed', 100);
A = randn(pars.basis_M,pars.basis_N,num_channels,pars.num_bases);

% normalize over the basis elements!
for m=1:pars.num_bases,
    A(:,:,:,m)=A(:,:,:,m)-mean(mean(mean(A(:,:,:,m)))); % mean 0
    A(:,:,:,m)=A(:,:,:,m)/sqrt(mean(mean(mean(A(:,:,:,m).^2)))); % norm 1
end


% Main loop
% if reshaped, S = patch_M x patch_N x (number of bases)
% in our case, S = 2500 x 1 x 20
% must be 2-D for sparse matrix representation
s_all = sparse(patch_M*patch_N*pars.num_bases,num_patches); 
total_time = 0;
lambda = ones(pars.num_bases,1);

% over loop
% load random raw data file structure
% take data field from large struct
% take data from inner struct

num_batches = floor(num_patches/batch_size);
if num_batches == 0, error('Not enough patches for that batch size.'); end;

for k = 1:num_trials
for i = 1:num_pre_data_files
    sprintf(strcat('Data File #', int2str(i)))
    
    if i < 10
        datatr = load(strcat(base_file_pre,'0',int2str(i)));
    else
        datatr = load(strcat(base_file_pre, int2str(i)));
    end
    
    test_i = num_pre_data_files+i;
    if test_i < 10
        datate = load(strcat(base_file_pre,'0',int2str(test_i)));
    else
        datate = load(strcat(base_file_pre,int2str(test_i)));
    end
    
    data = getfield(datatr, strcat('preictal_segment_', int2str(i)));
    testy = getfield(datate, strcat('preictal_segment_', int2str(test_i)));
    eeg = data.data';
    eeg_t = testy.data';
    
    % reshape so eeg is...
    % size(eeg) = (number of time points) x (1d time) x (16 channels) x (1
    % image)
    % i.e. size(eeg) = T x 1 x 16 x 1 
    eeg = reshape(eeg,[size(eeg,1) 1 size(eeg,2) 1]);
    eeg_t = reshape(eeg_t,[size(eeg,1) 1 size(eeg_t,2) 1]);
    [M, N, num_channels, num_data] = size(eeg);
    
    
    % make eeg data zero mean, unit variance
    eeg = eeg - mean(mean(mean(mean(eeg))));
    eeg = eeg/sqrt(mean(mean(mean(mean(eeg.^2)))));
    
    % make test data zero mean, unit variance
    eeg_t = eeg_t - mean(mean(mean(mean(eeg_t))));
    eeg_t = eeg_t/sqrt(mean(mean(mean(mean(eeg_t.^2)))));
    
    % extract 100 random patches from data
    % loop over range(num_trials=100) and choose random patches accordingly
    % patches should be X_ALL M x N x num_channels x num_patches
    % in our case ..........  T x 1 x 16 x (random large number)
    idxs = randi(size(eeg, 1) - patch_M,100);
    X_ALL_train = zeros(patch_M,patch_N,num_channels,num_patches);
    for nummy=1:num_patches;
        index = idxs(nummy);
        data = eeg(index:index + patch_M - 1, :, :, 1);
        X_ALL_train(:,:,:,nummy) = data;
    end
    
    idxt = randi(size(eeg_t, 1) - patch_M, num_test_patches);
    X_ALL_test = zeros(patch_M,patch_N,num_channels,num_test_patches);
    for nummy=1:num_test_patches;
        index = idxt(nummy);
        data = eeg_t(index:index + patch_M - 1, :, :, 1);
        X_ALL_test(:,:,:,nummy) = data;
    end
    
    %error('sleepeh brah')
    
    for tr=1:pars.num_trials,
        batch = mod(tr-1,num_batches);
        patches = batch_size*batch+(1:batch_size);
        [A,batch_stats,lambda,s_all(:,patches), rec1] = run_batcheeg(X_ALL_train(:,:,:,patches),A,pars,coef_pars,tr,true,lambda,s_all(:,patches));
        
        % Save timing information, compute objective function on test set
        total_time = total_time + batch_stats.coef_time_total + batch_stats.bases_time;
        train_time(i+tr) = total_time/60;
        
        dummy_pars = coef_pars;
        dummy_pars.exact = true;
        dummy_pars.coeff_iter = 10000;
        dummy_pars.num_coords = 100;
        [A,batch_stats, a, aa, rec2] = run_batcheeg(X_ALL_test,A,pars,dummy_pars,0,false);
        fobj_all(i+tr) = batch_stats.fobj_pre_total;
        figure(2); plot(fobj_all); title('Objective function by batch.');
        figure(3); plot(train_time,fobj_all); title('Objective function by time.');
        drawnow;
        
        if pars.display_bases_every ~= 0 && mod(i+tr,pars.display_bases_every)==0,
            figure(1000)
            for channel=1:16;
                subplot(4,4,channel);
                hold all
                for basis=1:pars.num_bases;
                    plot(A(:,1,channel,basis));
                end
            end
            suptitle('all bases on all channels')
            
            figure(1234)
            choice = 7;
            for channel=1:16;
                subplot(4,4,channel);
                plot(A(:,1,channel,choice))
            end
            suptitle(['all channels of basis element ',int2str(choice)]);
            
            
            figure(193487)
            S_vis = reshape(full(s_all(:,choice)),pars.patch_M,1,pars.num_bases);
            for arty=1:pars.num_bases
                subplot(5,4,arty)
                plot(S_vis(:,1,arty))
            end
            suptitle(['activations for basis ', int2str(choice)]);
        end
        
    end
end
end