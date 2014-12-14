function run_sisceeg()

% get default parameters
pars = default_pars(struct);
pars.basis_M = 500;
pars.basis_N = 1;
pars.num_bases = 16;
pars.patch_M = 2500;
pars.patch_N = 1;
num_pre_data_files = 1; %20; % ONE FOR TESTING!
coef_pars = default_coef_pars(struct);

patch_M = pars.patch_M;
patch_N = 1;%pars.patch_N;
num_channels = 16;
num_patches = 1;
num_trials = 100;

base_file_pre = './data/Dog_2/Dog_2_preictal_segment_00';

% Generate white noise bases, normalize
rand('seed', 100);
A = randn(pars.basis_M,pars.basis_N,num_channels,pars.num_bases);
for m=1:pars.num_bases,
    A(:,:,:,m)=A(:,:,:,m)-mean(mean(mean(A(:,:,:,m))));
    A(:,:,:,m)=A(:,:,:,m)/sqrt(mean(mean(mean(A(:,:,:,m).^2))));
end


% Main loop
s_all = sparse(patch_M*patch_N*pars.num_bases,num_patches);
total_time = 0;
lambda = ones(pars.num_bases,1);

% over loop
% load random raw data file structure
% take data field from large struct
% take data from inner struct

for i = 1:num_pre_data_files
    sprintf(strcat('Data File #', int2str(i)))
    
    if i < 10
        data = load(strcat(base_file_pre,'0',int2str(i)));
    else
        data = load(strcat(base_file_pre, int2str(i)));
    end
    data = getfield(data, strcat('preictal_segment_', int2str(i)));
    IMAGES = data.data;
    
    % Make images zero mean, unit variance
    for ind=1:size(IMAGES,3),
        IMAGES(:,:,ind) = IMAGES(:,:,ind)-mean(mean(IMAGES(:,:,ind)));
        IMAGES(:,:,ind) = IMAGES(:,:,ind)/sqrt(mean(mean(IMAGES(:,:,ind).^2)));
    end
    
    % extract 100 random patches from data
    % loop over range(num_trials=100) and choose random patches accordingly
    % patches should be X_ALL M x N x num_channels x num_patches
    % in our case ..........  T x 1 x 16 x (random large number)
    idxs = randi(size(IMAGES, 2) - 500 * 5,100);
    X_ALL = zeros(patch_M,patch_N,num_channels,num_patches);
    for i=1:num_patches;
        index = idxs(i);
        data = IMAGES(:, index:index + 500 * 5 - 1);
        X_ALL(:,:,:,i) = reshape(data,[2500 1 16]);
    end
    %error('sleepeh brah')
    
    for tr=1:num_trials,
        % solve for coefficients and bases in one batch
        patches = i.*1:num_patches;
        [A,batch_stats,lambda,s_all(:,patches)] = run_batch(X_ALL,A,pars,coef_pars,tr,true,lambda,s_all(:,patches));
        
        % Display and save
        if pars.display_bases_every ~= 0 && mod(tr,pars.display_bases_every)==0,
            display_bases(A,1);
            drawnow;
        end

    end
end