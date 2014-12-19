pars = default_pars(struct);
pars.basis_M = 40; % time window of basis (500hz -> 1s)
pars.basis_N = 1; % time of bases 1 dimensional
pars.num_bases = 20; % 20 basis elements to learn (made 20 so diff from 16)
pars.patch_M = 400; % time window of patches (500hz -> 5s)
pars.patch_N = 1; % time patches 1 dimensional
num_pre_data_files = 20; % ONE FOR TESTING!
coef_pars = default_coef_pars(struct);
pars.beta = 80; % 200 default

patch_M = pars.patch_M;
patch_N = pars.patch_N;
num_channels = 16; % 16 "channels" or electrodes
num_patches = 1000; % number of patches
num_test_patches = 1; % number of patches to test on
batch_size = 50;
num_trials = 10;

s_all = sparse(patch_M*patch_N*pars.num_bases,num_patches); 
total_time = 0;
lambda = ones(pars.num_bases,1);

AtA = get_AtA(A);

num_batches = floor(num_patches/batch_size);

[s,coef_stats, rec] = get_responses(X_ALL_train(:, :, :, 1),A,pars.beta,coef_pars,1,s_all,AtA);