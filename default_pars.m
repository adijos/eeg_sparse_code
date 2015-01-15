function pars = default_pars(pars,varargin)

% DEFAULT_PARS Retrieves default parameter settings for running the code on images.
%    Fields may be preset in argument PARS. These pre-settings may affect
%    the default values of other parameters. In addition, special parameters
%    may be passed as pairs. For example, to set pars.num_trials = 10000 and
%    set the special parameter 'mode' to 'learn', use
%
%        pars.num_trials = 10000;
%        pars = default_pars(pars, 'save', false);
%
%
%  List of regular parameters
%
%  Data parameters
%   batch_size               Number of random samples per batch
%   num_files                Number of data files to train over
%   T                        size of data samples (with remplacement), 3s
%   N                        number of electrodes
%   M                        number of dictionary elements
%   w                        lenght of basis window, 1s
%   
%   Noise and sparse prior params
%    noise_var
%    beta
%    sigma
%    eta_a
%
%   Learning parameters
%    eta                     Dictionary learning rate
%    VAR_GOAL
%    var_eta
%    alpha
%
% adapted from default_pars.m in sisc by ng

if mod(nargin,2) == 0, error('Optional arguments must come in pairs.'); end;

special_pars = struct;

for curr=2:nargin-1,
  field = varargin{curr-1};
  val = varargin{curr};
  special_pars.(field) = val;
end

if ~isfield(pars, 'batch_size'), pars.batch_size = 1; end;
if ~isfield(pars, 'num_files'), pars.num_files = 20; end;
if ~isfield(pars, 'T'), pars.T = 1200; end;
if ~isfield(pars, 'N'), pars.N = 16; end;
if ~isfield(pars, 'M'), pars.M = 32; end;
if ~isfield(pars, 'w'), pars.w = 401; end;

if ~isfield(pars, 'noise_var'), pars.noise_var = 0.1; end;
if ~isfield(pars, 'beta'), pars.beta = 2; end;
if ~isfield(pars, 'sigma'), pars.sigma = 1; end;
if ~isfield(pars, 'eta_a'), pars.eta_a = 1e-5; end;

if ~isfield(pars, 'eta'), pars.eta = 1; end;
if ~isfield(pars, 'VAR_GOAL'), pars.VAR_GOAL = 1; end;
if ~isfield(pars, 'var_eta'), pars.var_eta = 0.1; end;
if ~isfield(pars, 'alpha'), pars.alpha = 0.02; end;