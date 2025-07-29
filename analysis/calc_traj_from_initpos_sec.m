% 
% Calc particle traj from a section of the pre-calculated initial positions (*.nc).
%   The sub-section is chosen to avoid excessive memory use! 
% 
clear
homedir = getenv('HOME');
workdir = getenv('WORK');
campdir = getenv('CAMP');
addpath(genpath([homedir '/work_Lagr']));
addpath(genpath([homedir '/work_MOM']));
addpath(genpath([homedir '/mytoolbox']));
addpath(genpath([homedir '/MyFuncs']));

%% params config from shell script

tloop = lpSh; % lpS-h: time interval for the integration
isub_region = isubSh; % isubS-h: release ptcls from which part of domain 

% ------------ time interval
% Integrate ptcls for 120 days and then release a new set
%   e.g., D1-131; D121-D251
yr_s = 21;
day_interv = 130;   % length of integration [d]
% Start day [d] - release every 120 days, e.g., D1-131, D121-251
day_s = (tloop-1)*120 + 1;
day_e = day_s + day_interv;
fprintf(1,'Integrate trajs from D%03d to D%03d...\n\n', day_s, day_e);

%% dirs and traj params

% ------- all init pos
init_fnm = [campdir '/lagr_study/ptcl_initpos/ptclset1.nc'];
xin2d_al = ncread(init_fnm,'xin2d'); % [km]
yin2d_al = ncread(init_fnm,'yin2d');
[npx_al, npy_al] = size(xin2d_al);

% ------- grid
grid_dir = [workdir '/MOM6_exp/swm_spunup/'];  
[grid, ~, ~] = read_grid_MOM([grid_dir '']); 

% ------- u v dir
uvh_dir = [campdir '/mom_ptemp/sol_prog/'];
flx_file_preStr = 'prog__';

% ------- parameters for the traj integration
biliORcubic = 'cubic'; % "bilinear" or "cubic" interpolation
ik = 1;

% time step for uv [d]
dt_uv = 6/24;
tuv = day_s:dt_uv:day_e;
ntuv = length(tuv);
% time step for traj [d]
dt_traj = 2/24; 
ttraj = day_s:dt_traj:day_e;
nttraj = length(ttraj);
%  for save
dt_save = 12/24;
tsave = day_s:dt_save:day_e;
ntsave = length(tsave);
% 
it_save = 1:dt_save/dt_traj:nttraj; % id of elements in ttraj that = tsave

% for ODE solver
options = odeset('RelTol',10^(-6),'AbsTol',10^(-6)); 
% convert u/v [m/s] to [km/d]
factor = 3600*24/1e3; 
% coordinate 
[xu,yu] = deal(grid.lonq, grid.lath); 
[xv,yv] = deal(grid.lonh, grid.latq); 

%% read uv
u_t = zeros(grid.niu,grid.nju,ntuv);
v_t = zeros(grid.niv,grid.njv,ntuv);

for it = 1:ntuv
    % times
    [yrstr, dystr, hrstr] = get_timestr(tuv(it), yr_s); 
    % read
    flx_fnm = [uvh_dir flx_file_preStr yrstr '_' dystr '_' hrstr '.nc'];
    u = ncread(flx_fnm,'u');
    v = ncread(flx_fnm,'v');
    % u and v [m/s] --> [km/d]
    [u_t(:,:,it), v_t(:,:,it)] = deal(u(:,:,ik)*factor, v(:,:,ik)*factor);
end
fprintf('Readed uv.\n');

%% ptcls to be integrated in this script

% divide all into sections
[i_s,di,i_e] = deal(1, 90, npx_al); % 100 uses 2h (bilin, 1e-6)
i_pcels_al = cell( 1, ceil((i_e - i_s) / di) );
for ii = 1:size(i_pcels_al,2)
    [is, ie] = deal(i_s + di * (ii-1), i_s + di * ii - 1);
    if ie > npx_al
        ie = npx_al;
    end
    i_pcels_al{ii} = is:ie;
end
clearvars i_s di i_e ii is ie
% ptcls to be integrated
i_pcels = i_pcels_al{isub_region};
j_pcels = 1:npy_al;
[xin2d, yin2d] = deal(xin2d_al(i_pcels,j_pcels), yin2d_al(i_pcels,j_pcels));
[npx, npy] = size(xin2d);

% ----- dir for save
save_dir = [campdir '/lagr_study/exp1/trajs_' biliORcubic '/lp' num2str(tloop,'%02d')...
    '/Z' num2str(ik,'%02d') '/sec' num2str(isub_region,'%02d') '/'];
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
params_fnm = [save_dir 'params.mat'];
save(params_fnm, 'ik','init_fnm','tuv','ttraj','tsave','dt_*','options',...
    'i_pcels_al','j_pcels','xin2d','yin2d','isub_region')

%% int

% --- [km]
fprintf('Integrating trajs...\n');
tic;
[xtr, ytr] = calc_traj_func(xin2d,yin2d,u_t,v_t,tuv,ttraj,xu,yu,xv,yv,...
    options,biliORcubic);
toc;

% --- calc lagrangian velocity [m/s] from trajs
[ul, vl] = calc_Lagr_vel_CG(xtr,ytr,ttraj,1);
% [km/d] --> [m/s]
[ul, vl] = deal(ul/factor, vl/factor);

% --- reduce to tsave
[xtr_save, ytr_save] = deal(xtr(:,:,it_save), ytr(:,:,it_save));
[ul_save, vl_save] = deal(ul(:,:,it_save), vl(:,:,it_save));

%% save trajs (one file for one snapshot)

for it = 1:ntsave
    [yrstr, dystr, hrstr] = get_timestr(tsave(it), yr_s); 
    savefnm = [save_dir '/trajs__' yrstr '_' dystr '_' hrstr '.nc'];
    if exist(savefnm,'file') 
        disp(['NC file already exists, skip : ', savefnm]);
        continue
    end  
    fprintf(1,'\nPtcls pos saving to: %s...\n\n', savefnm);
    
    dim_name = {'Nx','Ny'};
    dim_length = [npx, npy];
    varname = {'xtr','ytr','ul','vl'};
    data = {xtr_save(:,:,it), ytr_save(:,:,it), ul_save(:,:,it), vl_save(:,:,it)};
    dimNum_of_var = {[1,2], [1,2], [1,2], [1,2]};
    %
    global_att  = [ 'Ptcl trajs [km] and vel [m/s] at instant time; ' ...
        'init_fnm=' init_fnm '; params_fnm=' params_fnm ...
        '; day_s=' num2str(day_s) '; day_e=' num2str(day_e)];
    FUN_nc_easywrite_enhanced( savefnm, dim_name, dim_length,...
        varname, dimNum_of_var, data, global_att )
end


