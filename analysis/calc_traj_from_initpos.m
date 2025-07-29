function calc_traj_from_initpos(lp)

% Calc particle traj from the pre-calculated initial positions (*.nc).
% 
clear
homedir = getenv('HOME');
workdir = getenv('WORK');
campdir = getenv('CAMP');
addpath(genpath([homedir '/work_MOM']));
addpath(genpath([homedir '/mytoolbox']));
addpath(genpath([homedir '/MyFuncs']));

%% dirs and traj params

lp = 0;

% ------- grid
grid_dir = [workdir '/MOM6_exp/swm_spunup/'];  
[grid, ~, ~] = read_grid_MOM([grid_dir '']); 

% ------- init pos
init_fnm = [campdir '/lagr_study/ptcl_initpos/ptclset1.nc'];
xin2d = ncread(init_fnm,'xin2d'); % [km]
yin2d = ncread(init_fnm,'yin2d');

% ------- u v dir
uvh_dir = [campdir '/mom_ptemp/sol_prog/'];
flx_file_preStr = 'prog__';

% ------- parameters for the traj integration
biliORcubic = 'bilinear'; % "bilinear" or "cubic" interpolation
ik = 1;
[yr_s, day_s, day_e] = deal(21, 1, 5);

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

it_save = 1:dt_save/dt_traj:nttraj; % id of elements in ttraj that = tsave

% for ODE solver
options = odeset('RelTol',10^(-6),'AbsTol',10^(-6)); 
% convert u/v [m/s] to [km/d]
factor = 3600*24/1e3; 
% coordinate 
[xu,yu] = deal(grid.lonq, grid.lath); 
[xv,yv] = deal(grid.lonh, grid.latq); 

% ----- dir for save
save_dir = [campdir '/lagr_study/trajs/run1_' biliORcubic '/Z' ...
    num2str(ik,'%02d') '/'];
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
params_fnm = [save_dir 'params.mat'];
save(params_fnm, 'ik','init_fnm','tuv','ttraj','tsave','dt_*','options')

%% read uv
u_t = zeros(grid.niu,grid.nju,ntuv);
v_t = zeros(grid.niv,grid.njv,ntuv);

fprintf('Reading uv...\n');
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

%% integration - loop over every ptcl

% # of ptcls
[npx, npy] = size(xin2d);
nptcls = npx*npy;
[xtr,ytr] = deal(NaN * zeros(npx, npy, ntsave));
timelen_valid = zeros(npx, npy);
% 
fprintf('Calc trajs of %d ptcls using %s interp...\n',nptcls,biliORcubic);

% loop over each ptcl
tic;

M = 10;
% for ip = 1:npx
parfor ip = 1:npx
    fprintf(1,'\nCalc ptcl traj at x = %d out of %d... \n',ip,npx);

    for jp = 1:npy
        
        % init pos of this ptcl [km] x-y (lon-lat)
        z0 = [xin2d(ip,jp); yin2d(ip,jp)];
        zz = []; 

        % intgrate
        if strcmp(biliORcubic,'bilinear') 
            [~,zz] = ode45(@(t,zz) ...
                HamEqSolver_BiLin_CGrid(t,zz,u_t,v_t,xu,yu,xv,yv,tuv),...
                ttraj, z0(:), options);
        elseif strcmp(biliORcubic,'cubic') 
            [~,zz] = ode45(@(t,zz) ...
                HamEqSolver_Cubic_CGrid(t,zz,u_t,v_t,xu,yu,xv,yv,tuv),...
                ttraj, z0(:), options);
        end
        [xtr_1ptcl, ytr_1ptcl] = deal(zz(:,1), zz(:,2));

        % stats
        timelen_valid(ip,jp) = find(isnan(xtr_1ptcl),1,'first');

        % 
        xtr(ip,jp,:) = xtr_1ptcl(it_save);
        ytr(ip,jp,:) = ytr_1ptcl(it_save);

    end
end
delete(gcp('nocreate'));
toc;

disp('Traj calculated.')
save(params_fnm, 'timelen_valid', "-append") 

%% save trajs (one single file for a snapshot)

for it = 1:ntsave

    [yrstr, dystr, hrstr] = get_timestr(tsave(it), yr_s); 
    savefnm = [save_dir 'trajs_' yrstr '_' dystr '_' hrstr '.nc'];
    if exist(savefnm,'file') 
        disp(['NC file already exists, skip : ', savefnm]);
        continue
    end  
    fprintf(1,'\nPtcls pos saving to: %s...\n\n', savefnm);

    [xtr_sv, ytr_sv] = deal(xtr(:,:,it), ytr(:,:,it));
    dim_name = {'Nx','Ny'};
    dim_length = [npx, npy];
    varname = {'xtr','ytr'};
    data = {xtr_sv, ytr_sv};
    dimNum_of_var = {[1,2], [1,2]};
    %
    global_att  = [ 'Ptcl trajs [km] saved at certain time(s); ' ...
        'init_fnm=' init_fnm '; params_fnm=' params_fnm ...
        '; day_s=' num2str(day_s) '; day_e=' num2str(day_s)];
    FUN_nc_easywrite_enhanced( savefnm, dim_name, dim_length,...
        varname, dimNum_of_var, data, global_att )
end

