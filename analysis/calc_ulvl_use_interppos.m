% 
% Calc particles' Lagr velocities by interpolating the Eulerian vel based
% on particles positions (trajectories).
% 
% Need:
%   trajs (npx,npy)
%   Eulerian vel fields (nx, ny)
% Output:
%   Lagragian vel fields (npx, npy)
% 
clear
homedir = getenv('HOME');
workdir = getenv('WORK');
campdir = getenv('CAMP');
addpath(genpath([homedir '/work_Lagr']));
addpath(genpath([homedir '/work_MOM']));
addpath(genpath([homedir '/mytoolbox']));
addpath(genpath([homedir '/MyFuncs']));

%% 

ik = 1;

for tloop = 5:6
cs_len = 32;

% --- times
yr_s = 21;
day_interv = 130;  
dt = 12/24;
day_s = (tloop-1)*120 + 1;
day_e = day_s + day_interv;
t_al = day_s:dt:day_e;
nt_al = length(t_al);

% --- dirs
traj_dir = [campdir '/lagr_study/exp1_new/trajs_bilinear/lp' num2str(tloop,'%02d')...
    '/Z' num2str(ik,'%02d') '/full'];
uv_eule_dir = [campdir '/mom_ptemp/sol_prog'];

save_dir = [campdir '/lagr_study/exp1_new/ulvl/lp' num2str(tloop,'%02d') '/Z' ...
    num2str(ik,'%02d')];
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end

% --- grid
grid_dir = [workdir '/MOM6_exp/swm_spunup/'];  
[grid, ~, ~] = read_grid_MOM([grid_dir '']); 

% convert u/v [m/s] to [km/d]
factor = 3600*24/1e3; 
% coordinate 
[xu,yu] = deal(grid.lonq, grid.lath); 
[xv,yv] = deal(grid.lonh, grid.latq); 

%% interp Eule vel to get Lagr one
tic;
for it = 1:nt_al
    % times
    [yrstr, dystr, hrstr] = get_timestr(t_al(it), yr_s); 

    save_fnm = [save_dir '/ulvl__' yrstr '_' dystr '_' hrstr '.nc'];
    if exist(save_fnm,'file')
        fprintf(1,'\nLagr uv exist, so exit! \n%s\n',save_fnm);
        continue
    end
    
    %---- read trajs
    traj_fnm = [traj_dir '/trajs__' yrstr '_' dystr '_' hrstr '.nc'];
    xtr = ncread(traj_fnm,'xtr');
    ytr = ncread(traj_fnm,'ytr');
    [npx, npy] = size(xtr);
    fprintf(1,'Trajs readed from: %s\n',traj_fnm);

    %---- read Eule uv
    flx_fnm = [uv_eule_dir '/prog__' yrstr '_' dystr '_' hrstr '.nc'];
    u = ncread(flx_fnm,'u');
    v = ncread(flx_fnm,'v');
    % u and v [m/s] --> [km/d]
    [u_eule, v_eule] = deal(u(:,:,ik)*factor, v(:,:,ik)*factor);
    fprintf(1,'Eulerian vel readed from: %s\n',flx_fnm);

    %---- interp 
    ul = interp2(xu,yu,u_eule',xtr,ytr,'cubic');
    vl = interp2(xv,yv,v_eule',xtr,ytr,'cubic');
    % [km/d] --> [m/s]
    [ul, vl] = deal(ul./factor, vl./factor);

    %--- save
    dim_name = {'Nx','Ny'};
    dim_length = [npx, npy];
    varname = {'ul','vl'};
    data = {ul, vl};
    dimNum_of_var = {[1,2], [1,2]};
    global_att  = [ 'Ptcl Lagr vel [m/s] at instant time; ' ...
        '; lp=' num2str(tloop)];
    FUN_nc_easywrite_enhanced( save_fnm, dim_name, dim_length,...
        varname, dimNum_of_var, data, global_att )
    fprintf(1,'Ptcl Lagr vel saved to: %s...\n\n', save_fnm);
end
toc;

end % loop
