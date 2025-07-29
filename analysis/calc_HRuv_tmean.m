clear
homedir = getenv('HOME');
workdir = getenv('WORK');
campdir = getenv('CAMP');
addpath(genpath([homedir '/work_Lagr']));
addpath(genpath([homedir '/work_MOM']));
addpath(genpath([homedir '/mytoolbox']));
addpath(genpath([homedir '/MyFuncs']));

%%
yr_s = 21;
[day_s, day_e, dt] = deal(1.0, 730.0, 6/24); 
t_al = day_s:dt:day_e;
nt_al = length(t_al);

grid_dir = [workdir '/MOM6_exp/swm_spunup/'];
[grid, ~, ~] = read_grid_MOM([grid_dir '']);

uv_dir = [campdir '/forc_uvh_sm33/prog_SM_nodecomp'];
save_fnm = [uv_dir '/prog__tm_D' num2str(t_al(1)) '_D' num2str(t_al(end)) '.nc'];

um = zeros(grid.niu,grid.nju);
vm = zeros(grid.niv,grid.njv);

for it = 1:nt_al
    [yrstr, dystr, hrstr] = get_timestr(t_al(it), yr_s);
    fnm = [uv_dir '/prog__' yrstr '_' dystr '_' hrstr '.nc'];
    fprintf(1,'u/v reading from: %s\n',fnm);
    u = ncread(fnm,'u');
    v = ncread(fnm,'v');

    um = um + u;
    vm = vm + v;
end

um = um ./ nt_al;
vm = vm ./ nt_al;

%% save

dim_name = {'xh','yh','xq','yq','zl'};
dim_length = [grid.nih, grid.njh, grid.niq, grid.njq,3];
varname = {'u', 'v'};
data = {um, vm};
dimNum_of_var = {[3,2,5], [1,4,5]};
global_att  = [ 'Time averaged u and v; ' ...
    'uv_dir=' uv_dir '; day_s=' num2str(day_s) '; day_e=' num2str(day_e) ...
    'dt=' num2str(dt)];
FUN_nc_easywrite_enhanced( save_fnm, dim_name, dim_length,...
    varname, dimNum_of_var, data, global_att )
fprintf(1,'Time mean u/v saved to: %s...\n\n', save_fnm);
