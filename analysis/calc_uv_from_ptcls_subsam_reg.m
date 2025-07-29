% 
% Calculate the coarse-grid u/v by subsampling the pre-calc trajs and Lagr velocities.
%  Note that u/v are defined on the coarse u-/v- grids!
% 
% Need: 
%   initial positions for the subsampled particles' init positions (.nc)
%   the indices of the subsampled particles (subsamp1.mat) 
%   trajectories of the original ptcls
%   Lagr velocities of the original ptcls
%   
% 
clear
homedir = getenv('HOME');
workdir = getenv('WORK');
campdir = getenv('CAMP');
addpath(genpath([homedir '/work_Lagr']));
addpath(genpath([homedir '/work_MOM']));
addpath(genpath([homedir '/mytoolbox']));
addpath(genpath([homedir '/MyFuncs']));

% + magic_args="params "
ik = 1;
expStr = '6'; % corresponds to subsampling ratio: 1, 0.5, 0.125...
subexpStr = num2str(subexpStrSh,'%02d'); % num2str(subexpStrSh,'%02d');
% -

% ------- indices of the subsampled ptcls
initmat_fnm = [campdir '/lagr_study/ptcl_initpos/region1_sub' expStr '_' subexpStr '.mat'];
ds = load(initmat_fnm);
indx_choose = ds.indx_choose_al; % indices of the subsampled ptcls in xtr_orig

for tloop = 1:6

% ------- original particl trajs (to be subsampled)
% 'exp2_01': sub ratio=0.5, randomly choosen - exp1
trajs_orig_dir = [campdir '/lagr_study/exp1/trajs_bilinear/lp' num2str(tloop,'%02d')...
    '/Z' num2str(ik,'%02d') '/full/'];
ulvl_orig_dir = [campdir '/lagr_study/exp1/ulvl/lp' num2str(tloop,'%02d')...
    '/Z' num2str(ik,'%02d') '/'];

% ------- time
yr_s = 21;
day_interv = 130;  
dt = 12/24;
day_s = (tloop-1)*120 + 1;
day_e = day_s + day_interv;
t_do = day_s:dt:day_e;
fprintf(1,'\nThis script does times D%d to D%d.\n\n',day_s, day_e);

save_dir = [campdir '/lagr_study/reg1_exp' expStr '_' subexpStr '/uv/lp' ...
    num2str(tloop,'%02d') '/Z' num2str(ik,'%02d')];
if ~exist(save_dir,'dir'); mkdir(save_dir); end

% % grid and bins

% --- HR grid on which ptcls are simulated
grid_dir = [workdir '/MOM6_exp/swm_spunup/'];  
[grid_h, ~, ~] = read_grid_MOM([grid_dir '']); 

%--- coarse bins
cells_in_bin = 32; % # of p-grid cells in one direction within a bin
% indices of the boundaries of bins (q-grid)
[xbins_bdry_id, ybins_bdry_id] = deal(1:cells_in_bin:grid_h.niq, 1:cells_in_bin:grid_h.njq);
if (xbins_bdry_id(end) ~= grid_h.niq) || (ybins_bdry_id(end) ~= grid_h.njq)
    warning('Coarse bins are not evenly distributed across the domain!')
end
% positions [km] of the boundaries of bins [nxbins_bdry-by-1]
[xbins_bdry_km, ybins_bdry_km] = deal(grid_h.lonq(xbins_bdry_id), grid_h.latq(ybins_bdry_id));
% # of P bins
[nxbins, nybins] = deal(length(xbins_bdry_id)-1, length(ybins_bdry_id)-1);

%----- CS grid (corresponding to the CS bins)
grid = build_grid_MOM(nxbins,nybins,grid_h.lonq([1 end]),grid_h.latq([1 end]));

% %

tic;
% 5s for one snap
for it = 1:length(t_do)

    [yrstr, dystr, hrstr] = get_timestr(t_do(it), yr_s);
    % --- savename
    save_fnm = [save_dir '/prog__' yrstr '_' dystr '_' hrstr '.nc'];
    if exist(save_fnm,'file')
        fprintf(1,'uv exist, so exit! \n%s\n',save_fnm);
        continue
    end

    % --- read trajs & ul
    traj_fnm = [trajs_orig_dir '/trajs__' yrstr '_' dystr '_' hrstr '.nc'];
    xtr_orig = ncread(traj_fnm,'xtr');
    ytr_orig = ncread(traj_fnm,'ytr');
    ulvl_fnm = [ulvl_orig_dir '/ulvl__' yrstr '_' dystr '_' hrstr '.nc'];
%     ulvl_fnm = traj_fnm;
    ul_orig = ncread(ulvl_fnm,'ul');
    vl_orig = ncread(ulvl_fnm,'vl');
    fprintf(1,'Lagr trajs readed from: %s\n',traj_fnm);
    fprintf(1,'Lagr vels readed from: %s\n',ulvl_fnm);
    % reshape
    xtr_orig = reshape(xtr_orig,[numel(xtr_orig) 1]);
    ytr_orig = reshape(ytr_orig,[numel(ytr_orig) 1]);
    ul_orig = reshape(ul_orig,[numel(ul_orig) 1]);
    vl_orig = reshape(vl_orig,[numel(vl_orig) 1]);

    % subsample (MUST be consistent with "subsample_initpos.ipynb")
    xtr = xtr_orig(indx_choose);
    ytr = ytr_orig(indx_choose);
    ul = ul_orig(indx_choose);
    vl = vl_orig(indx_choose);
    clearvars xtr_orig ytr_orig ul_orig vl_orig

    % --- average all ul within a u-bin; and all vl in a v-bin
    u_bin = zeros(grid.niu,grid.nju);
    v_bin = zeros(grid.niv,grid.njv);

    % --- loop over u-bins. Skip u_bin([1 end],:) due to solid boundaries.
    for iu = 2:grid.niu-1
        for ju = 1:grid.nju
            % bdry of each u bin: lon given by p-grids, lat by q-grids
            [wlon, elon] = deal(grid.lonh(iu-1), grid.lonh(iu));
            [slat, nlat] = deal(grid.latq(ju), grid.latq(ju+1));
            % Linear indices of ptcls within the cell at current inst.
            indx_temp = find(xtr>=wlon & xtr<elon & ytr>=slat & ytr<nlat );

            % average of ul
            if isempty(indx_temp)
                u_bin(iu,ju) = NaN;
            else
                u_bin(iu,ju) = mean(ul(indx_temp),'omitnan');
            end
        end
    end

    % --- loop over v-bins. Skip v_bin(:,[1 end]) due to solid boundaries.
    for iv = 1:grid.niv
        for jv = 2:grid.njv-1
            % bdry of each v bin: lon given by q-grids, lat by p-grids
            [wlon, elon] = deal(grid.lonq(iv), grid.lonq(iv+1));
            [slat, nlat] = deal(grid.lath(jv-1), grid.lath(jv));
            % Linear indices of ptcls within the cell at current inst.
            indx_temp = find(xtr>=wlon & xtr<elon & ytr>=slat & ytr<nlat );

            if isempty(indx_temp)
                v_bin(iv,jv) = NaN;
            else
                v_bin(iv,jv) = mean(vl(indx_temp),'omitnan');
            end
        end
    end
    
    %---------- fill NaN values in Eulerian flds caused by no-particles
    %u_bin = fillmissing2(u_bin,'linear');
    %v_bin = fillmissing2(v_bin,'linear');

    % ------- save
    dim_name = {'xh','yh','xq','yq'};
    dim_length = [grid.nih, grid.njh, grid.niq, grid.njq];
    varname = {'u', 'v'};
    data = {u_bin, v_bin};
    dimNum_of_var = {[3,2], [1,4]};
    global_att  = [ 'Euelerian u/v reconstructed from SUBSAMPLED particles; ' ...
        'traj_fnm=' traj_fnm '; ulvl_fnm=' ulvl_fnm '; tloop=' num2str(tloop) ...
        'initmat_fnm=' initmat_fnm];
    FUN_nc_easywrite_enhanced( save_fnm, dim_name, dim_length,...
        varname, dimNum_of_var, data, global_att )
    fprintf(1,'Ptcl-based u/v saved to: %s...\n\n', save_fnm);
end
toc;

end % lp

