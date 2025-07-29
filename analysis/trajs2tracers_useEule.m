% 
% Convert the particle trajectories to tracer fields by interpolating the
%   Eulerian HR tracer onto each particle position
% Unlike the original process, we do not convert trajs to pos_indx here!
% 
% Input:
%   tracers initially assigned to ptcls in each bin.
% 
clear
homedir = getenv('HOME');
workdir = getenv('WORK');
campdir = getenv('CAMP');
addpath(genpath([homedir '/work_Lagr']));
addpath(genpath([homedir '/work_MOM']));
addpath(genpath([homedir '/mytoolbox']));
addpath(genpath([homedir '/MyFuncs']));

%% params
ik = 1;
carry_al = 1:8;
ntr = numel(carry_al);

for tloop = 1:6

% ------- time do
yr_s = 21;
day_interv = 130;  
dt = 12/24;
day_s = (tloop-1)*120 + 1;
day_e = day_s + day_interv; 
t_do = day_s:dt:day_e;

% ------- original particl trajs (to be subsampled)
trajs_orig_dir = [campdir '/lagr_study/exp1_cat/trajs/lp' num2str(tloop,'%02d')  '/']; % /full
%trajs_orig_dir = [campdir '/lagr_study/exp1_cat/trajs/lp' num2str(tloop,'%02d')...
%    '/Z' num2str(ik,'%02d') '/']; % /full

% ------- save dir
exp_dir = [campdir '/lagr_study/exp1_cat/'];
save_c_dir = [exp_dir '/c_interp/lp' num2str(tloop,'%02d') '/Z' num2str(ik,'%02d')];

if ~exist(save_c_dir,'dir'); mkdir(save_c_dir); end

% ------- dirs of HR Eule tracers
c_eule_dir = [campdir '/mom_ptemp/sol_tr_idl/'];

%% CS bins

% --- HR grid on which ptcls are simulated
grid_dir = [workdir '/MOM6_exp/swm_spunup/'];  
[grid_h, ~, ~] = read_grid_MOM([grid_dir '']); 

%--- coarse bins
cells_in_bin = 32; % # of p-grid cells in one direction within a bin
% indices of the boundaries of bins (q-grid)
[xbins_bdry_id, ybins_bdry_id] = deal(1:cells_in_bin:grid_h.niq, 1:cells_in_bin:grid_h.njq);
% positions [km] of the boundaries of bins [nxbins_bdry-by-1]
[xbins_bdry_km, ybins_bdry_km] = deal(grid_h.lonq(xbins_bdry_id), grid_h.latq(ybins_bdry_id));
% # of bins
[nxbins_bdry, nybins_bdry] = deal(length(xbins_bdry_id), length(ybins_bdry_id));
[nxbins, nybins] = deal(nxbins_bdry-1, nybins_bdry-1);

%----- CS grid 
grid = build_grid_MOM(nxbins,nybins,grid_h.lonq([1 end]),grid_h.latq([1 end]));
cs_len = 1024/nxbins;

%% calc tracers

for it = 1:length(t_do)
    % time
    [yrstr, dystr, hrstr] = get_timestr(t_do(it), yr_s);
    %
    save_c_fnm = [save_c_dir '/tr__' yrstr '_' dystr '_' hrstr '.nc'];
    if exist(save_c_fnm,'file')
        fprintf(1,'\nPtcl-based flds exist, skip: %s\n',save_c_fnm);
        continue
    end

    % ------ read original trajs
    traj_fnm = [trajs_orig_dir '/trajs__' yrstr '_' dystr '_' hrstr '.nc'];
    xtr = ncread(traj_fnm,'xtr');
    ytr = ncread(traj_fnm,'ytr');
    nptcls = numel(xtr);
    xtr = reshape(xtr,[nptcls 1]);
    ytr = reshape(ytr,[nptcls 1]);
    fprintf(1,'Trajs readed from: %s\n',traj_fnm);

    % ------ read HR Eule tracer
    cEule_fnm = [c_eule_dir '/tr__' yrstr '_' dystr '_' hrstr '.nc'];
    ds_cEule = ncstruct(cEule_fnm);
    fprintf(1,'Eulerian tr readed from: %s\n',cEule_fnm);
    
    % ------ interp to get tracer per ptcl
    trPerPtcl_1d = cell(ntr,1);
    trPerPtcl_1d(:) = {zeros(size(xtr))};
    for itr = 1:ntr
        varname = ['tr' num2str(carry_al(itr))];
        c_eule = ds_cEule.(varname)(:,:,ik);
        trPerPtcl_1d{itr} = interp2(grid_h.lonh,grid_h.lath,c_eule',xtr,ytr,'cubic');
    end

    % ------ average interpolated ptcl' tracers in each cs bin
    trPerBin = cell(ntr,1); 
    trPerBin(:) = {zeros(nxbins,nybins)};

    for ibin = 1:nxbins
%         fprintf('Doing bin i = %d from nxbins=%d \n', ibin, nxbins)
        for jbin = 1:nybins
            % bdry of each bin
            [wlon, elon] = deal(xbins_bdry_km(ibin), xbins_bdry_km(ibin+1));
            [slat, nlat] = deal(ybins_bdry_km(jbin), ybins_bdry_km(jbin+1));
            % Linear indices of ptcls within the cell at current inst.
            indx_temp = find(xtr>=wlon & xtr<elon & ytr>=slat & ytr<nlat );

            % compute 
            if isempty(indx_temp)
                for itr = 1:ntr
                    trPerBin{itr}(ibin,jbin) = NaN;
                end
            else
                for itr = 1:ntr
                    trPerBin{itr}(ibin,jbin) = mean(trPerPtcl_1d{itr}(indx_temp));
                end
            end

        end % ibin
    end % jbin

    %---------- fill NaN values in Eulerian flds caused by no-particles
    for itr = 1:ntr
        trPerBin{itr} = fillmissing2(trPerBin{itr},'linear');
    end

    %----------- save all tracers
    dim_name = {'xh','yh'};
    dim_length = [nxbins, nybins];
    varname = {cellstr(num2str(carry_al(:),'tr%d'))};
    varname = cat(1, varname{:});
    data = {trPerBin};
    data = cat(1, data{:});
    dimNum_of_var = cell(size(data)); 
    dimNum_of_var(:) = {[1,2]};
    global_att  = [ 'Eulerian tracer flds reconstructed from particles with Eule tr; ' ...
        'cEule_fnm=' cEule_fnm ];
    FUN_nc_easywrite_enhanced( save_c_fnm, dim_name, dim_length,...
        varname, dimNum_of_var, data, global_att )
    fprintf(1,'Ptcl-based c saved to: %s...\n', save_c_fnm);

end

end % lp

