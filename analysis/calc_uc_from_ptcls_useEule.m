% 
% Calc particle-based c and uc, using Eulerian tracers to assign ptcl's
%  tracer fields.
% 
% Need: 
%   ptcles initial position
%   initial h and c
%   Particles' trajectories
%   Lagr velocities of particles
%   
clear
homedir = getenv('HOME');
workdir = getenv('WORK');
campdir = getenv('CAMP');
addpath(genpath([homedir '/work_Lagr']));
addpath(genpath([homedir '/work_MOM']));
addpath(genpath([homedir '/mytoolbox']));
addpath(genpath([homedir '/MyFuncs']));

%% dirs
ik = 1;
carry_al = 1:8;
ntr = numel(carry_al);

exp_dir = [campdir '/lagr_study/exp1_cat'];

% ------- dirs of HR Eule tracers
c_eule_dir = [campdir '/mom_ptemp/sol_tr_idl/'];

for tloop = 1:6

% ------- particl trajs & lagr vel
%trajs_dir = [exp_dir '/trajs_bilinear/lp' num2str(tloop,'%02d')...
%    '/Z' num2str(ik,'%02d') '/full/'];
%ulvl_dir = [exp_dir '/ulvl/lp' num2str(tloop,'%02d') ...
%    '/Z' num2str(ik,'%02d')];

trajs_dir = [exp_dir '/trajs/lp' num2str(tloop,'%02d') '/'];
ulvl_dir = [exp_dir '/ulvl/lp' num2str(tloop,'%02d') '/'];

% ------- time do
yr_s = 21;
day_interv = 130;  
dt = 12/24;
day_s = (tloop-1)*120 + 1;
day_e = day_s + day_interv;
t_do = day_s:dt:day_e;

% --- HR grid on which ptcls are simulated
grid_dir = [workdir '/MOM6_exp/swm_spunup/'];  
[grid_h, ~, ~] = read_grid_MOM([grid_dir '']); 

%--- coarse bins
cells_in_bin = 32; % # of p-grid cells in one direction within a bin
% indices of the boundaries of bins (q-grid)
[xbins_bdry_id, ybins_bdry_id] = deal(1:cells_in_bin:grid_h.niq, 1:cells_in_bin:grid_h.njq);
% # of bins
[nxbins_bdry, nybins_bdry] = deal(length(xbins_bdry_id), length(ybins_bdry_id));
[nxbins, nybins] = deal(nxbins_bdry-1, nybins_bdry-1);
% CS grid 
grid = build_grid_MOM(nxbins,nybins,grid_h.lonq([1 end]),grid_h.latq([1 end]));
cs_len = 1024/nxbins;

%----- save dir
save_dir = [exp_dir '/uc_interp/lp' num2str(tloop,'%02d') '/Z' num2str(ik,'%02d')];
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end

%% calc u*h*c
tic
for it = 1:length(t_do)
    % time
    [yrstr, dystr, hrstr] = get_timestr(t_do(it), yr_s);
    %
    save_fnm = [save_dir '/C05/uc__' yrstr '_' dystr '_' hrstr '.nc'];
    if exist(save_fnm,'file')
        fprintf(1,'\nPtcl-based uc exist, skip: %s\n',save_fnm);
        continue
    end

    % --- read trajs and ul
    traj_fnm = [trajs_dir '/trajs__' yrstr '_' dystr '_' hrstr '.nc'];
    xtr = ncread(traj_fnm,'xtr');
    ytr = ncread(traj_fnm,'ytr');
    nptcls = numel(xtr);
    xtr = reshape(xtr,[nptcls 1]);
    ytr = reshape(ytr,[nptcls 1]);

    ulvl_fnm = [ulvl_dir '/ulvl__' yrstr '_' dystr '_' hrstr '.nc'];
    ul = ncread(ulvl_fnm,'ul');
    vl = ncread(ulvl_fnm,'vl');
    fprintf(1,'Trajs readed from: %s\n',traj_fnm);
    fprintf(1,'Lagr vels readed from: %s\n',ulvl_fnm);

    % ------ use HR Eule tracer to get tracer per ptcl
    cEule_fnm = [c_eule_dir '/tr__' yrstr '_' dystr '_' hrstr '.nc'];
    ds_cEule = ncstruct(cEule_fnm);
    fprintf(1,'Eulerian tr readed from: %s\n',cEule_fnm);
    %
    trPerPtcl_1d = cell(ntr,1);
    trPerPtcl_1d(:) = {zeros(size(xtr))};
    for itr = 1:ntr
        varname = ['tr' num2str(carry_al(itr))];
        c_eule = ds_cEule.(varname)(:,:,ik);
        trPerPtcl_1d{itr} = interp2(grid_h.lonh,grid_h.lath,c_eule',xtr,ytr,'cubic');
    end

    % ------ goal
    [uc_bin, vc_bin] = deal(cell(ntr,1));
    uc_bin(:) = {zeros(grid.niu,grid.nju)};
    vc_bin(:) = {zeros(grid.niv,grid.njv)};
    
    % --- loop over u-bins. Skip u_bin([1 end],:) due to solid boundaries.
    fprintf(1,'Doing uc...\n');
    for iu = 2:grid.niu-1
        for ju = 1:grid.nju
            % bdry of each u bin: lon given by p-grids, lat by q-grids
            [wlon, elon] = deal(grid.lonh(iu-1), grid.lonh(iu));
            [slat, nlat] = deal(grid.latq(ju), grid.latq(ju+1));
            % Linear indices of ptcls within the bin at current inst.
            indx_temp = find(xtr>=wlon & xtr<elon & ytr>=slat & ytr<nlat );

            % average of ptcls in the bin
            if isempty(indx_temp)
                for itr = 1:ntr
                    uc_bin{itr}(iu,ju) = NaN;
                end
            else
                % properties carried by ptcls [vectors]
                ul_ptcls = ul(indx_temp);
                % calc flux
                for itr = 1:ntr
                    c_ptcls = trPerPtcl_1d{itr}(indx_temp);
                    uc_bin{itr}(iu,ju) = mean(ul_ptcls.*c_ptcls,'omitnan');
                end 
            end
        end
    end % iu

    % --- loop over v-bins. Skip v_bin(:,[1 end]) due to solid boundaries.
    fprintf(1,'Doing vc...\n');
    for iv = 1:grid.niv
        for jv = 2:grid.njv-1
            % bdry of each v bin: lon given by q-grids, lat by p-grids
            [wlon, elon] = deal(grid.lonq(iv), grid.lonq(iv+1));
            [slat, nlat] = deal(grid.lath(jv-1), grid.lath(jv));
            % Linear indices of ptcls within the cell at current inst.
            indx_temp = find(xtr>=wlon & xtr<elon & ytr>=slat & ytr<nlat );

            if isempty(indx_temp)
                for itr = 1:ntr
                    vc_bin{itr}(iv,jv) = NaN;
                end
            else
                % properties carried by ptcls [vectors]
                vl_ptcls = vl(indx_temp);
                % calc flux
                for itr = 1:ntr
                    c_ptcls = trPerPtcl_1d{itr}(indx_temp);
                    vc_bin{itr}(iv,jv) = mean(vl_ptcls.*c_ptcls,'omitnan');
                end 
            end
        end
    end % iv

    % ------- save
    dim_name = {'xh','yh','xq','yq'};
    dim_length = [grid.nih, grid.njh, grid.niq, grid.njq];
    for itr = 1:ntr
        varname = {'uc', 'vc'};
        data = {uc_bin{itr}, vc_bin{itr}};
        dimNum_of_var = {[3,2], [1,4]};
        global_att  = [ 'Euelerian <uc> reconstructed from particles (interp from Eule c); ' ...
            'traj_fnm=' traj_fnm '; ulvl_fnm=' ulvl_fnm '; tloop=' num2str(tloop) ...
            '; C=' ['tr' num2str(carry_al(itr))] ];
        
        save_uc_dir = [save_dir '/C' num2str(carry_al(itr),'%02d')];
        if ~exist(save_uc_dir,'dir')
            mkdir(save_uc_dir);
        end
        save_fnm = [save_uc_dir '/uc__' yrstr '_' dystr '_' hrstr '.nc'];
        FUN_nc_easywrite_enhanced( save_fnm, dim_name, dim_length,...
            varname, dimNum_of_var, data, global_att )
        fprintf(1,'Ptcl-based uc saved to: %s...\n', save_fnm);
    end

end % it
toc

end
