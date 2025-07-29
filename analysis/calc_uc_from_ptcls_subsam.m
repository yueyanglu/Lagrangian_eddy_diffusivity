% 
% Calculate the coarse-grid <uc> by subsampling the pre-calc trajs and Lagr velocities.
%   Note that this does not use pos-indx since the indx is for P-bins but 
%   the fluxes here are for U/V-bins.
% 
% Need: 
%   initial positions for the subsampled particles' init positions (.nc)
%   the indices of the subsampled particles (subsamp1.mat) 
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
expStr = '4'; % corresponds to subsampling ratio: 1, 0.5, 0.125...
subexpStr = num2str(subexpStrSh,'%02d');
useWichInit = 2; % 1: use ref initial c&h for ptcls' properties at beginning of tloop
                 % 2: use the ptcl-based c&h for ptcls' properties (consecutive)

for tloop = 1:6

% ------- particl trajs & lagr vel
trajs_orig_dir = [campdir '/lagr_study/exp1/trajs_bilinear/lp' num2str(tloop,'%02d')...
    '/Z' num2str(ik,'%02d') '/full/'];
ulvl_orig_dir = [campdir '/lagr_study/exp1/ulvl/lp' num2str(tloop,'%02d')...
    '/Z' num2str(ik,'%02d') '/'];

% ------- time do
yr_s = 21;
day_interv = 130;  
dt = 12/24;
day_s = (tloop-1)*120 + 1;
day_e = day_s + day_interv;
t_do = day_s:dt:day_e;
% initial time 
[yrstr0, dystr0, hrstr0] = get_timestr(t_do(1), yr_s); 

% ------- indices of the subsampled ptcls
init_fnm = [campdir '/lagr_study/ptcl_initpos/ptclset1_sub' expStr '_' subexpStr '.nc'];
xin1d = ncread(init_fnm,'xin1d');
yin1d = ncread(init_fnm,'yin1d');
initmat_fnm = [campdir '/lagr_study/ptcl_initpos/sub' expStr '_' subexpStr '.mat'];
ds = load(initmat_fnm);
indx_choose = ds.indx_choose_al; % indices of the subsampled ptcls in xtr_orig

%----- save dir
exp_dir = [campdir '/lagr_study/exp' expStr '_' subexpStr '/'];
save_dir = [exp_dir 'uc/lp' num2str(tloop,'%02d') '/Z' num2str(ik,'%02d')];
if ~exist(save_dir,'dir'); mkdir(save_dir); end

%% ------- CS bins
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
cs_len = grid_h.nih/nxbins;

%% initial layer thickness & tr conc (for the unevenly distributed subsampled ptcls)

% FOR a certain time inverval!
if tloop == 1
    tr0_fnm = [campdir '/csflds_len' num2str(cs_len,'%02d') '/c/tr__' ...
        yrstr0 '_' dystr0 '_' hrstr0 '.nc'];
else % for successive time intervals
    if useWichInit == 1
        tr0_fnm = [campdir '/csflds_len' num2str(cs_len,'%02d') '/c/tr__' ...
            yrstr0 '_' dystr0 '_' hrstr0 '.nc'];
    elseif useWichInit == 2
        tr0_fnm = [exp_dir '/c/lp' num2str(tloop-1,'%02d') '/Z' num2str(ik,'%02d') ...
            '/tr__' yrstr0 '_' dystr0 '_' hrstr0 '.nc'];
    end
end

%----- initial tr conc (cell array of all tracers)
[trPerBin0, trPerPtcl0_1d] = deal(cell(ntr,1));
ds_tr0 = ncstruct(tr0_fnm);
fprintf(1,'trs(t=0) %s readed from: %s\n',mat2str(carry_al), tr0_fnm);
for itr = 1:ntr
    varnm = ['tr' num2str(carry_al(itr))];
    trPerBin0{itr} = ds_tr0.(varnm)(:,:,ik);
end

% %---------- interpolate hPerBin0 and trPerBin0, if they have NaN (missing) values
% hPerBin0 = fillmissing2(hPerBin0,'linear');
% for itr = 1:ntr
%     trPerBin0{itr} = fillmissing2(trPerBin0{itr},'linear');
% end
% 
%----- assign each ptcl with h/c (THIS IS THE DIFFERENT PART FROM original
% calculation)
fprintf(1,'\nAssigning properties to each ptcl...\n\n');
% All ptcls in a bin have the same tracer conc
for itr = 1:ntr
    trPerPtcl0_1d{itr} = perBin2perPtcl_gen(trPerBin0{itr},...
        xbins_bdry_km,ybins_bdry_km,xin1d,yin1d);
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

    %---- goal
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
                    c_ptcls = trPerPtcl0_1d{itr}(indx_temp);
                    uc_bin{itr}(iu,ju) = mean(ul_ptcls.*c_ptcls,'omitnan');
                end 
            end
        end
    end % iu

    % --- loop over v-bins. Skip v_bin(:,[1 end]) due to solid boundaries.
    fprintf(1,'Doing vhc...\n');
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
                    c_ptcls = trPerPtcl0_1d{itr}(indx_temp);
                    vc_bin{itr}(iv,jv) = mean(vl_ptcls.*c_ptcls,'omitnan');
                end 
            end
        end
    end % iv

    %---------- fill NaN values in Eulerian flds caused by no-particles
    for itr = 1:ntr
        uc_bin{itr} = fillmissing2(uc_bin{itr},'linear');
        vc_bin{itr} = fillmissing2(vc_bin{itr},'linear');
    end

    % ------- save
    dim_name = {'xh','yh','xq','yq'};
    dim_length = [grid.nih, grid.njh, grid.niq, grid.njq];
    for itr = 1:ntr
        varname = {'uc', 'vc'};
        data = {uc_bin{itr}, vc_bin{itr}};
        dimNum_of_var = {[3,2], [1,4]};
        global_att  = [ '<uc> reconstructed from SUBSAMPLED particles; ' ...
            'traj_orig_fnm=' traj_fnm '; ulvl_orig_fnm=' ulvl_fnm '; tloop=' num2str(tloop) ...
            '; C=' ['tr' num2str(carry_al(itr))] '; initmat_fnm=' initmat_fnm];
        
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

%% function
function PerPtcl1d = perBin2perPtcl_gen(fPerBin,xbins_bdry,ybins_bdry,xtr1d,ytr1d)
% Convert a property (e.g., h or tracer conc) per bin to property per particle
%   fPerBin [nxbin, nybin]: property like h and tr, per bin (Eulerian)
%   xbins_bdry [nxbins_bdry-by-1]: boundaries of CS bins
%   xtr1d: positions of ptcls [unit must be same with xbins_bdry]
% 
%   PerPtcl [nptcls-1]:     property per particle (Lagrangian)
% 
% This also works for unevenly distributed particles !!!!
% 

% # of bins in each direction
[nxbins, nybins] = size(fPerBin);
% # of ptcls 
nptcls = length(xtr1d);
if nptcls ~= length(ytr1d)
    error('xin1d must have same size with yin1d !!!')
end

% goal
PerPtcl1d = NaN*zeros(nptcls,1);

% loop over each CS bin, find ptcls in it, and assign the ptcls with fld
% proterty of the bin
for jbin = 1:nxbins
    for ibin = 1:nybins
        
        % bdry of each bin
        [wlon, elon] = deal(xbins_bdry(ibin), xbins_bdry(ibin+1));
        [slat, nlat] = deal(ybins_bdry(jbin), ybins_bdry(jbin+1));
        
        % Linear indices of ptcls in the bin
        indx_temp = (xtr1d>=wlon & xtr1d<elon & ytr1d>=slat & ytr1d<nlat);

        % ptcls within the bib have the same tracer conc
        if ~isempty(indx_temp)
            PerPtcl1d(indx_temp) = fPerBin(ibin,jbin);
        end
        
    end
end

% ---- if there is any value not 
if ~isempty(find(isnan(PerPtcl1d), 1))
    warning('Some particles are NOT located in any bin!!!');
end

end
