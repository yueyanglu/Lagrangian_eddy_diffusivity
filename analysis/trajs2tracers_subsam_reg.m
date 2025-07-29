% 
% Convert the particle trajectories to tracer fields by subsampling the
% original trajs
% Unlike the original process, we do not convert trajs to pos_indx here!
% 
% Input:
%   initial positions for the subsampled particles' init positions (.nc)
%   the indices of the subsampled particles (subsamp1.mat) 
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

% + magic_args="params"
ik = 1;
carry_al = 1:8;
ntr = numel(carry_al);
expStr = '6'; % corresponds to subsampling ratio: 1, 0.5, 0.125...
subexpStr = num2str(subexpStrSh, '%02d'); % num2str(subexpStrSh, '%02d'); 02 03,...

useWichInit = 2; % 1: use ref initial c&h for ptcls' properties at beginning of tloop
                 % 2: use the ptcl-based c&h for ptcls' properties (consecutive)

for tloop = 1:6

% ------- time do
yr_s = 21;
day_interv = 130;  
dt = 12/24;
day_s = (tloop-1)*120 + 1;
day_e = day_s + day_interv; 
t_do = day_s:dt:day_e;
% initial time of this loop
[yrstr0, dystr0, hrstr0] = get_timestr(t_do(1), yr_s); 

% ------- original particl trajs (to be subsampled)
% 'exp2_01': sub ratio=0.5, randomly choosen - exp1
trajs_orig_dir = [campdir '/lagr_study/exp1/trajs_bilinear/lp' num2str(tloop,'%02d')...
    '/Z' num2str(ik,'%02d') '/full/'];

% ------- indices of the subsampled ptcls
init_fnm = [campdir '/lagr_study/ptcl_initpos/region1_ptclset1_sub' expStr '_' subexpStr '.nc'];  
xin1d = ncread(init_fnm,'xin1d');
yin1d = ncread(init_fnm,'yin1d');
initmat_fnm = [campdir '/lagr_study/ptcl_initpos/region1_sub' expStr '_' subexpStr '.mat'];
ds = load(initmat_fnm);
indx_choose = ds.indx_choose_al; % indices of the subsampled ptcls in xtr_orig

% ------- save dir
exp_dir = [campdir '/lagr_study/reg1_exp' expStr '_' subexpStr '/'];
if ~exist(exp_dir,'dir');  mkdir(exp_dir); end
save_np_dir = [exp_dir '/nptcls/lp' num2str(tloop,'%02d') '/Z' num2str(ik,'%02d')];
save_h_dir = [exp_dir '/h/lp' num2str(tloop,'%02d') '/Z' num2str(ik,'%02d')];
save_c_dir = [exp_dir '/c/lp' num2str(tloop,'%02d') '/Z' num2str(ik,'%02d')];

if ~exist(save_np_dir,'dir'); mkdir(save_np_dir); end
if ~exist(save_h_dir,'dir'); mkdir(save_h_dir); end
if ~exist(save_c_dir,'dir'); mkdir(save_c_dir); end

% % CS bins

% --- HR grid on which ptcls are simulated
grid_dir = [workdir '/MOM6_exp/swm_spunup/'];  
[grid_h, ~, ~] = read_grid_MOM([grid_dir '']); 

%--- coarse bins, on which the tracer and h fields are defined
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

% + magic_args="initial layer thickness & tr conc (for the unevenly distributed subsampled ptcls)"
% FOR a certain time inverval!
if tloop == 1
    h0_fnm = [campdir '/csflds_len' num2str(cs_len,'%02d') '/h/h_snap__' ...
        yrstr0 '_' dystr0 '_' hrstr0 '.nc'];
    tr0_fnm = [campdir '/csflds_len' num2str(cs_len,'%02d') '/c/tr__' ...
        yrstr0 '_' dystr0 '_' hrstr0 '.nc'];
else % for successive time intervals
    if useWichInit == 1
        h0_fnm = [campdir '/csflds_len' num2str(cs_len,'%02d') '/h/h_snap__' ...
            yrstr0 '_' dystr0 '_' hrstr0 '.nc'];
        tr0_fnm = [campdir '/csflds_len' num2str(cs_len,'%02d') '/c/tr__' ...
            yrstr0 '_' dystr0 '_' hrstr0 '.nc'];
    elseif useWichInit == 2
        h0_fnm = [save_h_dir(1:end-6) num2str(tloop-1,'%02d') '/Z' num2str(ik,'%02d') ...
            '/h_snap__' yrstr0 '_' dystr0 '_' hrstr0 '.nc'];
        tr0_fnm = [save_c_dir(1:end-6) num2str(tloop-1,'%02d') '/Z' num2str(ik,'%02d') ...
            '/tr__' yrstr0 '_' dystr0 '_' hrstr0 '.nc'];
    end
end
% -

%----- initial h
fprintf(1,'h(t=0) reading from: %s...\n',h0_fnm);
ds_h0 = ncstruct(h0_fnm);
hPerBin0 = ds_h0.h(:,:,ik);

%----- initial tr conc (cell array of all tracers)
fprintf(1,'trs(t=0) %s reading from: %s...\n',mat2str(carry_al), tr0_fnm);
[trPerBin0] = deal(cell(ntr,1));
ds_tr0 = ncstruct(tr0_fnm);
for itr = 1:ntr
    varnm = ['tr' num2str(carry_al(itr))];
    trPerBin0{itr} = ds_tr0.(varnm)(:,:,ik);
end

%---------- interpolate hPerBin0 and trPerBin0, if they have NaN (missing) values
hPerBin0 = fillmissing2(hPerBin0,'linear');
for itr = 1:ntr
   trPerBin0{itr} = fillmissing2(trPerBin0{itr},'linear');
end

%----- assign each ptcl with h/c (THIS IS THE DIFFERENT PART FROM original
volPerBin = hPerBin0 .* grid.Ah; % [m^3]
% calculation)
fprintf(1,'\nAssigning properties to each ptcl...\n\n');
% All ptcls in a bin have the same h
hPerPtcl0_1d = perBin2perPtcl_gen(hPerBin0,xbins_bdry_km,ybins_bdry_km,xin1d,yin1d);
% All ptcls in a bin equally divide the total vol
volPerBin_norm = volPerBin./ds.numptcls_perbin; volPerBin_norm(isinf(volPerBin_norm)) = NaN;
volPerPtcl0_1d = perBin2perPtcl_gen(volPerBin_norm ...
    ,xbins_bdry_km,ybins_bdry_km,xin1d,yin1d);
% All ptcls in a bin have the same tracer conc
trPerPtcl0_1d = cell(ntr,1);
for itr = 1:ntr
    trPerPtcl0_1d{itr} = perBin2perPtcl_gen(trPerBin0{itr},...
        xbins_bdry_km,ybins_bdry_km,xin1d,yin1d);
end

% + magic_args="calc h & tracers"
trPerPtcl_1d = trPerPtcl0_1d;
hPerPtcl_1d = hPerPtcl0_1d;
volPerPtcl_1d = volPerPtcl0_1d;
% -

for it = 1:length(t_do)
    % time
    [yrstr, dystr, hrstr] = get_timestr(t_do(it), yr_s);
    %
    save_np_fnm = [save_np_dir '/nptcls__' yrstr '_' dystr '_' hrstr '.nc'];
    save_h_fnm = [save_h_dir '/h_snap__' yrstr '_' dystr '_' hrstr '.nc'];
    save_c_fnm = [save_c_dir '/tr__' yrstr '_' dystr '_' hrstr '.nc'];
    if exist(save_np_fnm,'file') && exist(save_h_fnm,'file') && exist(save_c_fnm,'file')
        fprintf(1,'\nPtcl-based flds exist, skip: %s\n',save_np_fnm);
        continue
    end

    % ------ read original trajs & SUBSAMPLE
    traj_fnm = [trajs_orig_dir '/trajs__' yrstr '_' dystr '_' hrstr '.nc'];
    xtr_orig = ncread(traj_fnm,'xtr');
    ytr_orig = ncread(traj_fnm,'ytr');
    fprintf(1,'Trajs readed from: %s\n',traj_fnm);
    % reshape (this is important since it affects the indices!)
    xtr_orig = reshape(xtr_orig,[numel(xtr_orig) 1]);
    ytr_orig = reshape(ytr_orig,[numel(ytr_orig) 1]);
    % subsample (MUST be consistent with "subsample_initpos.ipynb")
    xtr = xtr_orig(indx_choose);
    ytr = ytr_orig(indx_choose);
    clearvars xtr_orig ytr_orig
    
    % ------ goal
    [hPerBin, volPerBin, numPerBin] = deal(zeros(nxbins,nybins));
    trPerBin = cell(ntr,1); 
    trPerBin(:) = {zeros(nxbins,nybins)};

    % ---- loop over each bin
    % for 32*32 bins, 0.5 sub-ratio, ~3s
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
                numPerBin(ibin,jbin) = 0;
                hPerBin(ibin,jbin) = NaN;
                volPerBin(ibin,jbin) = NaN;
                for itr = 1:ntr
                    trPerBin{itr}(ibin,jbin) = NaN;
                end
            else
                % num of ptcls
                numPerBin(ibin,jbin) = numel(indx_temp);
                % h of the bin = average of h's carried by ptcls
                hPerBin(ibin,jbin) = mean(hPerPtcl_1d(indx_temp));
                % vol of bin = sum of vol carried by ptcls
                volPerBin(ibin,jbin) = sum(volPerPtcl_1d(indx_temp));
                % tr conc = h-weighted mean of c's = sum(c*h)/sum(h)
                for itr = 1:ntr
                    trPerBin{itr}(ibin,jbin) = sum( trPerPtcl_1d{itr}(indx_temp)...
                        .*hPerPtcl_1d(indx_temp) ) / sum(hPerPtcl_1d(indx_temp));
                    %  trPerBin{itr}(ibin,jbin) = mean(trPerPtcl_1d{itr}(indx_temp));
                end
            end
            
        end % ibin
    end % jbin

    %---------- fill NaN values in Eulerian flds caused by no-particles
    hPerBin = fillmissing2(hPerBin,'linear');
    volPerBin = fillmissing2(volPerBin,'linear');
    for itr = 1:ntr
        trPerBin{itr} = fillmissing2(trPerBin{itr},'linear');
    end

    %----------- save # of ptcls
    dim_name = {'xh','yh'};
    dim_length = [nxbins, nybins];
    varname = {'numptcls'};
    data = {numPerBin};
    dimNum_of_var = {[1,2]};
    global_att  = [ 'Euelerian num of ptcls reconstructed from SUBSAMPLED particles; ' ...
        'initmat_fnm=' initmat_fnm '; init_fnm=' init_fnm];
    FUN_nc_easywrite_enhanced( save_np_fnm, dim_name, dim_length,...
        varname, dimNum_of_var, data, global_att )
    fprintf(1,'\nNum of ptcls saved to: %s...\n', save_np_fnm);

    %----------- save h, vol
    dim_name = {'xh','yh'};
    dim_length = [nxbins, nybins];
    varname = {'h', 'vol'};
    data = {hPerBin, volPerBin};
    dimNum_of_var = {[1,2], [1,2]};
    global_att  = [ 'Euelerian h, vol reconstructed from SUBSAMPLED particles; ' ...
        'initmat_fnm=' initmat_fnm '; h0_fnm=' h0_fnm '; init_fnm=' init_fnm];
    FUN_nc_easywrite_enhanced( save_h_fnm, dim_name, dim_length,...
        varname, dimNum_of_var, data, global_att )
    fprintf(1,'Ptcl-based h saved to: %s...\n', save_h_fnm);

    %----------- save all tracers
    dim_name = {'xh','yh'};
    dim_length = [nxbins, nybins];
    varname = {cellstr(num2str(carry_al(:),'tr%d'))};
    varname = cat(1, varname{:});
    data = {trPerBin};
    data = cat(1, data{:});
    dimNum_of_var = cell(size(data)); 
    dimNum_of_var(:) = {[1,2]};
    global_att  = [ 'Eulerian tracer flds reconstructed from particles; ' ...
        'initmat_fnm=' initmat_fnm '; tr0_fnm=' tr0_fnm ...
        '; h0_fnm=' h0_fnm '; init_fnm=' init_fnm];
    FUN_nc_easywrite_enhanced( save_c_fnm, dim_name, dim_length,...
        varname, dimNum_of_var, data, global_att )
    fprintf(1,'Ptcl-based c saved to: %s...\n\n\n', save_c_fnm);

end

end % lp

% end % ensembles: sub01, sub02,...

% + magic_args="function"
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
% -

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
