% 
% Calculate the coarse-grid uhc using pre-recorded trajs and Lagr vels.
%   Note that this does not use pos-indx since the indx is for P-bins but 
%   the fluxes here are for U/V-bins.
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
useWichInit = 2; % 1: use ref initial c&h for ptcls' properties at beginning of tloop
                 % 2: use the ptcl-based c&h for ptcls' properties (consecutive)

exp_dir = [campdir '/lagr_study/exp1'];

for tloop = 1:6

% ------- particl trajs & lagr vel
trajs_dir = [exp_dir '/trajs_bilinear/lp' num2str(tloop,'%02d')...
    '/Z' num2str(ik,'%02d') '/full/'];
ds_trajparams = load([trajs_dir 'params.mat']);
ulvl_dir = [exp_dir '/ulvl/lp' num2str(tloop,'%02d') ...
    '/Z' num2str(ik,'%02d')];

% ------- time do
yr_s = 21;
day_interv = 130;  
dt = 12/24;
day_s = (tloop-1)*120 + 1;
day_e = day_s + day_interv;
t_do = day_s:dt:day_e;

% ------- init pos
init_fnm = ds_trajparams.init_fnm;
xin2d = ncread(init_fnm,'xin2d'); 
yin2d = ncread(init_fnm,'yin2d');
[npx, npy] = size(xin2d);
% initial time 
[yrstr0, dystr0, hrstr0] = get_timestr(t_do(1), yr_s); 

% ------- CS bins and CS grid
posind_dir = [exp_dir '/pos_ind/lp' num2str(tloop,'%02d') '/Z' num2str(ik,'%02d')];
posind0_fnm = [posind_dir '/posind__' yrstr0 '_' dystr0 '_' hrstr0 '.mat'];
ds_pos0 = load(posind0_fnm);
[xbins_bdry_km,ybins_bdry_km] = deal(ds_pos0.xbins_bdry_km,ds_pos0.ybins_bdry_km);
[nxbins, nybins] = deal(length(xbins_bdry_km)-1, length(ybins_bdry_km)-1);
fprintf(1,'pos_ind(t=0) readed from: %s\n',posind0_fnm);
grid = build_grid_MOM(nxbins,nybins,xbins_bdry_km([1 end]),ybins_bdry_km([1 end]));

%----- save dir
save_dir = [exp_dir '/uhc/lp' num2str(tloop,'%02d') '/Z' num2str(ik,'%02d')];
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end

%% initial layer thickness & tracer conc
cs_len = 32;

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
        h0_fnm = [exp_dir '/h/full/Z' num2str(ik,'%02d') ...
            '/h_snap__' yrstr0 '_' dystr0 '_' hrstr0 '.nc'];
        tr0_fnm = [exp_dir '/c/full/Z' num2str(ik,'%02d') ...
            '/tr__' yrstr0 '_' dystr0 '_' hrstr0 '.nc'];
    end
end

%----- initial h
ds_h0 = ncstruct(h0_fnm);
hPerBin0 = ds_h0.h(:,:,ik);
fprintf(1,'h(t=0) readed from: %s\n',h0_fnm);

%----- initial tr conc (cell array of all tracers)
[trPerBin0, trPerPtcl0] = deal(cell(ntr,1));
ds_tr0 = ncstruct(tr0_fnm);
fprintf(1,'trs(t=0) %s readed from: %s\n',mat2str(carry_al), tr0_fnm);
for itr = 1:ntr
    varnm = ['tr' num2str(carry_al(itr))];
    trPerBin0{itr} = ds_tr0.(varnm)(:,:,ik);
end

%----- assign each ptcl with the initial h/c  
% All ptcls in a bin have the same h
hPerPtcl0 = perBin2perPtcl(hPerBin0,npx,npy);
% All ptcls in a bin have the same tracer conc
for itr = 1:ntr
    trPerPtcl0{itr} = perBin2perPtcl(trPerBin0{itr},npx,npy);
end

%% calc u*h*c
tic
for it = 1:length(t_do)
    % time
    [yrstr, dystr, hrstr] = get_timestr(t_do(it), yr_s);
    %
    save_fnm = [save_dir '/C05/uhc__' yrstr '_' dystr '_' hrstr '.nc'];
    if exist(save_fnm,'file')
        fprintf(1,'\nPtcl-based uhc exist, skip: %s\n',save_fnm);
        continue
    end

    % --- read trajs and ul
    traj_fnm = [trajs_dir '/trajs__' yrstr '_' dystr '_' hrstr '.nc'];
    xtr = ncread(traj_fnm,'xtr');
    ytr = ncread(traj_fnm,'ytr');
    ulvl_fnm = [ulvl_dir '/ulvl__' yrstr '_' dystr '_' hrstr '.nc'];
    ul = ncread(ulvl_fnm,'ul');
    vl = ncread(ulvl_fnm,'vl');
    fprintf(1,'Trajs readed from: %s\n',traj_fnm);
    fprintf(1,'Lagr vels readed from: %s\n',ulvl_fnm);

    % goal
    [uhc_bin, vhc_bin] = deal(cell(ntr,1));
    uhc_bin(:) = {zeros(grid.niu,grid.nju)};
    vhc_bin(:) = {zeros(grid.niv,grid.njv)};
    
    % --- loop over u-bins. Skip u_bin([1 end],:) due to solid boundaries.
    fprintf(1,'Doing uhc...\n');
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
                    uhc_bin{itr}(iu,ju) = NaN;
                end
            else
                % properties carried by ptcls [vectors]
                h_ptcls = hPerPtcl0(indx_temp);
                ul_ptcls = ul(indx_temp);
                % calc flux
                for itr = 1:ntr
                    c_ptcls = trPerPtcl0{itr}(indx_temp);
                    uhc_bin{itr}(iu,ju) = mean(ul_ptcls.*h_ptcls.*c_ptcls,'omitnan');
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
                    vhc_bin{itr}(iv,jv) = NaN;
                end
            else
                % properties carried by ptcls [vectors]
                h_ptcls = hPerPtcl0(indx_temp);
                vl_ptcls = vl(indx_temp);
                % calc flux
                for itr = 1:ntr
                    c_ptcls = trPerPtcl0{itr}(indx_temp);
                    vhc_bin{itr}(iv,jv) = mean(vl_ptcls.*h_ptcls.*c_ptcls,'omitnan');
                end 
            end
        end
    end % iv

    % ------- save
    dim_name = {'xh','yh','xq','yq'};
    dim_length = [grid.nih, grid.njh, grid.niq, grid.njq];
    for itr = 1:ntr
        varname = {'uhc', 'vhc'};
        data = {uhc_bin{itr}, vhc_bin{itr}};
        dimNum_of_var = {[3,2], [1,4]};
        global_att  = [ 'Euelerian <uhc> reconstructed from particles; ' ...
            'traj_fnm=' traj_fnm '; ulvl_fnm=' ulvl_fnm '; tloop=' num2str(tloop) ...
            '; C=' ['tr' num2str(carry_al(itr))] ];
        
        save_uhc_dir = [save_dir '/C' num2str(carry_al(itr),'%02d')];
        if ~exist(save_uhc_dir,'dir')
            mkdir(save_uhc_dir);
        end
        save_fnm = [save_uhc_dir '/uhc__' yrstr '_' dystr '_' hrstr '.nc'];
        FUN_nc_easywrite_enhanced( save_fnm, dim_name, dim_length,...
            varname, dimNum_of_var, data, global_att )
        fprintf(1,'Ptcl-based uhc saved to: %s...\n', save_fnm);
    end

end % it
toc

end

%% functions
function PerPtcl = perBin2perPtcl(PerBin,npx,npy)
% Convert a property (e.g., h or tracer conc) per bin to property per particle
%   PerCell [nxbin, nybin]: property like h and tr, per bin (Eulerian)
%   npx, npy: number of particles in each direction in the whole domain
%   PerPtcl [npx, npy]:     property per particle (Lagrangian)
% 
% ONLY works for evenly distributed particles !!!!

% # of bins in each direction
[nxbins, nybins] = size(PerBin);
% # of ptcls per bin in each direction
[npx_perbin, npy_perbin] = deal(npx/nxbins, npy/nybins);
if floor(npx_perbin) ~= npx_perbin
    error('# of ptcls per bin must be integer!!!')
end

% goal
PerPtcl = zeros(npx,npy);

% Inner loop must be along x direction, consistent with pos indices !!!!!!
for jbin = 1:nxbins
    for ibin = 1:nybins
        
        % Indices of ptcls within current cell
        [ii_ptcl,jj_ptcl] = deal( (ibin-1)*npx_perbin+1 : ibin*npx_perbin,...
            (jbin-1)*npy_perbin+1 : jbin*npy_perbin );
        
        % ptcls within the p-cell have the same tracer conc
        PerPtcl(ii_ptcl,jj_ptcl) = PerBin(ibin,jbin);
        
        clear ii_ptcl jj_ptcl
    end
end

end
