% 
% Convert particle trajetories to an cell array. Each entry of the array is
% the indices of the particles that locate in the bin correponding to the
% entry.
% This script perform the calculation at an instant time.
% 
% posIndPtcls{ibin,jbin} = {indices of ptcles} at the instant.
% 
clear
homedir = getenv('HOME');
workdir = getenv('WORK');
campdir = getenv('CAMP');
addpath(genpath([homedir '/work_Lagr']));
addpath(genpath([homedir '/work_MOM']));
addpath(genpath([homedir '/mytoolbox']));
addpath(genpath([homedir '/MyFuncs']));

%% grid and bins

% --- HR grid on which ptcls are simulated
grid_dir = [workdir '/MOM6_exp/swm_spunup/'];  
[grid, ~, ~] = read_grid_MOM([grid_dir '']); 

% --- configure the bins (normally consistent with eddy scale)
% the bin can be HR grid cells! 
% See "plot_cs_bins.ipynb"

%--- coarse bins
cells_in_bin = 32; % # of p-grid cells in one direction within a bin

% indices of the boundaries of bins (q-grid)
[xbins_bdry_id, ybins_bdry_id] = deal(1:cells_in_bin:grid.niq, 1:cells_in_bin:grid.njq);
if (xbins_bdry_id(end) ~= grid.niq) || (ybins_bdry_id(end) ~= grid.njq)
    warning('Coarse bins are not evenly distributed across the domain!')
end

% positions [km] of the boundaries of bins [nxbins_bdry-by-1]
[xbins_bdry_km, ybins_bdry_km] = deal(grid.lonq(xbins_bdry_id), grid.latq(ybins_bdry_id));

% # of bins
[nxbins_bdry, nybins_bdry] = deal(length(xbins_bdry_id), length(ybins_bdry_id));
[nxbins, nybins] = deal(nxbins_bdry-1, nybins_bdry-1);

%% params 
ik = 1;
tloop = 2;

yr_s = 21;
day_interv = 130;  
dt = 12/24;
day_s = (tloop-1)*120 + 1;
day_e = day_s + day_interv;
t_do = day_s:dt:day_e;

save_dir = [campdir '/lagr_study/pos_ind/lp' num2str(tloop,'%02d') '/Z' num2str(ik,'%02d')];
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end

tic;
for it = 1:length(t_do)

    [yrstr, dystr, hrstr] = get_timestr(t_do(it), yr_s);
    % --- read trajs
    traj_fnm = [campdir '/lagr_study/trajs_bilinear/lp' num2str(tloop,'%02d') ...
        '/Z' num2str(ik,'%02d') '/full/trajs__' yrstr '_' dystr '_' hrstr '.nc'];
    xtr = ncread(traj_fnm,'xtr');
    ytr = ncread(traj_fnm,'ytr');
    [npx, npy] = size(xtr);
    fprintf(1,'Trajs readed from: %s\n',traj_fnm);

    % --- save
    save_fnm = [save_dir '/posind__' yrstr '_' dystr '_' hrstr '.mat'];
    if exist(save_fnm,'file')
        fprintf(1,'Pos indices exist, so exit exc! \n%s\n',save_fnm);
        continue
    end

    %% calc
    posIndPtcls = cell(nxbins,nybins);
    numptcls = zeros(nxbins,nybins);

    % ---- reshape (this is important since it affects the indices!)
    xtr = reshape(xtr,[npx*npy 1]);
    ytr = reshape(ytr,[npx*npy 1]);

    % ---- count
    % for 32*32 bins, ~3.75s
    for ibin = 1:nxbins

%         fprintf('Doing bin i = %d from nxbins=%d \n', ibin, nxbins)

        for jbin = 1:nybins
            % bdry of each bin
            [wlon, elon] = deal(xbins_bdry_km(ibin), xbins_bdry_km(ibin+1));
            [slat, nlat] = deal(ybins_bdry_km(jbin), ybins_bdry_km(jbin+1));

            % Linear indices of ptcls within the cell at current inst.
            indx_temp = find(xtr>=wlon & xtr<elon & ytr>=slat & ytr<nlat );

            % Calc num of ptcls, vol and concentraions
            if isempty(indx_temp)
                posIndPtcls{ibin,jbin} = [];
            else
                posIndPtcls{ibin,jbin} = indx_temp;
            end
            numptcls(ibin,jbin) = numel(indx_temp);

        end
    end

    % ---- save
    save(save_fnm,'posIndPtcls','numptcls','traj_fnm','xbins_bdry_km','ybins_bdry_km');
    fprintf(1,'Pos indices of ptcls saved to: \n%s\n',save_fnm);
end % t
toc;
