%
% Calc the Lagrangian diffusivity using the particles trajs and/or velocities
%
% Need:
%   particle's trajectories
%   particle's Lagr velocities
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

%% params
ik = 1;
expStr = '3'; 
subexpStr = num2str(1, '%02d'); 
exp_dir = [campdir '/lagr_study/exp' expStr '_' subexpStr '/'];
tloop = tloopSh; % tloopSh

factor = 3600*24/1e3; % [m/s] to [km/d]

% ------- autocorrelation of lagr vel (1) or disp of trajs (2)
wichmethd = 2;
if wichmethd == 1
    wichmethdStr = 'auto';
elseif wichmethd == 2
    wichmethdStr = 'disp';
end
disp(['wichmethd = ' wichmethdStr]);

%% CS bins
% --- HR grid on which ptcls are simulated
grid_dir = [workdir '/MOM6_exp/swm_spunup/'];
[grid_h, ~, ~] = read_grid_MOM([grid_dir '']);
% [km]
[xu,yu] = deal(grid_h.lonq, grid_h.lath);
[xv,yv] = deal(grid_h.lonh, grid_h.latq);

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
grid = build_grid_MOM(nxbins,nybins,grid_h.lonq([1 end]),grid_h.latq([1 end]));

%% initial pos of particle sets, & divide them into each cs bins

% ------- indices of the subsampled ptcls
init_fnm = [campdir '/lagr_study/ptcl_initpos/ptclset1_sub' expStr '_' subexpStr '.nc'];
xin1d = ncread(init_fnm,'xin1d');
yin1d = ncread(init_fnm,'yin1d');
initmat_fnm = [campdir '/lagr_study/ptcl_initpos/sub' expStr '_' subexpStr '.mat'];
ds = load(initmat_fnm);
indx_choose = ds.indx_choose_al; % indices of the subsampled ptcls in xtr_orig

% ---- indices for the subsampled trajs that are within each CS bin
posInd0 = cell(nxbins,nybins);
for ibin = 1:nxbins
    for jbin = 1:nybins
        [wlon, elon] = deal(xbins_bdry_km(ibin), xbins_bdry_km(ibin+1));
        [slat, nlat] = deal(ybins_bdry_km(jbin), ybins_bdry_km(jbin+1));
        indx_temp = find(xin1d>=wlon & xin1d<elon & yin1d>=slat & yin1d<nlat );
        %
        if isempty(indx_temp)
            posInd0{ibin,jbin} = [];
        else
            posInd0{ibin,jbin} = indx_temp;
        end
    end
end

%% -------- for each interval, calc spreading of ptcls originating from each bin -------

% ------- particl trajs & lagr vel
trajs_orig_dir = [campdir '/lagr_study/exp1/trajs_bilinear/lp' num2str(tloop,'%02d')...
    '/Z' num2str(ik,'%02d') '/full/'];
ulvl_orig_dir = [campdir '/lagr_study/exp1/ulvl/lp' num2str(tloop,'%02d')...
    '/Z' num2str(ik,'%02d') '/'];

% Eulerian spatial mean vel, only used for "auto" method
uvm_fnm = [campdir '/forc_uvh_sm' num2str(cells_in_bin+1) ...
    '/prog_SM_nodecomp/prog__tm_D1_D730.nc'];

uv_dir = [campdir '/mom_ptemp/sol_prog'];

% ------- time steps for trajs & vels in this interval
yr_s = 21;
day_interv = 128; % originally 130, but we skip the last few steps
dt = 12/24;
day_s = (tloop-1)*120 + 1;
day_e = day_s + day_interv;
tuv = day_s:dt:day_e;
ntuv = length(tuv);
dtspan = 6/24;
tspan = day_s:dtspan:day_e;  % times at which trajs are calculated [d]

% ------- savename
save_dir = [exp_dir '/lagr_diffus/Z' num2str(ik,'%02d')];
if ~exist(save_dir,'dir'); mkdir(save_dir); end
save_fnm = [save_dir '/' wichmethdStr '__lp' num2str(tloop,'%02d') '.mat'];
if exist(save_fnm,'file')
    fprintf(1,'file exist, so exit! \n%s\n',save_fnm);
    exit
end

%----- read trajs&vel at all time in all bin
[xtr_al, ytr_al, ul_al, vl_al] = deal(zeros(ds.nptcl,ntuv));
[u_t,ue_t] = deal(zeros(grid_h.niu,grid_h.nju,ntuv));
[v_t,ve_t] = deal(zeros(grid_h.niv,grid_h.njv,ntuv));

um3d = ncread(uvm_fnm,'u');
vm3d = ncread(uvm_fnm,'v');
um = um3d(:,:,ik) * factor; % [km/d]
vm = vm3d(:,:,ik) * factor;
fprintf(1,'Eul vel-mean readed from: %s\n',uvm_fnm);

for it = 1:ntuv
    [yrstr, dystr, hrstr] = get_timestr(tuv(it), yr_s);
    % there are nan's in xtr and vel
    traj_fnm = [trajs_orig_dir '/trajs__' yrstr '_' dystr '_' hrstr '.nc'];
    xtr_orig = ncread(traj_fnm,'xtr');
    ytr_orig = ncread(traj_fnm,'ytr');
    ulvl_fnm = [ulvl_orig_dir '/ulvl__' yrstr '_' dystr '_' hrstr '.nc'];
    ul_orig = ncread(ulvl_fnm,'ul');
    vl_orig = ncread(ulvl_fnm,'vl');
    fprintf(1,'Orig trajs readed from: %s\n',traj_fnm);
    fprintf(1,'Lagr vel readed from: %s\n',ulvl_fnm);
    % reshape
    xtr_orig = reshape(xtr_orig,[numel(xtr_orig) 1]);
    ytr_orig = reshape(ytr_orig,[numel(ytr_orig) 1]);
    ul_orig = reshape(ul_orig,[numel(ul_orig) 1]);
    vl_orig = reshape(vl_orig,[numel(vl_orig) 1]);

    % subsample (MUST be consistent with "subsample_initpos.ipynb")
    xtr_al(:,it) = xtr_orig(indx_choose); % [km]
    ytr_al(:,it) = ytr_orig(indx_choose);
    ul_al(:,it) = ul_orig(indx_choose);   % [m/s]
    vl_al(:,it) = vl_orig(indx_choose);
    clearvars xtr_orig ytr_orig ul_orig vl_orig

    if wichmethd == 2 % disp
        uv_fnm = [uv_dir '/prog__' yrstr '_' dystr '_' hrstr '.nc'];
        u = ncread(uv_fnm,'u');
        v = ncread(uv_fnm,'v');
        u_t(:,:,it) = u(:,:,ik) * factor; % [km/d]
        v_t(:,:,it) = v(:,:,ik) * factor;
        fprintf(1,'Eul full vel readed from: %s\n',uv_fnm);
        % Eulerian eddy vel [km/d] 
        ue_t(:,:,it) = u_t(:,:,it) - um;
        ve_t(:,:,it) = v_t(:,:,it) - vm;
    end
end

% ------- target fields
if wichmethd == 1
    % velocity auto-covariance R [m2/s2]
    % calc from full & residual (eddy) vels
    [ruu,ruv,rvu,rvv] = deal(zeros(nxbins,nybins,ntuv));
    [ruu_e,ruv_e,rvu_e,rvv_e] = deal(zeros(nxbins,nybins,ntuv));
elseif wichmethd == 2
    % dispersion [m2]
    [dxx,dxy,dyy] = deal(zeros(nxbins,nybins,ntuv));
    [dxx_efft,dxy_efft,dyy_efft] = deal(zeros(nxbins,nybins,ntuv));
end

% -------------- loop over each bin --------------
% auto: ~22s per bin for 1024 ptcls per bin
% disp: ~100s per bin for 1024 ptcls
[np2d_use,mlen2d] = deal(zeros(nxbins,nybins));

for ibin = 1:nxbins
    fprintf('Doing bin i = %d from nxbins=%d \n', ibin, nxbins)
    for jbin = 1:nybins
        tic;

        np_bin = length(posInd0{ibin,jbin});

        %-- trajs, vel of particles orig from this bin [np,nttraj]
        [xtr,ytr] = deal(xtr_al(posInd0{ibin,jbin},:), ytr_al(posInd0{ibin,jbin},:));
        [ul,vl] = deal(ul_al(posInd0{ibin,jbin},:), vl_al(posInd0{ibin,jbin},:));
        Ltraj = size(xtr,2);

        %----------- deal with ptcls that hit bdry
        NotNaNmatrx = ~isnan(xtr.*ytr); %[np-t]
        % time leng of valid traj for each ptcl, min=1 max=nt [np-1]
        tlen1d = sum(NotNaNmatrx,2);
        
        % If tlen is too short, ignore the corresponding partc!
        tlen1d(tlen1d<Ltraj*.4) = NaN;
        [xtr(isnan(tlen1d),:),ytr(isnan(tlen1d),:)] = deal(NaN);
        [ul(isnan(tlen1d),:),vl(isnan(tlen1d),:)] = deal(NaN);
        % # of ptcls actually used for calc in this bin
        np2d_use(ibin,jbin) = length(find(~isnan(xtr(:,1))));

        %%----------- time len that all the rest of ptcls have experienced
        % applied to 'auto' method only
%         mlen = max( min(tlen1d(~isnan(tlen1d))), 2); % mlen>=Ltraj*.4
        mlen = Ltraj;
        mlen2d(ibin,jbin) = mlen;
        
        %--
        if wichmethd == 1
            % ensemble-mean of the auto-covariance of each ptcl's vel, a function of
            %  lag time & bins
            
            %-- Lagr mean vel [np,nttraj], interp from Eulerian mean vel
            [ulm, vlm] = deal(zeros(size(ul))); % [m/s]
            for it = 1:ntuv
                ulm(:,it) = interp2(xu,yu,um',xtr(:,it),ytr(:,it),'cubic'); % km/d
                vlm(:,it) = interp2(xv,yv,vm',xtr(:,it),ytr(:,it),'cubic');
                % to [m/s]
                [ulm(:,it), vlm(:,it)] = deal(ulm(:,it)./factor, vlm(:,it)./factor); 
            end
            [ule, vle] = deal(ul - ulm, vl - vlm); % m/s

            % time lag 'tau': range from 0 to mlem-2
            for taut = 0:mlen-2
                % auto-covariance [m2/s2] for each lag time: full & eddy vel
                [cuu,cvv,cuv,cvu] = deal(zeros(1,np_bin));
                [cuu_e,cvv_e,cuv_e,cvu_e] =  deal(zeros(1,np_bin));
                for ip = 1:np_bin
                    % full
                    CV = cov(ul(ip,1:mlen-taut),ul(ip,1+taut:mlen));
                    cuu(ip) = CV(1,2);
                    CV = cov(vl(ip,1:mlen-taut),vl(ip,1+taut:mlen));
                    cvv(ip) = CV(1,2);
                    CV = cov(ul(ip,1:mlen-taut),vl(ip,1+taut:mlen));
                    cuv(ip) = CV(1,2);
                    CV = cov(vl(ip,1:mlen-taut),ul(ip,1+taut:mlen));
                    cvu(ip) = CV(1,2);
                    % residual
                    CV = cov(ule(ip,1:mlen-taut),ule(ip,1+taut:mlen));
                    cuu_e(ip) = CV(1,2);
                    CV = cov(vle(ip,1:mlen-taut),vle(ip,1+taut:mlen));
                    cvv_e(ip) = CV(1,2);
                    CV = cov(ule(ip,1:mlen-taut),vle(ip,1+taut:mlen));
                    cuv_e(ip) = CV(1,2);
                    CV = cov(vle(ip,1:mlen-taut),ule(ip,1+taut:mlen));
                    cvu_e(ip) = CV(1,2);
                end

                % ensemble-mean over ptcls
                ruu(ibin,jbin,taut+1) = mean(cuu,'omitnan');
                rvv(ibin,jbin,taut+1) = mean(cvv,'omitnan');
                ruv(ibin,jbin,taut+1) = mean(cuv,'omitnan');
                rvu(ibin,jbin,taut+1) = mean(cvu,'omitnan');
                %
                ruu_e(ibin,jbin,taut+1) = mean(cuu_e,'omitnan');
                rvv_e(ibin,jbin,taut+1) = mean(cvv_e,'omitnan');
                ruv_e(ibin,jbin,taut+1) = mean(cuv_e,'omitnan');
                rvu_e(ibin,jbin,taut+1) = mean(cvu_e,'omitnan');
            end % tau

        elseif wichmethd == 2

            % ---- full dispersion [km^2]
            % displacements relative to initial positions  [km]
            % note xtr/ytr sizes [Nptcl-nt]
            dx = xtr - repmat(xtr(:,1), [1,ntuv]);
            dy = ytr - repmat(ytr(:,1), [1,ntuv]);
            % ensamble-mean dispersions (spread from CoM) [km^2]
            dxx(ibin,jbin,:) = mean(dx.*dx,1,'omitnan') - mean(dx,1,'omitnan').^2;
            dyy(ibin,jbin,:) = mean(dy.*dy,1,'omitnan') - mean(dy,1,'omitnan').^2;
            dxy(ibin,jbin,:) = mean(dx.*dy,1,'omitnan') - mean(dx,1,'omitnan').*mean(dy,1,'omitnan');

            % ---- EFFT dispersions [km^2]
            % Initialize ptcls along trajs and find disp due to eddy vel.
            [dxf,dyf] = deal(zeros(size(xtr))); % [np,ntuv]
            % Advect ptcls along trajs by eddy vel
            for nt = 2 : ntuv
                % ptcl pos [km] [2*nptcl-1]
                z0 = [xtr(:,nt-1)'; ytr(:,nt-1)'];
                % calc eddy disp along the trj
                clear zz
                ttf = tuv(nt-1:nt); % integration interval (one time step)
                options = odeset('RelTol',10^(-6),'AbsTol',10^(-6));
                % zz [nt-2*nptcl]
                [~,zz] = ode45(@(t,zz) ...
                    HamEqSolver_BiLin_CGrid(t,zz,ue_t,ve_t,xu,yu,xv,yv,tuv),...
                    ttf, z0(:), options);
                % x/ytr [nptcls-nttf]
                [xtr_temp, ytr_temp] = deal(zz(:,1:2:end)', zz(:,2:2:end)');

                % add the disp together to get pseudo-tracks
                dxf(:,nt) = dxf(:,nt-1) + (xtr_temp(:,end)-xtr_temp(:,1));
                dyf(:,nt) = dyf(:,nt-1) + (ytr_temp(:,end)-ytr_temp(:,1));
            end
            %   Calculate dispersion in [km^2]
            dxx_efft(ibin,jbin,:) = mean(dxf.*dxf,1,'omitnan') - mean(dxf,1,'omitnan').^2;
            dyy_efft(ibin,jbin,:) = mean(dyf.*dyf,1,'omitnan') - mean(dyf,1,'omitnan').^2;
            dxy_efft(ibin,jbin,:) = mean(dxf.*dyf,1,'omitnan') - mean(dxf,1,'omitnan').*mean(dyf,1,'omitnan');

        end
        toc;
    end % nbin
end % nbin

%---------- save
if wichmethd == 1
    save(save_fnm,'ruu*','ruv*','rvu*','rvv*','day_interv',...
        '*_dir','dt','tloop','tuv','mlen2d','np2d_use');
elseif wichmethd == 2
    save(save_fnm,'dxx*','dxy*','dyy*','day_interv',...
        '*_dir','dt','tloop','tuv','mlen2d','np2d_use');
end

