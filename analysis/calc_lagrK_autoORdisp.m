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
exp_dir = [campdir '/lagr_study/exp1'];
tloop = tloopSh;

% autocorrelation of lagr vel (1) or disp of trajs (2)
wichmethd = 1;
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
init_fnm = [campdir '/lagr_study/ptcl_initpos/ptclset1.nc'];
xin2d = ncread(init_fnm,'xin2d');
yin2d = ncread(init_fnm,'yin2d');
% for original set, reshape
[npx, npy] = size(xin2d);
xin1d = reshape(xin2d,[npx*npy 1]);
yin1d = reshape(yin2d,[npx*npy 1]);

% ----
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
trajs_dir = [exp_dir '/trajs_bilinear/lp' num2str(tloop,'%02d')...
    '/Z' num2str(ik,'%02d') '/full/'];
ds_trajparams = load([trajs_dir 'params.mat']);
ulvl_dir = [exp_dir '/ulvl/lp' num2str(tloop,'%02d') ...
    '/Z' num2str(ik,'%02d')];
% Eulerian spatial mean vel, only used for "auto" method
% uvm_dir = [exp_dir '/uv/full/Z' num2str(ik,'%02d')];
uvm_dir = [campdir '/forc_uvh_sm' num2str(cells_in_bin+1) '/prog_SM_nodecomp'];
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
[xtr_al, ytr_al, ul_al, vl_al] = deal(zeros(npx*npy,ntuv));
[um_t,u_t] = deal(zeros(grid_h.niu,grid_h.nju,ntuv));
[vm_t,v_t] = deal(zeros(grid_h.niv,grid_h.njv,ntuv));

for it = 1:ntuv
    [yrstr, dystr, hrstr] = get_timestr(tuv(it), yr_s);
    % there are nan's in xtr and vel
    traj_fnm = [trajs_dir '/trajs__' yrstr '_' dystr '_' hrstr '.nc'];
    xtr = ncread(traj_fnm,'xtr');
    ytr = ncread(traj_fnm,'ytr');
    ulvl_fnm = [ulvl_dir '/ulvl__' yrstr '_' dystr '_' hrstr '.nc'];
    ul = ncread(ulvl_fnm,'ul');
    vl = ncread(ulvl_fnm,'vl');
    uvm_fnm = [uvm_dir '/prog__' yrstr '_' dystr '_' hrstr '.nc'];
    um = ncread(uvm_fnm,'u');
    vm = ncread(uvm_fnm,'v');
    fprintf(1,'Trajs readed from: %s\n',traj_fnm);
    fprintf(1,'Lagr vel readed from: %s\n',ulvl_fnm);
    fprintf(1,'Eul vel-mean readed from: %s\n',uvm_fnm);
    %
    xtr_al(:,it) = reshape(xtr,[npx*npy 1]);
    ytr_al(:,it) = reshape(ytr,[npx*npy 1]);
    ul_al(:,it) = reshape(ul,[npx*npy 1]);
    vl_al(:,it) = reshape(vl,[npx*npy 1]);
    um_t(:,:,it) = um(:,:,ik);
    vm_t(:,:,it) = vm(:,:,ik);

    if wichmethd == 2
        uv_fnm = [uv_dir '/prog__' yrstr '_' dystr '_' hrstr '.nc'];
        u = ncread(uv_fnm,'u');
        v = ncread(uv_fnm,'v');
        u_t(:,:,it) = u(:,:,ik);
        v_t(:,:,it) = v(:,:,ik);
        fprintf(1,'Eul vel-full readed from: %s\n',uv_fnm);
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
% auto: ~22s per bin for 1024 ptcls
% disp: ~100s per bin for 1024 ptcls
for ibin = 1:nxbins
    fprintf('Doing bin i = %d from nxbins=%d \n', ibin, nxbins)
    for jbin = 1:nybins
        tic;

        nptcl = length(posInd0{ibin,jbin});

        %-- trajs, vel of particles orig from this bin [np,nttraj]
        [xtr,ytr] = deal(xtr_al(posInd0{ibin,jbin},:), ytr_al(posInd0{ibin,jbin},:));
        [ul,vl] = deal(ul_al(posInd0{ibin,jbin},:), vl_al(posInd0{ibin,jbin},:));

        %--
        if wichmethd == 1
            % ensemble-mean of the auto-covariance of each ptcl's vel, a function of
            %  lag time & bins

            %-- Lagr mean vel [np,nttraj], interp from Eulerian mean vel
            [ulm, vlm] = deal(zeros(size(ul)));
            factor = 3600*24/1e3; % [m/s] to [km/d]
            for it = 1:ntuv
                ulm(:,it) = interp2(xu,yu,um_t(:,:,it)',xtr(:,it),ytr(:,it),'cubic');
                vlm(:,it) = interp2(xv,yv,vm_t(:,:,it)',xtr(:,it),ytr(:,it),'cubic');
                [ulm(:,it), vlm(:,it)] = deal(ulm(:,it)./factor, vlm(:,it)./factor); % to m/s
            end
            [ule, vle] = deal(ul - ulm, vl - vlm);

            %-- convert nan's to 0
            ul(isnan(ul)) = 0;  vl(isnan(vl)) = 0;
            ule(isnan(ule)) = 0;  vle(isnan(vle)) = 0;

            % length of auto
            mlen = ntuv; %?
            mlen2d = zeros(nxbins,nybins);
            mlen2d(ibin,jbin) = mlen;

            % time lag 'tau': range from 0 to mlem-2
            for taut = 0:mlen-2
                % auto-covariance [m2/s2] for each lag time: full & eddy vel
                [cuu,cvv,cuv,cvu] = deal(zeros(1,nptcl));
                [cuu_e,cvv_e,cuv_e,cvu_e] =  deal(zeros(1,nptcl));
                for ip = 1:nptcl
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

            % ---- full dispersion [m^2]
            % displacements [m] relative to initial positions  [km]
            % note xtr/ytr sizes [Nptcl-nt]
            dx = xtr - repmat(xtr(:,1), [1,ntuv]);
            dy = ytr - repmat(ytr(:,1), [1,ntuv]);
            % ensamble-mean dispersions (spread from CoM) [m^2]
            dxx(ibin,jbin,:) = mean(dx.*dx,1,'omitnan') - mean(dx,1,'omitnan').^2;
            dyy(ibin,jbin,:) = mean(dy.*dy,1,'omitnan') - mean(dy,1,'omitnan').^2;
            dxy(ibin,jbin,:) = mean(dx.*dy,1,'omitnan') - mean(dx,1,'omitnan').*mean(dy,1,'omitnan');

            % ---- EFFT dispersions [m^2]
            % Initialize ptcls along trajs and find disp due to eddy vel.
            [dxf,dyf] = deal(zeros(size(xtr))); % [np,ntuv]
            % Eulerian eddy vel [m/s] [x,y,tuv]
            ue_t = u_t - um_t;
            ve_t = v_t - vm_t;
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
            %   Calculate dispersion in [m^2]
            dxx_efft(ibin,jbin,:) = mean(dxf.*dxf,1,'omitnan') - mean(dxf,1,'omitnan').^2;
            dyy_efft(ibin,jbin,:) = mean(dyf.*dyf,1,'omitnan') - mean(dyf,1,'omitnan').^2;
            dxy_efft(ibin,jbin,:) = mean(dxf.*dyf,1,'omitnan') - mean(dxf,1,'omitnan').*mean(dyf,1,'omitnan');

        end
        toc;
    end % nbin
end % nbin

%---------- save
if wichmethd == 1
    save(save_fnm,'ruu*','ruv*','rvu*','rvv*','day_interv','mlen2d',...
        '*_dir','dt','tloop','tuv');
elseif wichmethd == 2
    save(save_fnm,'dxx*','dxy*','dyy*','day_interv',...
        '*_dir','dt','tloop','tuv');
end

