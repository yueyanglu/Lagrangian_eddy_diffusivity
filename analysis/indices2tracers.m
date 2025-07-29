% 
% Convert the position indices (posIndPtcls{ibin,jbin} = {ptcles}) to tracer
%  evolutions, based on the tracers initially assigned to ptcls in each bin.
% 
% Need: 
%   Initial ptcl positions
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

%% dirs
ik = 1;
carry_al = 1:8;
ntr = numel(carry_al);
useWichInit = 2; % 1: use ref initial c&h for ptcls' properties at beginning of tloop
                 % 2: use the ptcl-based c&h for ptcls' properties (consecutive)

tloop = 6;

% ------- particl trajs
exp_dir = [campdir '/lagr_study/exp1'];
trajs_dir = [exp_dir '/trajs_bilinear/lp' num2str(tloop,'%02d')...
    '/Z' num2str(ik,'%02d') '/full/'];
ds_trajparams = load([trajs_dir 'params.mat']);

% ------- time do
yr_s = 21;
day_interv = 130;  
dt = 12/24;
day_s = (tloop-1)*120 + 1;
day_e = day_s + day_interv;
% 
t_do = day_s:dt:day_e;

% ------- init pos
init_fnm = ds_trajparams.init_fnm;
xin2d = ncread(init_fnm,'xin2d'); 
yin2d = ncread(init_fnm,'yin2d');
[npx, npy] = size(xin2d);
% initial time 
[yrstr0, dystr0, hrstr0] = get_timestr(t_do(1), yr_s); 
%     get_timestr(ds_trajparams.ttraj(1), ds_trajparams.yr_s); 

% ------- pos indices & CS bins
posind_dir = [exp_dir '/pos_ind/lp' num2str(tloop,'%02d') '/Z' num2str(ik,'%02d')];
posind0_fnm = [posind_dir '/posind__' yrstr0 '_' dystr0 '_' hrstr0 '.mat'];
ds_pos0 = load(posind0_fnm);
[xbins_bdry_km,ybins_bdry_km] = deal(ds_pos0.xbins_bdry_km,ds_pos0.ybins_bdry_km);
[nxbins, nybins] = deal(length(xbins_bdry_km)-1, length(ybins_bdry_km)-1);
fprintf(1,'pos_ind(t=0) readed from: %s\n',posind0_fnm);

%----- save dir
save_h_dir = [exp_dir '/h_new/lp' num2str(tloop,'%02d') '/Z' num2str(ik,'%02d')];
save_c_dir = [exp_dir '/c_new/lp' num2str(tloop,'%02d') '/Z' num2str(ik,'%02d')];
if ~exist(save_h_dir,'dir')
    mkdir(save_h_dir);
end
if ~exist(save_c_dir,'dir')
    mkdir(save_c_dir);
end
%% initial layer thickness & tracer conc

%----- CS grid (corresponding to the CS bins)
grid = build_grid_MOM(nxbins,nybins,[0 3840],[0 3840]);
cs_len = 1024/nxbins;

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

%----- initial h
ds_h0 = ncstruct(h0_fnm);
hPerBin0 = ds_h0.h(:,:,ik);
volPerBin = hPerBin0 .* grid.Ah; % [m^3]
fprintf(1,'h(t=0) readed from: %s\n',h0_fnm);

%----- initial tr conc (cell array of all tracers)
[trPerBin0, trPerPtcl0, trPerPtcl0_1d] = deal(cell(ntr,1));
ds_tr0 = ncstruct(tr0_fnm);
fprintf(1,'trs(t=0) %s readed from: %s\n',mat2str(carry_al), tr0_fnm);
for itr = 1:ntr
    varnm = ['tr' num2str(carry_al(itr))];
    trPerBin0{itr} = ds_tr0.(varnm)(:,:,ik);
end

%----- assign each ptcl with the initial h/c 
% 
% Also convert to 1d, to be consistent with the "pos index" which counts the
% ptcls based on their 1d index
% 
% All ptcls in a bin have the same h
hPerPtcl0 = perBin2perPtcl(hPerBin0,npx,npy);
hPerPtcl0_1d = reshape(hPerPtcl0, [npx*npy 1]); 
% All ptcls in a bin equally divide the total vol
volPerPtcl0 = perBin2perPtcl(volPerBin./ds_pos0.numptcls,npx,npy);
volPerPtcl0_1d = reshape(volPerPtcl0, [npx*npy 1]);
% All ptcls in a bin have the same tracer conc
for itr = 1:ntr
    trPerPtcl0{itr} = perBin2perPtcl(trPerBin0{itr},npx,npy);
    trPerPtcl0_1d{itr} = reshape(trPerPtcl0{itr}, [npx*npy 1]);
end

%% calc tracers
trPerPtcl_1d = trPerPtcl0_1d;
hPerPtcl_1d = hPerPtcl0_1d;
volPerPtcl_1d = volPerPtcl0_1d;

for it = 1:length(t_do)
    % time
    [yrstr, dystr, hrstr] = get_timestr(t_do(it), yr_s);
    %
    save_h_fnm = [save_h_dir '/h_snap__' yrstr '_' dystr '_' hrstr '.nc'];
    save_c_fnm = [save_c_dir '/tr__' yrstr '_' dystr '_' hrstr '.nc'];
    if exist(save_h_fnm,'file') && exist(save_c_fnm,'file')
        fprintf(1,'\nPtcl-based flds exist, skip: %s\n',save_h_fnm);
        continue
    end

    % pos indices
    posind_fnm = [posind_dir '/posind__' yrstr '_' dystr '_' hrstr '.mat'];
    ds_pos = load(posind_fnm);
    fprintf(1,'\npos_ind readed from: %s\n',posind_fnm);

    % goal
    [hPerBin, volPerBin] = deal(zeros(nxbins,nybins));
    trPerBin = cell(ntr,1); 
    trPerBin(:) = {zeros(nxbins,nybins)};

    % loop over each bin
    for jbin = 1:nxbins
        for ibin = 1:nybins
            % Linear indices of ptcls that are located in the bin(ib,jb)
            indx_temp = ds_pos.posIndPtcls{ibin,jbin};
            %
            if isempty(indx_temp)
                hPerBin(ibin,jbin) = NaN;
                volPerBin(ibin,jbin) = NaN;
                for itr = 1:ntr
                    trPerBin{itr}(ibin,jbin) = NaN;
                end
            else
                % h of the bin = average of h's carried by ptcls
                hPerBin(ibin,jbin) = mean(hPerPtcl_1d(indx_temp));
                % vol of bin = sum of vol carried by ptcls
                volPerBin(ibin,jbin) = sum(volPerPtcl_1d(indx_temp));
                % tr conc = h-weighted mean of c's = sum(c*h)/sum(h)
                for itr = 1:ntr
                    trPerBin{itr}(ibin,jbin) = sum( trPerPtcl_1d{itr}(indx_temp)...
                       .*hPerPtcl_1d(indx_temp) ) / sum(hPerPtcl_1d(indx_temp));
                    %trPerBin{itr}(ibin,jbin) = mean(trPerPtcl_1d{itr}(indx_temp));
                end
            end
        end % ibin
    end %jbin

    %------------------------- re-assign the quantities carried
    %                          by ptcls
%     if ifupdate
%         hPerPtcl_1d(indx_temp) = c(jc,ic,it,lp);
%         volPerPtcl_1d(indx_temp) = ch(jc,ic,it,lp);
%         trPerPtcl_1d(indx_temp) = cv(jc,ic,it,lp); %  / nptcl(jc,ic,it,lp)
%     end

    %----------- save h, vol
    dim_name = {'xh','yh'};
    dim_length = [nxbins, nybins];
    varname = {'h', 'vol'};
    data = {hPerBin, volPerBin};
    dimNum_of_var = {[1,2], [1,2]};
    global_att  = [ 'Euelerian h reconstructed from particles; ' ...
        'posind0_fnm=' posind0_fnm '; h0_fnm=' h0_fnm ...
        '; trajs_dir=' trajs_dir '; useWichInit=' num2str(useWichInit)];
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
        'posind0_fnm=' posind0_fnm '; tr0_fnm=' tr0_fnm ...
        '; h0_fnm=' h0_fnm '; trajs_dir=' trajs_dir];
    FUN_nc_easywrite_enhanced( save_c_fnm, dim_name, dim_length,...
        varname, dimNum_of_var, data, global_att )
    fprintf(1,'Ptcl-based c saved to: %s...\n', save_c_fnm);

end % t

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
