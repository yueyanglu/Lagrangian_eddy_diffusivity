function calc_eddyflux_ptcls(tloop)
%Calculate tracer eddy flux = <uhc>-<u><h><c>
% 
%
homedir = getenv('HOME');
campdir = getenv('CAMP');

addpath(genpath([homedir '/work_MOM']));
addpath(genpath([homedir '/mytoolbox']));
addpath(genpath([homedir '/MyFuncs']));

%%
% ----- which tracer
ik = 1;
tloop = 1;
carry_al = [1:8]; % 2 3 4 5 6 7 8
ntr = numel(carry_al);
fprintf(1,'CS-graining tracers: %s ...\n',mat2str(carry_al));

%---------------- dirs
exp_dir = [campdir '/lagr_study/'];  
uhc_dir = [exp_dir 'uhc/lp' num2str(tloop,'%02d') '/Z' num2str(ik,'%02d')];
uv_dir = [exp_dir 'uv/lp' num2str(tloop,'%02d') '/Z' num2str(ik,'%02d')];
h_dir = [exp_dir 'h/full/Z' num2str(ik,'%02d')];
tr_dir = [exp_dir 'c/full/Z' num2str(ik,'%02d')];

%-------- filter len, time interval
yr_s = 21;
day_interv = 130;  
dt = 12/24;
day_s = (tloop-1)*120 + 1;
day_e = day_s + day_interv;
t_al = day_s:dt:day_e;
nt_al = length(t_al);

% save dir
save_varname_fu = 'fu'; save_varname_fv = 'fv';
save_file_preStr = 'flx__';
save_dir = [exp_dir '/eddyflux/Z' num2str(ik,'%02d')];
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end

%%
for it = 1:nt_al
    
    %----------------------- current time
    nday = t_al(it);
    [yrstr, dystr, hrstr] = get_timestr(nday, yr_s); % for tend
    fprintf(1,'\n time = Y%s-D%s-H%s (it=%d) of %d snapshots ...\n',...
        yrstr,dystr,hrstr,it,nt_al);

    %----------------------- read <u> and <h>
    h_fnm = [h_dir '/h_snap__' yrstr '_' dystr '_' hrstr '.nc'];
    uv_fnm = [uv_dir '/prog__' yrstr '_' dystr '_' hrstr '.nc'];
    h2d = ncread(h_fnm, 'h');
    u2d = ncread(uv_fnm, 'u');
    v2d = ncread(uv_fnm, 'v');
    tr_fnm = [tr_dir '/tr__' yrstr '_' dystr '_' hrstr '.nc'];

    %----------------------- calc eddy flux
    for itr = 1:ntr

        wichtr = carry_al(itr);
        varname = ['tr' num2str(wichtr)];

        %----------- check save file
        save_eddyflux_dir = [save_dir '/C' num2str(wichtr,'%02d') '/'];
        if ~exist(save_eddyflux_dir,'dir'); mkdir(save_eddyflux_dir); end
        savefnm = [save_eddyflux_dir save_file_preStr yrstr '_' dystr '_' hrstr '.nc'];
        if exist(savefnm,'file')
            disp(['NC file already exists, skip : ', savefnm]);
            continue
        end
        
        %----------- read <uhc> and <c>
        uhc_fnm = [uhc_dir '/C' num2str(wichtr,'%02d') '/uhc__' yrstr '_' dystr '_' hrstr '.nc'];
        c = ncread(tr_fnm,varname);
        [uhc, vhc] = deal(ncread(uhc_fnm,'uhc'), ncread(uhc_fnm,'vhc'));
        
        %----------- <uhc> - <u><h><c>
        % <u><h><c>
        [umhmcm, vmhmcm] = calc_TFluxes_CG(u2d, v2d, h2d.*c(:,:,ik), '4th_order');
        % tracer eddy flux
        [eflxu, eflxv] = deal(uhc - umhmcm, vhc - vmhmcm);

        %----------- save
        dim_name = {'xu','yu','xv','yv'};
        dim_length = [size(eflxu,[1 2]), size(eflxv,[1 2])];
        savevarname = { save_varname_fu, save_varname_fv };
        data = {eflxu, eflxv};
        dimNum_of_var = {[1,2], [3,4]};
        %
        global_att  = [ 'Tracer eddy flx:<uhc>-<u><h><c>: ' ...
            'tracer=' varname '; uhc_fnm=' uhc_fnm '; ' ...
            'h_fnm=' h_fnm '; uv_fnm=' uv_fnm '; tr_fnm=' tr_fnm];
        FUN_nc_easywrite_enhanced( savefnm, dim_name, dim_length,...
            savevarname, dimNum_of_var, data, global_att )
        fprintf(1,'\nEddy fluxes saved to: %s...\n\n', savefnm);

    end
end