% 
% Calc K-tensor from the eddy flux F={<uhc>}-{<u>}{<h>}{<c>} and {<h>}del{<c>}
%   where <> is coarse-graining and {} is time mean.
%  <uhc>, <u>, <h>, <c> are readed from files. In this script they will be
%   time averaged {}!
% 
%   Fields needed: <uhc>, <u>, <h>, <c> 
% 
clear
homedir = getenv('HOME');
workdir = getenv('WORK');
campdir = getenv('CAMP');
addpath(genpath([homedir '/work_MOM']));
addpath(genpath([homedir '/mytoolbox']));
addpath(genpath([homedir '/MyFuncs']));

%%

ik = 1;

%---- set 
carries = [1:8];
ntr = length(carries);
fprintf(1,'Using tracers: %s ...\n',mat2str(carries));

for subexpStrSh = 1:20

subexpStr = num2str(subexpStrSh,'%02d');

%-------- dir
exp_dir = [campdir '/lagr_study/exp10_' subexpStr '/'];  
uhc_dir = [exp_dir 'uhc/full/Z' num2str(ik,'%02d')];
uv_dir = [exp_dir 'uv/full/Z' num2str(ik,'%02d')];
h_dir = [exp_dir 'h/full/Z' num2str(ik,'%02d')];
tr_dir = [exp_dir 'c/full/Z' num2str(ik,'%02d')];

%-------- time average length
yr_s = 21;
[day_s, day_e, dt] = deal(1, 365*2, 1);
t_al = day_s:dt:day_e;
nt_al = length(t_al);

%-------- grid
cs_len = 32;
grid = build_grid_MOM(1024/cs_len,1024/cs_len,[0 3840],[0 3840]);

%--------- dir for saving
save_dir = [exp_dir '/ktensor/C' num2str(carries,'%02d') '/Z' num2str(ik,'%02d')];
 
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
fprintf(1,'K-tensor will be saved to: %s\n',save_dir);

%% read

u_sm = zeros(grid.niu,grid.nju);
v_sm = zeros(grid.niv,grid.njv);
h_sm = zeros(grid.nih,grid.njh);
c_sm = zeros(grid.nih,grid.njh,ntr);
[hcx_sm, uhc_sm, umhmcm] = deal(zeros(grid.niu,grid.nju,ntr));
[hcy_sm, vhc_sm, vmhmcm] = deal(zeros(grid.niv,grid.njv,ntr));

for it = 1:nt_al

    %----------------------- current time
    [yrstr, dystr, hrstr] = get_timestr(t_al(it), yr_s);
    fprintf(1,'\n time = Y%s-D%s-H%s (it=%d) of %d snapshots ...\n',...
        yrstr,dystr,hrstr,it,nt_al);

    %----------------------- read <u> <h> <c>
    uv_fnm = [uv_dir '/prog__' yrstr '_' dystr '_' hrstr '.nc'];
    fprintf(1,'\nReading uv from: %s...\n', uv_fnm);
    ds_uv = ncstruct(uv_fnm);
    [u, v] = deal(ds_uv.u, ds_uv.v);
    u(isnan(u)) = 0; v(isnan(v)) = 0;
    u_sm = u_sm + u; 
    v_sm = v_sm + v; 
    %
    h_fnm = [h_dir '/h_snap__' yrstr '_' dystr '_' hrstr '.nc'];
    fprintf(1,'\nReading h from: %s...\n', h_fnm);
    ds_h = ncstruct(h_fnm);
    h = ds_h.h;
    h(isnan(h)) = 0;
    h_sm = h_sm + h;
    %
    tr_fnm = [tr_dir '/tr__' yrstr '_' dystr '_' hrstr '.nc'];
    fprintf(1,'\nReading tr from: %s...\n\n', tr_fnm);
    ds_tr = ncstruct(tr_fnm);

    %----------------------- read all c and <uhc>
    fprintf(1,'\nReading uhc from: %s...\n\n', uhc_dir);
    for itr = 1:ntr
        wichtr = carries(itr);
        varname = ['tr' num2str(wichtr)];
        uhc_fnm = [uhc_dir '/C' num2str(wichtr,'%02d') '/uhc__' ...
            yrstr '_' dystr '_' hrstr '.nc'];
        ds_uhc = ncstruct(uhc_fnm);
        % convert nan as 0
        c = ds_tr.(varname);
        uhc = ds_uhc.uhc;
        vhc = ds_uhc.vhc;
        c(isnan(c)) = 0;
        uhc(isnan(uhc)) = 0;
        vhc(isnan(vhc)) = 0;

        %
        c_sm(:,:,itr) = c_sm(:,:,itr) + c; 
        uhc_sm(:,:,itr) = uhc_sm(:,:,itr) + uhc;
        vhc_sm(:,:,itr) = vhc_sm(:,:,itr) + vhc;
    end % ntr
end % nt

% -------- time mean [xyzC]
u_sm = u_sm / nt_al;
v_sm = v_sm / nt_al;
h_sm = h_sm / nt_al;
c_sm = c_sm / nt_al;
uhc_sm = uhc_sm / nt_al;
vhc_sm = vhc_sm / nt_al;
% eddy flux: {<uhc>}-{<u>}{<h>}{<c>} [m2/s*c] and h*delC [c]
for itr = 1:ntr
    [umhmcm(:,:,itr), vmhmcm(:,:,itr)] = calc_TFluxes_CG(u_sm, ...
        v_sm, h_sm.*c_sm(:,:,itr), '4th_order');
    [hcx_sm(:,:,itr),hcy_sm(:,:,itr)]  = calc_GxGy_CG(c_sm(:,:,itr),...
        h_sm,grid.dxCu,grid.dyCu,grid.dxCv,grid.dyCv,0);
end
hcx_sm([1 end],:,:) = 0;
hcy_sm(:,[1 end],:) = 0;
% 
edflxu = uhc_sm - umhmcm;
edflxv = vhc_sm - vmhmcm;

% -------- calc K [m2/s] 
fprintf('\nCalc K-tensor ...\n');
[Kxx,Kxy,Kyx,Kyy] = flx2K_MOM( hcx_sm, hcy_sm, edflxu, edflxv );
% Kxx = filter_extreme(Kxx,1,99); Kyx = filter_extreme(Kyx,1,99);
% Kxy = filter_extreme(Kxy,1,99); Kyy = filter_extreme(Kyy,1,99);

% -------- save --------
savename = [save_dir '/K__tmflx_D' num2str(day_s) '_' num2str(day_e) '.nc'];
dim_name = {'xu','yu','xv','yv', 'tr'};
dim_length = [grid.niu,grid.nju, grid.niv,grid.njv, ntr];
savevarname = { 'Kxx', 'Kxy', 'Kyx', 'Kyy', ...
    'hcx','hcy','fu','fv' };
data = {Kxx, Kxy, Kyx, Kyy, ...
    hcx_sm, hcy_sm, edflxu, edflxv};
dimNum_of_var = {[1,2], [1,2], [3,4], [3,4],...
    [1,2,5], [3,4,5], [1,2,5], [3,4,5]};
%
global_att  = [ 'K = -F/C from ptcls, F = {<uhc>}-{<u>}{<h>}{<c>}: ' ...
    'carries=' carries '; uhc_dir=' uhc_dir '; ' '; uv_dir=' uv_dir '; ' ...
    'h_fnm=' h_fnm '; tr_fnm=' tr_fnm '; day_s=' num2str(day_s) ...
    '; day_e=' num2str(day_e) '; dt=' num2str(dt)];
FUN_nc_easywrite_enhanced( savename, dim_name, dim_length,...
    savevarname, dimNum_of_var, data, global_att )
fprintf(1,'\nK saved to: %s...\n\n', savename);

end % subexpStr
