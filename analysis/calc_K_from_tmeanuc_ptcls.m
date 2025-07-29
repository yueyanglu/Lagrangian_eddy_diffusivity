% 
% Calc K-tensor from the eddy flux F={<uc>}-{<u>}{<c>} and del{<c>}
%   where <> is coarse-graining and {} is time mean.
%  <uc>, <u>, <c> are readed from files. In this script they will be
%   time averaged {}!
% 
%   Fields needed: <uc>, <u>, <c> 
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

for subexpStrSh = 1:10

subexpStr = num2str(subexpStrSh,'%02d');

%-------- dir
exp_dir = [campdir '/lagr_study/reg1_exp4_' subexpStr '/'];  
%exp_dir = [campdir '/lagr_study/exp1/'];  
%exp_dir = [campdir '/lagr_study/reg1_exp1_0/'];  
uc_dir = [exp_dir 'uc/full/Z' num2str(ik,'%02d')];
uv_dir = [exp_dir 'uv/full/Z' num2str(ik,'%02d')];
tr_dir = [exp_dir 'c/full/Z' num2str(ik,'%02d')];
%tr_dir = [exp_dir 'c_interp/full/Z' num2str(ik,'%02d')];

%-------- time average length
yr_s = 21;
[day_s, day_e, dt] = deal(1, 365*2, 1);
t_al = day_s:dt:day_e;
nt_al = length(t_al);

%-------- grid
cs_len = 32;
grid = build_grid_MOM(1024/cs_len,1024/cs_len,[0 3840],[0 3840]);

%--------- dir for saving
save_dir = [exp_dir '/ktensor_noh/C' num2str(carries,'%02d') '/Z' num2str(ik,'%02d')];
 
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
fprintf(1,'K-tensor will be saved to: %s\n',save_dir);

%% read

u_sm = zeros(grid.niu,grid.nju);
v_sm = zeros(grid.niv,grid.njv);
c_sm = zeros(grid.nih,grid.njh,ntr);
[cx_sm, uc_sm, umcm] = deal(zeros(grid.niu,grid.nju,ntr));
[cy_sm, vc_sm, vmcm] = deal(zeros(grid.niv,grid.njv,ntr));

% count nan values
ctnan_u = 0*u_sm;
ctnan_v = 0*v_sm;
ctnan_c = 0*c_sm;
ctnan_uc = 0*uc_sm;
ctnan_vc = 0*vc_sm;


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
    u(isnan(u)) = 0; ctnan_u = ctnan_u + isnan(u);
    v(isnan(v)) = 0; ctnan_v = ctnan_v + isnan(v);
    u_sm = u_sm + u; 
    v_sm = v_sm + v; 
    %
    tr_fnm = [tr_dir '/tr__' yrstr '_' dystr '_' hrstr '.nc'];
    fprintf(1,'\nReading tr from: %s...\n\n', tr_fnm);
    ds_tr = ncstruct(tr_fnm);

    %----------------------- read all c and <uc>
    fprintf(1,'\nReading uc from: %s...\n\n', uc_dir);
    for itr = 1:ntr
        wichtr = carries(itr);
        varname = ['tr' num2str(wichtr)];
        uc_fnm = [uc_dir '/C' num2str(wichtr,'%02d') '/uc__' ...
            yrstr '_' dystr '_' hrstr '.nc'];
        ds_uc = ncstruct(uc_fnm);
        c = ds_tr.(varname);
        uc = ds_uc.uc;
        vc = ds_uc.vc;
        % convert nan to 0
        c(isnan(c)) = 0; ctnan_c = ctnan_c + isnan(c);
        uc(isnan(uc)) = 0; ctnan_uc = ctnan_uc + isnan(uc);
        vc(isnan(vc)) = 0; ctnan_vc = ctnan_vc + isnan(vc);

        %
        c_sm(:,:,itr) = c_sm(:,:,itr) + c; 
        uc_sm(:,:,itr) = uc_sm(:,:,itr) + uc;
        vc_sm(:,:,itr) = vc_sm(:,:,itr) + vc;
    end % ntr
end % nt

% -------- time mean [xyzC]
u_sm = u_sm ./ (nt_al - ctnan_u);
v_sm = v_sm ./ (nt_al - ctnan_v);
c_sm = c_sm ./ (nt_al - ctnan_c);
uc_sm = uc_sm ./ (nt_al - ctnan_uc);
vc_sm = vc_sm ./ (nt_al - ctnan_vc);
% eddy flux: {<uc>}-{<u>}{<c>} [m/s*c] and delC [c/m]
for itr = 1:ntr
    [umcm(:,:,itr), vmcm(:,:,itr)] = calc_TFluxes_CG(u_sm, ...
        v_sm, c_sm(:,:,itr), '4th_order');
    [cx_sm(:,:,itr),cy_sm(:,:,itr)]  = calc_GxGy_CG(c_sm(:,:,itr),...
        ones(grid.nih,grid.njh),grid.dxCu,grid.dyCu,grid.dxCv,grid.dyCv,0);
end
cx_sm([1 end],:,:) = 0;
cy_sm(:,[1 end],:) = 0;
% 
edflxu = uc_sm - umcm;
edflxv = vc_sm - vmcm;

% 
cx_sm(cx_sm==0) = NaN;
cy_sm(cy_sm==0) = NaN;
edflxu(edflxu==0) = NaN;
edflxv(edflxv==0) = NaN;

% -------- calc K [m2/s] 
fprintf('\nCalc K-tensor ...\n');
[Kxx,Kxy,Kyx,Kyy] = flx2K_MOM( cx_sm, cy_sm, edflxu, edflxv );
% Kxx = filter_extreme(Kxx,1,99); Kyx = filter_extreme(Kyx,1,99);
% Kxy = filter_extreme(Kxy,1,99); Kyy = filter_extreme(Kyy,1,99);

% -------- save --------
savename = [save_dir '/K__tmflx_D' num2str(day_s) '_' num2str(day_e) '.nc'];
dim_name = {'xu','yu','xv','yv', 'tr'};
dim_length = [grid.niu,grid.nju, grid.niv,grid.njv, ntr];
savevarname = { 'Kxx', 'Kxy', 'Kyx', 'Kyy', ...
    'cx','cy','fu','fv' };
data = {Kxx, Kxy, Kyx, Kyy, ...
    cx_sm, cy_sm, edflxu, edflxv};
dimNum_of_var = {[1,2], [1,2], [3,4], [3,4],...
    [1,2,5], [3,4,5], [1,2,5], [3,4,5]};
%
global_att  = [ 'K = -F/C from ptcls, F = {<uc>}-{<u>}{<c>}: ' ...
    'carries=' carries '; uc_dir=' uc_dir '; ' '; uv_dir=' uv_dir '; ' ...
    '; tr_fnm=' tr_fnm '; day_s=' num2str(day_s) ...
    '; day_e=' num2str(day_e) '; dt=' num2str(dt)];
FUN_nc_easywrite_enhanced( savename, dim_name, dim_length,...
    savevarname, dimNum_of_var, data, global_att )
fprintf(1,'\nK saved to: %s...\n\n', savename);

end % subexpStr
