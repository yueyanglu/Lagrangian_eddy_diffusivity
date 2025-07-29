clear

%% cp uhc
dir_parent = '/glade/campaign/univ/umia0037/lagr_study/uhc';

for lp = 1:6
    for itr = 1:8
        dir_aim = [dir_parent '/full/Z01/C' num2str(itr,'%02d') '/'];
        if ~exist(dir_aim,'dir')
            mkdir(dir_aim);
        end

        dir_orig = [dir_parent '/lp' num2str(lp,'%02d') '/Z01/C' num2str(itr,'%02d') '/'];
        fnms_orig = [dir_orig '*'];

        % cp
        fprintf(1,'Copying %s to %s\n',fnms_orig, dir_aim);
        comd = ['cp' 32 fnms_orig 32 dir_aim];  % 32 is a space
        system(comd);
    end
end