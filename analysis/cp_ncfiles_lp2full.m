clear

subexpStr_al = {'02','03','04','05','06','07','08','09','10'};
M = 2;

parfor (isub = 1:length(subexpStr_al), M)

    subexpStr = subexpStr_al{isub};
    dir_parent = ['/glade/campaign/univ/umia0037/lagr_study/exp4_' subexpStr '/'];

    %% mv # ptcls

    for lp = 1:6
        dir_aim = [dir_parent 'nptcls/full/Z01/'];
        if ~exist(dir_aim,'dir')
            mkdir(dir_aim);
        end

        dir_orig = [dir_parent 'nptcls/lp' num2str(lp,'%02d') '/Z01/'];
        fnms_orig = [dir_orig '*'];

        % cp
        fprintf(1,'Copying %s to %s\n',fnms_orig, dir_aim);
        comd = ['mv' 32 fnms_orig 32 dir_aim];  % 32 is a space
        system(comd);
    end

    %% cp c

    for lp = 1:6
        dir_aim = [dir_parent 'c/full/Z01/'];
        if ~exist(dir_aim,'dir')
            mkdir(dir_aim);
        end

        dir_orig = [dir_parent 'c/lp' num2str(lp,'%02d') '/Z01/'];
        fnms_orig = [dir_orig '*'];

        % cp
        fprintf(1,'Copying %s to %s\n',fnms_orig, dir_aim);
        comd = ['cp' 32 fnms_orig 32 dir_aim];  % 32 is a space
        system(comd);
    end

    %% cp h

    for lp = 1:6
        dir_aim = [dir_parent 'h/full/Z01/'];
        if ~exist(dir_aim,'dir')
            mkdir(dir_aim);
        end

        dir_orig = [dir_parent 'h/lp' num2str(lp,'%02d') '/Z01/'];
        fnms_orig = [dir_orig '*'];

        % cp
        fprintf(1,'Copying %s to %s\n',fnms_orig, dir_aim);
        comd = ['cp' 32 fnms_orig 32 dir_aim];  % 32 is a space
        system(comd);
    end


    %% cp uv
    for lp = 1:6
        dir_aim = [dir_parent 'uv/full/Z01/'];
        if ~exist(dir_aim,'dir')
            mkdir(dir_aim);
        end

        dir_orig = [dir_parent 'uv/lp' num2str(lp,'%02d') '/Z01/'];
        fnms_orig = [dir_orig '*'];

        % cp
        fprintf(1,'Copying %s to %s\n',fnms_orig, dir_aim);
        comd = ['cp' 32 fnms_orig 32 dir_aim];  % 32 is a space
        system(comd);
    end

    %% cp uhc
    for lp = 1:6
        for itr = 1:8
            dir_aim = [dir_parent 'uhc/full/Z01/C' num2str(itr,'%02d') '/'];
            if ~exist(dir_aim,'dir')
                mkdir(dir_aim);
            end

            dir_orig = [dir_parent 'uhc/lp' num2str(lp,'%02d') '/Z01/C' num2str(itr,'%02d') '/'];
            fnms_orig = [dir_orig '*'];

            % cp
            fprintf(1,'Copying %s to %s\n',fnms_orig, dir_aim);
            comd = ['cp' 32 fnms_orig 32 dir_aim];  % 32 is a space
            system(comd);
        end
    end

end % sub
delete(gcp('nocreate'))

