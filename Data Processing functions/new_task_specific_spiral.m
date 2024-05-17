function [subject] = new_task_specific_spiral(subject)

restoredefaultpath
cd '/headnode2/akar5239/BrainVortexToolbox-main'
addpath(genpath([pwd]))
main_folder = pwd;


No_of_Subject = 1; % number of subjects used for analysis, randomly selected from HCP database (S1200)


hemisphere = 1; % 1 for left hemisphere, 2 for right hemisphere, keep this as 1 for demonstration purpose

% flagTask:
% 1 = language task, original 100 subjects;
% 2 = language task, additional 100 subjects;
% 3 = working memory task;
% 4 = Motor Task
flagTask = 4;




disp(['loading task label...'])
if flagTask == 4 % Motor Task
    % load task label of each subject
    foldername = [main_folder,'/Sample Data/Motor Task/TaskLabel'];
    cd(foldername)
    name = dir(pwd) ;
    file_name2 = ['MotorTaskLabelAllSubject.mat'];
    load (file_name2);
    for isubject = 1:150
        fullTime_allsubject{isubject} = TaskLabel_AllSubject_motor{isubject};
    end
end



%% Motor Task
if flagTask == 4 % Motor Task
    % define parameters
    session_duration = 17;   % (currently this is the number of frames 12 seconds/0.72 )


    disp('[Loading Raw Phase Signals...]')
    tic
    if hemisphere == 1
        % load spatiotemporally bandpass filtered (smoothed) real data
        folder_name = [main_folder,'/Sample Data/Motor Task/Preprocessed Data'];
        cd(folder_name)
        filename = ['Preprocessed_spatiotemporalbandpass_data_motor_task_LEFT_sub',num2str(subject),'.mat'];
        load(filename)

        DataIn_smooth = DataOut(:,:,1:283);
    elseif hemisphere == 2
        % load spatiotemporally bandpass filtered (smoothed) real data
        folder_name = [main_folder,'/Sample Data/Motor Task/Preprocessed Data'];
        cd(folder_name)
        filename = ['Preprocessed_spatiotemporalbandpass_data_motor_task_RIGHT_sub',num2str(subject),'.mat'];
        load(filename)

        DataIn_smooth = DataOut(:,:,1:283);
    end
    toc

    %  Find phase map of fMRI signal
    smooth_phase_map_real = nan(size(DataIn_smooth(1:175,:,:)));
    for irow = 1:size(smooth_phase_map_real,1)
        for icol = 1:size(smooth_phase_map_real,2)
            temp1 = DataIn_smooth(irow,icol,:);
            if nansum(temp1(:))~=0
                smooth_phase_map_real(irow,icol,:) = angle(hilbert(temp1(:)));
            end
        end
    end

    disp(['extracting task specific spiral centre distributions...'])
    tic
    foldername = [main_folder,'/Sample Data/Motor Task/Spiral Detected'];
    cd(foldername)
    if hemisphere == 1
        filename = ['Spiral_detected_surfilt_motor_task_LEFT_sub',num2str(subject),'.mat'];
        load(filename)
    elseif hemisphere == 2
        filename = ['Spiral_detected_surfilt_motor_task_RIGHT_sub',num2str(subject),'.mat'];
        load(filename)
    end
    toc


    spiral_filt_nega = spiral_filt_nega_real_95perc_extend;      
    spiral_filt_pos = spiral_filt_pos_real_95perc_extend;

    % extract task label for each subject
    temp1 = fullTime_allsubject{subject};
    temp1_time = temp1(:,3);

    % find time points for "Right" Tasks (even for hand odd for foot)
    count = find(temp1_time==4); % 4 for Right Hand Tasks
    start_end_time_right_hand = temp1(count,:); 

    count = find(temp1_time==3); % 3 for Right Foot Tasks
    start_end_time_right_foot = temp1(count,:);

    count = find(temp1_time==5); % 5 for tongue 
    start_end_time_tongue = temp1(count,:);

    % time points for "Left" tasks (even for hand odd for foot)
    count = find(temp1_time==1); % 1 for Left Foot Tasks
    start_end_time_left_foot = temp1(count,:);

    count = find(temp1_time==2); % 2 for Left Hand Tasks
    start_end_time_left_hand = temp1(count,:);


    % extract task specifc spiral distribution data: "Right" tasks

    % Movement periods for the right hand
    if size(start_end_time_right_hand,1) == 2
        block_1 = start_end_time_right_hand(1,1):start_end_time_right_hand(1,1)+session_duration;
        block_2 = start_end_time_right_hand(2,1):start_end_time_right_hand(2,1)+session_duration;
     

        [generated_vortices, modulated_vortices] = task_onset_id(spiral_filt_pos,spiral_filt_pos_centreONLY_real, [block_1; block_2]);
        
        %identify task relevant pos spiral filters 

        rh_block1_pos_generated_vortex = generated_vortices{1};
        rh_block2_pos_generated_vortex = generated_vortices{2};

        rh_block1_pos_modulated_vortex = modulated_vortices{1};
        rh_block2_pos_modulated_vortex = modulated_vortices{2};


        [generated_vortices, modulated_vortices] = task_onset_id(spiral_filt_nega,spiral_filt_nega_centreONLY_real, [block_1; block_2]);

        
        %identify task relevant neg spiral filters 


        rh_block1_nega_generated_vortex = generated_vortices{1};
        rh_block2_nega_generated_vortex = generated_vortices{2};

        rh_block1_nega_modulated_vortex = modulated_vortices{1};
        rh_block2_nega_modulated_vortex = modulated_vortices{2};
    end

    % Movement periods for the right foot
    if size(start_end_time_right_foot,1) == 2
        block_1 = start_end_time_right_foot(1,1):start_end_time_right_foot(1,1)+session_duration;
        block_2 = start_end_time_right_foot(2,1):start_end_time_right_foot(2,1)+session_duration;

        [generated_vortices, modulated_vortices] = task_onset_id(spiral_filt_pos, spiral_filt_pos_centreONLY_real, [block_1; block_2]);


        rf_block1_pos_generated_vortex = generated_vortices{1};
        rf_block2_pos_generated_vortex = generated_vortices{2};

        rf_block1_pos_modulated_vortex = modulated_vortices{1};
        rf_block2_pos_modulated_vortex = modulated_vortices{2};

        [generated_vortices, modulated_vortices] = task_onset_id(spiral_filt_nega,spiral_filt_nega_centreONLY_real, [block_1; block_2]);


        rf_block1_nega_generated_vortex = generated_vortices{1};
        rf_block2_nega_generated_vortex = generated_vortices{2};

        rf_block1_nega_modulated_vortex = modulated_vortices{1};
        rf_block2_nega_modulated_vortex = modulated_vortices{2};

    end

    % Movement periods for the tongue
    if size(start_end_time_tongue,1) == 2
        block_1 = start_end_time_tongue(1,1):start_end_time_tongue(1,1)+session_duration;
        block_2 = start_end_time_tongue(2,1):start_end_time_tongue(2,1)+session_duration;

        [generated_vortices, modulated_vortices] = task_onset_id(spiral_filt_pos,spiral_filt_pos_centreONLY_real , [block_1; block_2]);

        t_block1_pos_generated_vortex = generated_vortices{1};
        t_block2_pos_generated_vortex = generated_vortices{2};

        t_block1_pos_modulated_vortex = modulated_vortices{1};
        t_block2_pos_modulated_vortex = modulated_vortices{2};

        [generated_vortices, modulated_vortices] = task_onset_id(spiral_filt_nega, spiral_filt_nega_centreONLY_real, [block_1; block_2]);

        t_block1_nega_generated_vortex = generated_vortices{1};
        t_block2_nega_generated_vortex = generated_vortices{2};

        t_block1_nega_modulated_vortex = modulated_vortices{1};
        t_block2_nega_modulated_vortex = modulated_vortices{2};
    end

    % extract task specifc spiral distribution data: "Left" tasks


    % Movement Periods for the Left Foot
    if size(start_end_time_left_foot,1) == 2
        block_1 = start_end_time_left_foot(1,1):start_end_time_left_foot(1,1)+session_duration;
        block_2 = start_end_time_left_foot(2,1):start_end_time_left_foot(2,1)+session_duration;

        [generated_vortices, modulated_vortices] = task_onset_id(spiral_filt_pos, spiral_filt_pos_centreONLY_real, [block_1; block_2]);

        lf_block1_pos_generated_vortex = generated_vortices{1};
        lf_block2_pos_generated_vortex = generated_vortices{2};

        lf_block1_pos_modulated_vortex = modulated_vortices{1};
        lf_block2_pos_modulated_vortex = modulated_vortices{2};

        [generated_vortices, modulated_vortices] = task_onset_id(spiral_filt_nega,spiral_filt_nega_centreONLY_real , [block_1; block_2]);


        lf_block1_nega_generated_vortex = generated_vortices{1};
        lf_block2_nega_generated_vortex = generated_vortices{2};

        lf_block1_nega_modulated_vortex = modulated_vortices{1};
        lf_block2_nega_modulated_vortex = modulated_vortices{2};

    end

    % Movement Periods for the Left Hand
    if size(start_end_time_left_hand,1) == 2
        block_1 = start_end_time_left_hand(1,1):start_end_time_left_hand(1,1)+session_duration;
        block_2 = start_end_time_left_hand(2,1):start_end_time_left_hand(2,1)+session_duration;

        [generated_vortices, modulated_vortices] = task_onset_id(spiral_filt_pos, spiral_filt_pos_centreONLY_real, [block_1; block_2]);


        lh_block1_pos_generated_vortex = generated_vortices{1};
        lh_block2_pos_generated_vortex = generated_vortices{2};

        lh_block1_pos_modulated_vortex = modulated_vortices{1};
        lh_block2_pos_modulated_vortex = modulated_vortices{2};

        [generated_vortices, modulated_vortices] = task_onset_id(spiral_filt_nega,spiral_filt_nega_centreONLY_real, [block_1; block_2]);

        lh_block1_nega_generated_vortex = generated_vortices{1};
        lh_block2_nega_generated_vortex = generated_vortices{2};

        lh_block1_nega_modulated_vortex = modulated_vortices{1};
        lh_block2_nega_modulated_vortex = modulated_vortices{2};

    end
    



    % FIRST movement blocks, a distinction is made between earlier and
    % later recording periods     
    
    block1_rightHand_pos_neg{1,1} = rh_block1_pos_modulated_vortex;
  
    block1_rightHand_pos_neg{1,2} = rh_block1_pos_generated_vortex;

    block1_rightHand_pos_neg{1,3} = rh_block1_nega_modulated_vortex;

    block1_rightHand_pos_neg{1,4} = rh_block1_nega_generated_vortex;


    
    block1_leftHand_pos_neg{1,1} = lh_block1_pos_modulated_vortex;

    block1_leftHand_pos_neg{1,2} = lh_block1_pos_generated_vortex;

    block1_leftHand_pos_neg{1,3} = lh_block1_nega_modulated_vortex;

    block1_leftHand_pos_neg{1,4} = lh_block1_nega_generated_vortex;



    block1_rightFoot_pos_neg{1,1} = rf_block1_pos_modulated_vortex;

    block1_rightFoot_pos_neg{1,2} = rf_block1_pos_generated_vortex;

    block1_rightFoot_pos_neg{1,3} = rf_block1_nega_modulated_vortex;

    block1_rightFoot_pos_neg{1,4} = rf_block1_nega_generated_vortex;



    block1_leftFoot_pos_neg{1,1} = lf_block1_pos_modulated_vortex;

    block1_leftFoot_pos_neg{1,2} = lf_block1_pos_generated_vortex;

    block1_leftFoot_pos_neg{1,3} = lf_block1_nega_modulated_vortex;

    block1_leftFoot_pos_neg{1,4} = lf_block1_nega_generated_vortex;



    block1_tongue_pos_neg{1,1} = t_block1_pos_modulated_vortex;

    block1_tongue_pos_neg{1,2} = t_block1_pos_generated_vortex;

    block1_tongue_pos_neg{1,3} = t_block1_nega_modulated_vortex;

    block1_tongue_pos_neg{1,4} = t_block1_nega_generated_vortex;

    
    % SECOND Movement block

    block2_rightHand_pos_neg{1,1} = rh_block2_pos_modulated_vortex;

    block2_rightHand_pos_neg{1,2} = rh_block2_pos_generated_vortex;

    block2_rightHand_pos_neg{1,3} = rh_block2_nega_modulated_vortex;

    block2_rightHand_pos_neg{1,4} = rh_block2_nega_generated_vortex;


    
    block2_leftHand_pos_neg{1,1} = lh_block2_pos_modulated_vortex;

    block2_leftHand_pos_neg{1,2} = lh_block2_pos_generated_vortex;

    block2_leftHand_pos_neg{1,3} = lh_block2_nega_modulated_vortex;

    block2_leftHand_pos_neg{1,4} = lh_block2_nega_generated_vortex;



    block2_rightFoot_pos_neg{1,1} = rf_block2_pos_modulated_vortex;

    block2_rightFoot_pos_neg{1,2} = rf_block2_pos_generated_vortex;

    block2_rightFoot_pos_neg{1,3} = rf_block2_nega_modulated_vortex;

    block2_rightFoot_pos_neg{1,4} = rf_block2_nega_generated_vortex;



    block2_leftFoot_pos_neg{1,1} = lf_block2_pos_modulated_vortex;

    block2_leftFoot_pos_neg{1,2} = lf_block2_pos_generated_vortex;

    block2_leftFoot_pos_neg{1,3} = lf_block2_nega_modulated_vortex;
 
    block2_leftFoot_pos_neg{1,4} = lf_block2_nega_generated_vortex;



    block2_tongue_pos_neg{1,1} = t_block2_pos_modulated_vortex;

    block2_tongue_pos_neg{1,2} = t_block2_pos_generated_vortex;

    block2_tongue_pos_neg{1,3} = t_block2_nega_modulated_vortex;

    block2_tongue_pos_neg{1,4} = t_block2_nega_generated_vortex;
       


        foldername = [main_folder,'/Sample Data/Motor Task/Analysis/'];
        if hemisphere == 1
        filename = ['task_specific_spiral_ID_leftHem_sub',num2str(subject),'.mat'];
        elseif hemisphere == 2
            filename = ['task_specific_spiral_ID_rightHem_sub',num2str(subject),'.mat'];
        end
        save([foldername,filename],'block1_rightHand_pos_neg','block1_leftHand_pos_neg','block1_rightFoot_pos_neg','block1_leftFoot_pos_neg','block1_tongue_pos_neg',...
                                   'block2_rightHand_pos_neg','block2_leftHand_pos_neg','block2_rightFoot_pos_neg','block2_leftFoot_pos_neg','block2_tongue_pos_neg');

    cd(main_folder)



 
%% Spiral Trial averaging May be relevant as part of preliminary results
% 
% Right Finger Movement: Subject Averaged spiral distribution
    RH_block1_sub_dist = struct();
    RH_block2_sub_dist = struct();

%--> Modulated Spirals (exist prior to task onset and persist for more than 3 frames)
    RH_block1_sub_dist.modulated.pos = zeros(175,251,283);
    RH_block1_sub_dist.modulated.neg = zeros(175,251,283);

    RH_block1_sub_dist.generated.pos = zeros(175,251,283);
    RH_block1_sub_dist.generated.neg = zeros(175,251,283);
    
    
    for ipatt_index = 1:size(rh_block1_pos_modulated_vortex,1)
        if ~isempty(rh_block1_pos_modulated_vortex{1})
            ipatt = rh_block1_pos_modulated_vortex{ipatt_index,1};
            for time = rh_block1_pos_modulated_vortex{ipatt_index,3}
                temp1 = full(spiral_filt_pos{ipatt,time});
                RH_block1_sub_dist.modulated.pos(:,:,time) =  RH_block1_sub_dist.modulated.pos(:,:,time) + temp1;
            end
        end
    end


    for ipatt_index = 1:size(rh_block1_nega_modulated_vortex,1)
        if ~isempty(rh_block1_nega_modulated_vortex{1})
        ipatt = rh_block1_nega_modulated_vortex{ipatt_index,1};
        for time = rh_block1_nega_modulated_vortex{ipatt_index,3}
            temp1 = full(spiral_filt_nega{ipatt,time});
            RH_block1_sub_dist.modulated.neg(:,:,time) = RH_block1_sub_dist.modulated.neg(:,:,time) + temp1;
        end
        end
    end

%--> Generated Spirals (Appear anytime within the task period)
    for ipatt_index = 1:size(rh_block1_pos_generated_vortex,1)
        if ~isempty(rh_block1_pos_generated_vortex{1})
        ipatt = rh_block1_pos_generated_vortex{ipatt_index,1};
        for time = rh_block1_pos_generated_vortex{ipatt_index,3}
            temp1 = full(spiral_filt_pos{ipatt,time});
             RH_block1_sub_dist.generated.pos(:,:,time) =  RH_block1_sub_dist.generated.pos(:,:,time) + temp1;
        end
        end
    end


    for ipatt_index = 1:size(rh_block1_nega_generated_vortex,1)
        if ~isempty(rh_block1_nega_generated_vortex{1})
        ipatt = rh_block1_nega_generated_vortex{ipatt_index,1};
        for time = rh_block1_nega_generated_vortex{ipatt_index,3}
            temp1 = full(spiral_filt_nega{ipatt,time});
            RH_block1_sub_dist.generated.neg(:,:,time) = RH_block1_sub_dist.generated.neg(:,:,time) + temp1;
        end
        end
    end


      % BLOCK 2

 %--> Modulated Spirals (exist prior to task onset and persist for more than 3 frames)
    RH_block2_sub_dist.modulated.pos = zeros(175,251,283);
    RH_block2_sub_dist.modulated.neg = zeros(175,251,283);

    RH_block2_sub_dist.generated.pos = zeros(175,251,283);
    RH_block2_sub_dist.generated.neg = zeros(175,251,283);


    for ipatt_index = 1:size(rh_block2_pos_modulated_vortex,1)
        if ~isempty(rh_block2_pos_modulated_vortex{1})
        ipatt = rh_block2_pos_modulated_vortex{ipatt_index,1};
        for time = rh_block2_pos_modulated_vortex{ipatt_index,3}
            temp1 = full(spiral_filt_pos{ipatt,time});
            RH_block2_sub_dist.modulated.pos(:,:,time) =  RH_block2_sub_dist.modulated.pos(:,:,time) + temp1;
        end
        end
    end


    for ipatt_index = 1:size(rh_block2_nega_modulated_vortex,1)
        if ~isempty(rh_block2_nega_modulated_vortex{1})
        ipatt = rh_block2_nega_modulated_vortex{ipatt_index,1};
        for time = rh_block2_nega_modulated_vortex{ipatt_index,3}
            temp1 = full(spiral_filt_nega{ipatt,time});
            RH_block2_sub_dist.modulated.neg(:,:,time) = RH_block2_sub_dist.modulated.neg(:,:,time) + temp1;
        end
        end
    end

  
%--> Generated Spirals (Appear anytime within the task period)
for ipatt_index = 1:size(rh_block2_pos_generated_vortex,1)
    if ~isempty(rh_block2_pos_generated_vortex{1})
        ipatt = rh_block2_pos_generated_vortex{ipatt_index,1};
        for time = rh_block2_pos_generated_vortex{ipatt_index,3}
            temp1 = full(spiral_filt_pos{ipatt,time});
            RH_block2_sub_dist.generated.pos(:,:,time) =  RH_block2_sub_dist.generated.pos(:,:,time) + temp1;
        end
    end
end


for ipatt_index = 1:size(rh_block2_nega_generated_vortex,1)
    if ~isempty(rh_block2_nega_generated_vortex{1})
        ipatt = rh_block2_nega_generated_vortex{ipatt_index,1};
        for time = rh_block2_nega_generated_vortex{ipatt_index,3}
            temp1 = full(spiral_filt_nega{ipatt,time});
            RH_block2_sub_dist.generated.neg(:,:,time) = RH_block2_sub_dist.generated.neg(:,:,time) + temp1;
        end
    end
end

    %% LEFT Finger Movement: Subject Averaged spiral distribution
    LH_block1_sub_dist = struct();
    LH_block2_sub_dist = struct();

%--> Modulated Spirals (exist prior to task onset and persist for more than 3 frames)
    LH_block1_sub_dist.modulated.pos = zeros(175,251,283);
    LH_block1_sub_dist.modulated.neg = zeros(175,251,283);

    LH_block1_sub_dist.generated.pos = zeros(175,251,283);
    LH_block1_sub_dist.generated.neg = zeros(175,251,283);

    
    for ipatt_index = 1:size(lh_block1_pos_modulated_vortex,1)
        if ~isempty(lh_block1_pos_modulated_vortex{1})
            ipatt = lh_block1_pos_modulated_vortex{ipatt_index,1};
            for time = lh_block1_pos_modulated_vortex{ipatt_index,3}
                temp1 = full(spiral_filt_pos{ipatt,time});
                LH_block1_sub_dist.modulated.pos(:,:,time) =  LH_block1_sub_dist.modulated.pos(:,:,time) + temp1;
            end
        end
    end

    
    for ipatt_index = 1:size(lh_block1_nega_modulated_vortex,1)
        if ~isempty(lh_block1_nega_modulated_vortex{1})
            ipatt = lh_block1_nega_modulated_vortex{ipatt_index,1};
            for time = lh_block1_nega_modulated_vortex{ipatt_index,3}
                temp1 = full(spiral_filt_nega{ipatt,time});
                LH_block1_sub_dist.modulated.neg(:,:,time) = LH_block1_sub_dist.modulated.neg(:,:,time) + temp1;
            end
        end
    end
    
    %--> Generated Spirals (Appear anytime within the task period)
    for ipatt_index = 1:size(lh_block1_pos_generated_vortex,1)
        if ~isempty(lh_block1_pos_generated_vortex{1})
            ipatt = lh_block1_pos_generated_vortex{ipatt_index,1};
            for time = lh_block1_pos_generated_vortex{ipatt_index,3}
                temp1 = full(spiral_filt_pos{ipatt,time});
                LH_block1_sub_dist.generated.pos(:,:,time) =  LH_block1_sub_dist.generated.pos(:,:,time) + temp1;
            end
        end
    end
    
    
    for ipatt_index = 1:size(lh_block1_nega_generated_vortex,1)
        if ~isempty(lh_block1_nega_generated_vortex{1})
            ipatt = lh_block1_nega_generated_vortex{ipatt_index,1};
            for time = lh_block1_nega_generated_vortex{ipatt_index,3}
                temp1 = full(spiral_filt_nega{ipatt,time});
                LH_block1_sub_dist.generated.neg(:,:,time) = LH_block1_sub_dist.generated.neg(:,:,time) + temp1;
            end
        end
    end


      % BLOCK 2

 %--> Modulated Spirals (exist prior to task onset and persist for more than 3 frames)
    LH_block2_sub_dist.modulated.pos = zeros(175,251,283);
    LH_block2_sub_dist.modulated.neg = zeros(175,251,283);

    LH_block2_sub_dist.generated.pos = zeros(175,251,283);
    LH_block2_sub_dist.generated.neg = zeros(175,251,283);


    for ipatt_index = 1:size(lh_block2_pos_modulated_vortex,1)
        if ~isempty(lh_block2_pos_modulated_vortex{1})
            ipatt = lh_block2_pos_modulated_vortex{ipatt_index,1};
            for time = lh_block2_pos_modulated_vortex{ipatt_index,3}
                temp1 = full(spiral_filt_pos{ipatt,time});
                LH_block2_sub_dist.modulated.pos(:,:,time) =  LH_block2_sub_dist.modulated.pos(:,:,time) + temp1;
            end
        end
    end


    for ipatt_index = 1:size(lh_block2_nega_modulated_vortex,1)
        if ~isempty(lh_block2_nega_modulated_vortex{1})
            ipatt = lh_block2_nega_modulated_vortex{ipatt_index,1};
            for time = lh_block2_nega_modulated_vortex{ipatt_index,3}
                temp1 = full(spiral_filt_nega{ipatt,time});
                LH_block2_sub_dist.modulated.neg(:,:,time) = LH_block2_sub_dist.modulated.neg(:,:,time) + temp1;
            end
        end
    end


    %--> Generated Spirals (Appear anytime within the task period)
    for ipatt_index = 1:size(lh_block2_pos_generated_vortex,1)
        if ~isempty(lh_block2_pos_generated_vortex{1})
            ipatt = lh_block2_pos_generated_vortex{ipatt_index,1};
            for time = lh_block2_pos_generated_vortex{ipatt_index,3}
                temp1 = full(spiral_filt_pos{ipatt,time});
                LH_block2_sub_dist.generated.pos(:,:,time) =  LH_block2_sub_dist.generated.pos(:,:,time) + temp1;
            end
        end
    end


    for ipatt_index = 1:size(lh_block2_nega_generated_vortex,1)
        if ~isempty(lh_block2_nega_generated_vortex{1})
            ipatt = lh_block2_nega_generated_vortex{ipatt_index,1};
            for time = lh_block2_nega_generated_vortex{ipatt_index,3}
                temp1 = full(spiral_filt_nega{ipatt,time});
                LH_block2_sub_dist.generated.neg(:,:,time) = LH_block2_sub_dist.generated.neg(:,:,time) + temp1;
            end
        end
    end
%% % RIGHT toe Movement: Subject Averaged spiral distribution
    RF_block1_sub_dist = struct();
    RF_block2_sub_dist = struct();

%--> Modulated Spirals (exist prior to task onset and persist for more than 3 frames)
    RF_block1_sub_dist.modulated.pos = zeros(175,251,283);
    RF_block1_sub_dist.modulated.neg = zeros(175,251,283);

    RF_block1_sub_dist.generated.pos = zeros(175,251,283);
    RF_block1_sub_dist.generated.neg = zeros(175,251,283);

    
    for ipatt_index = 1:size(rf_block1_pos_modulated_vortex,1)
        if ~isempty(rf_block1_pos_modulated_vortex{1})
            ipatt = rf_block1_pos_modulated_vortex{ipatt_index,1};
            for time = rf_block1_pos_modulated_vortex{ipatt_index,3}
                temp1 = full(spiral_filt_pos{ipatt,time});
                RF_block1_sub_dist.modulated.pos(:,:,time) =  RF_block1_sub_dist.modulated.pos(:,:,time) + temp1;
            end
        end
    end

 
    for ipatt_index = 1:size(rf_block1_nega_modulated_vortex,1)
        if ~isempty(rf_block1_nega_modulated_vortex{1})
            ipatt = rf_block1_nega_modulated_vortex{ipatt_index,1};
            for time = rf_block1_nega_modulated_vortex{ipatt_index,3}
                temp1 = full(spiral_filt_nega{ipatt,time});
                RF_block1_sub_dist.modulated.neg(:,:,time) = RF_block1_sub_dist.modulated.neg(:,:,time) + temp1;
            end
        end
    end
    %--> Generated Spirals (Appear anytime within the task period)
    for ipatt_index = 1:size(rf_block1_pos_generated_vortex,1)
        if ~isempty(rf_block1_pos_generated_vortex{1})
            ipatt = rf_block1_pos_generated_vortex{ipatt_index,1};
            for time = rf_block1_pos_generated_vortex{ipatt_index,3}
                temp1 = full(spiral_filt_pos{ipatt,time});
                RF_block1_sub_dist.generated.pos(:,:,time) =  RF_block1_sub_dist.generated.pos(:,:,time) + temp1;
            end
        end
    end


    for ipatt_index = 1:size(rf_block1_nega_generated_vortex,1)
        if ~isempty(rf_block1_nega_generated_vortex{1})
            ipatt = rf_block1_nega_generated_vortex{ipatt_index,1};
            for time = rf_block1_nega_generated_vortex{ipatt_index,3}
                temp1 = full(spiral_filt_nega{ipatt,time});
                RF_block1_sub_dist.generated.neg(:,:,time) = RF_block1_sub_dist.generated.neg(:,:,time) + temp1;
            end
        end
    end


      % BLOCK 2

 %--> Modulated Spirals (exist prior to task onset and persist for more than 3 frames)
    RF_block2_sub_dist.modulated.pos = zeros(175,251,283);
    RF_block2_sub_dist.modulated.neg = zeros(175,251,283);

    RF_block2_sub_dist.generated.pos = zeros(175,251,283);
    RF_block2_sub_dist.generated.neg = zeros(175,251,283);


    for ipatt_index = 1:size(rf_block2_pos_modulated_vortex,1)
        if ~isempty(rf_block2_pos_modulated_vortex{1})
            ipatt = rf_block2_pos_modulated_vortex{ipatt_index,1};
            for time = rf_block2_pos_modulated_vortex{ipatt_index,3}
                temp1 = full(spiral_filt_pos{ipatt,time});
                RF_block2_sub_dist.modulated.pos(:,:,time) =  RF_block2_sub_dist.modulated.pos(:,:,time) + temp1;
            end
        end
    end


    for ipatt_index = 1:size(rf_block2_nega_modulated_vortex,1)
        if ~isempty(rf_block2_nega_modulated_vortex{1})
            ipatt = rf_block2_nega_modulated_vortex{ipatt_index,1};
            for time = rf_block2_nega_modulated_vortex{ipatt_index,3}
                temp1 = full(spiral_filt_nega{ipatt,time});
                RF_block2_sub_dist.modulated.neg(:,:,time) = RF_block2_sub_dist.modulated.neg(:,:,time) + temp1;
            end
        end
    end


    %--> Generated Spirals (Appear anytime within the task period)
    for ipatt_index = 1:size(rf_block2_pos_generated_vortex,1)
        if ~isempty(rf_block2_pos_generated_vortex{1})
            ipatt = rf_block2_pos_generated_vortex{ipatt_index,1};
            for time = rf_block2_pos_generated_vortex{ipatt_index,3}
                temp1 = full(spiral_filt_pos{ipatt,time});
                RF_block2_sub_dist.generated.pos(:,:,time) =  RF_block2_sub_dist.generated.pos(:,:,time) + temp1;
            end
        end
    end


    for ipatt_index = 1:size(rf_block2_nega_generated_vortex,1)
        if ~isempty(rf_block2_nega_generated_vortex{1})
            ipatt = rf_block2_nega_generated_vortex{ipatt_index,1};
            for time = rf_block2_nega_generated_vortex{ipatt_index,3}
                temp1 = full(spiral_filt_nega{ipatt,time});
                RF_block2_sub_dist.generated.neg(:,:,time) = RF_block2_sub_dist.generated.neg(:,:,time) + temp1;
            end
        end
    end

%%

% LEFT toe Movement: Subject Averaged spiral distribution
    LF_block1_sub_dist = struct();
    LF_block2_sub_dist = struct();

%--> Modulated Spirals (exist prior to task onset and persist for more than 3 frames)
    LF_block1_sub_dist.modulated.pos = zeros(175,251,283);
    LF_block1_sub_dist.modulated.neg = zeros(175,251,283);

    LF_block1_sub_dist.generated.pos = zeros(175,251,283);
    LF_block1_sub_dist.generated.neg = zeros(175,251,283);


    for ipatt_index = 1:size(lf_block1_pos_modulated_vortex,1)
        if ~isempty(lf_block1_pos_modulated_vortex{1})
            ipatt = lf_block1_pos_modulated_vortex{ipatt_index,1};
            for time = lf_block1_pos_modulated_vortex{ipatt_index,3}
                temp1 = full(spiral_filt_pos{ipatt,time});
                LF_block1_sub_dist.modulated.pos(:,:,time) =  LF_block1_sub_dist.modulated.pos(:,:,time) + temp1;
            end
        end
    end


    for ipatt_index = 1:size(lf_block1_nega_modulated_vortex,1)
        if ~isempty(lf_block1_nega_modulated_vortex{1})
            ipatt = lf_block1_nega_modulated_vortex{ipatt_index,1};
            for time = lf_block1_nega_modulated_vortex{ipatt_index,3}
                temp1 = full(spiral_filt_nega{ipatt,time});
                LF_block1_sub_dist.modulated.neg(:,:,time) = LF_block1_sub_dist.modulated.neg(:,:,time) + temp1;
            end
        end
    end

    %--> Generated Spirals (Appear anytime within the task period)
    for ipatt_index = 1:size(lf_block1_pos_generated_vortex,1)
        if ~isempty(lf_block1_pos_generated_vortex{1})
            ipatt = lf_block1_pos_generated_vortex{ipatt_index,1};
            for time = lf_block1_pos_generated_vortex{ipatt_index,3}
                temp1 = full(spiral_filt_pos{ipatt,time});
                LF_block1_sub_dist.generated.pos(:,:,time) =  LF_block1_sub_dist.generated.pos(:,:,time) + temp1;
            end
        end
    end


    for ipatt_index = 1:size(lf_block1_nega_generated_vortex,1)
        if ~isempty(lf_block1_nega_generated_vortex{1})
            ipatt = lf_block1_nega_generated_vortex{ipatt_index,1};
            for time = lf_block1_nega_generated_vortex{ipatt_index,3}
                temp1 = full(spiral_filt_nega{ipatt,time});
                LF_block1_sub_dist.generated.neg(:,:,time) = LF_block1_sub_dist.generated.neg(:,:,time) + temp1;
            end
        end
    end


      % BLOCK 2

 %--> Modulated Spirals (exist prior to task onset and persist for more than 3 frames)
    LF_block2_sub_dist.modulated.pos = zeros(175,251,283);
    LF_block2_sub_dist.modulated.neg = zeros(175,251,283);

    LF_block2_sub_dist.generated.pos = zeros(175,251,283);
    LF_block2_sub_dist.generated.neg = zeros(175,251,283);

    
    for ipatt_index = 1:size(lf_block2_pos_modulated_vortex,1)
        if ~isempty(lf_block2_pos_modulated_vortex{1})
            ipatt = lf_block2_pos_modulated_vortex{ipatt_index,1};
            for time = lf_block2_pos_modulated_vortex{ipatt_index,3}
                temp1 = full(spiral_filt_pos{ipatt,time});
                LF_block2_sub_dist.modulated.pos(:,:,time) =  LF_block2_sub_dist.modulated.pos(:,:,time) + temp1;
            end
        end
    end

    
    for ipatt_index = 1:size(lf_block2_nega_modulated_vortex,1)
        if ~isempty(lf_block2_nega_modulated_vortex{1})
            ipatt = lf_block2_nega_modulated_vortex{ipatt_index,1};
            for time = lf_block2_nega_modulated_vortex{ipatt_index,3}
                temp1 = full(spiral_filt_nega{ipatt,time});
                LF_block2_sub_dist.modulated.neg(:,:,time) = LF_block2_sub_dist.modulated.neg(:,:,time) + temp1;
            end
        end
    end


    %--> Generated Spirals (Appear anytime within the task period)
    for ipatt_index = 1:size(lf_block2_pos_generated_vortex,1)
        if ~isempty(lf_block2_pos_generated_vortex{1})
            ipatt = lf_block2_pos_generated_vortex{ipatt_index,1};
            for time = lf_block2_pos_generated_vortex{ipatt_index,3}
                temp1 = full(spiral_filt_pos{ipatt,time});
                LF_block2_sub_dist.generated.pos(:,:,time) =  LF_block2_sub_dist.generated.pos(:,:,time) + temp1;
            end
        end
    end


    for ipatt_index = 1:size(lf_block2_nega_generated_vortex,1)
        if ~isempty(lf_block2_nega_generated_vortex{1})
            ipatt = lf_block2_nega_generated_vortex{ipatt_index,1};
            for time = lf_block2_nega_generated_vortex{ipatt_index,3}
                temp1 = full(spiral_filt_nega{ipatt,time});
                LF_block2_sub_dist.generated.neg(:,:,time) = LF_block2_sub_dist.generated.neg(:,:,time) + temp1;
            end
        end
    end

%% 
% Tongue Movement: Subject Averaged spiral distribution
    T_block1_sub_dist = struct();
    T_block2_sub_dist = struct();

%--> Modulated Spirals (exist prior to task onset and persist for more than 3 frames)
    T_block1_sub_dist.modulated.pos = zeros(175,251,283);
    T_block1_sub_dist.modulated.neg = zeros(175,251,283);

    T_block1_sub_dist.generated.pos = zeros(175,251,283);
    T_block1_sub_dist.generated.neg = zeros(175,251,283);
    
    
    for ipatt_index = 1:size(t_block1_pos_modulated_vortex,1)
        if ~isempty(t_block1_pos_modulated_vortex{1})
            ipatt = t_block1_pos_modulated_vortex{ipatt_index,1};
            for time = t_block1_pos_modulated_vortex{ipatt_index,3}
                temp1 = full(spiral_filt_pos{ipatt,time});
                T_block1_sub_dist.modulated.pos(:,:,time) =  T_block1_sub_dist.modulated.pos(:,:,time) + temp1;
            end
        end
    end
    
    
    for ipatt_index = 1:size(t_block1_nega_modulated_vortex,1)
        if ~isempty(t_block1_nega_modulated_vortex{1})
            ipatt = t_block1_nega_modulated_vortex{ipatt_index,1};
            for time = t_block1_nega_modulated_vortex{ipatt_index,3}
                temp1 = full(spiral_filt_nega{ipatt,time});
                T_block1_sub_dist.modulated.neg(:,:,time) = T_block1_sub_dist.modulated.neg(:,:,time) + temp1;
            end
        end
    end

    %--> Generated Spirals (Appear anytime within the task period)
    for ipatt_index = 1:size(t_block1_pos_generated_vortex,1)
        if ~isempty(t_block1_pos_generated_vortex{1})
            ipatt = t_block1_pos_generated_vortex{ipatt_index,1};
            for time = t_block1_pos_generated_vortex{ipatt_index,3}
                temp1 = full(spiral_filt_pos{ipatt,time});
                T_block1_sub_dist.generated.pos(:,:,time) =  T_block1_sub_dist.generated.pos(:,:,time) + temp1;
            end
        end
    end


    for ipatt_index = 1:size(t_block1_nega_generated_vortex,1)
        if ~isempty(t_block1_nega_generated_vortex{1})
            ipatt = t_block1_nega_generated_vortex{ipatt_index,1};
            for time = t_block1_nega_generated_vortex{ipatt_index,3}
                temp1 = full(spiral_filt_nega{ipatt,time});
                T_block1_sub_dist.generated.neg(:,:,time) = T_block1_sub_dist.generated.neg(:,:,time) + temp1;
            end
        end
    end


      % BLOCK 2

 %--> Modulated Spirals (exist prior to task onset and persist for more than 3 frames)
    T_block2_sub_dist.modulated.pos = zeros(175,251,283);
    T_block2_sub_dist.modulated.neg = zeros(175,251,283);

    T_block2_sub_dist.generated.pos = zeros(175,251,283);
    T_block2_sub_dist.generated.neg = zeros(175,251,283);
    
    
    for ipatt_index = 1:size(t_block2_pos_modulated_vortex,1)
        if ~isempty(t_block2_pos_modulated_vortex{1})
            ipatt = t_block2_pos_modulated_vortex{ipatt_index,1};
            for time = t_block2_pos_modulated_vortex{ipatt_index,3}
                temp1 = full(spiral_filt_pos{ipatt,time});
                T_block2_sub_dist.modulated.pos(:,:,time) =  T_block2_sub_dist.modulated.pos(:,:,time) + temp1;
            end
        end
    end

    
    for ipatt_index = 1:size(t_block2_nega_modulated_vortex,1)
        if ~isempty(t_block2_nega_modulated_vortex{1})
            ipatt = t_block2_nega_modulated_vortex{ipatt_index,1};
            for time = t_block2_nega_modulated_vortex{ipatt_index,3}
                temp1 = full(spiral_filt_nega{ipatt,time});
                T_block2_sub_dist.modulated.neg(:,:,time) = T_block2_sub_dist.modulated.neg(:,:,time) + temp1;
            end
        end
    end

  
    %--> Generated Spirals (Appear anytime within the task period)
    for ipatt_index = 1:size(t_block2_pos_generated_vortex,1)
        if ~isempty(t_block2_pos_generated_vortex{1})
            ipatt = t_block2_pos_generated_vortex{ipatt_index,1};
            for time = t_block2_pos_generated_vortex{ipatt_index,3}
                temp1 = full(spiral_filt_pos{ipatt,time});
                T_block2_sub_dist.generated.pos(:,:,time) =  T_block2_sub_dist.generated.pos(:,:,time) + temp1;
            end
        end
    end


    for ipatt_index = 1:size(t_block2_nega_generated_vortex,1)
        if ~isempty(t_block2_nega_generated_vortex{1})
            ipatt = t_block2_nega_generated_vortex{ipatt_index,1};
            for time = t_block2_nega_generated_vortex{ipatt_index,3}
                temp1 = full(spiral_filt_nega{ipatt,time});
                T_block2_sub_dist.generated.neg(:,:,time) = T_block2_sub_dist.generated.neg(:,:,time) + temp1;
            end
        end
    end


        foldername = [main_folder,'/Sample Data/Motor Task/Analysis/'];
        if hemisphere == 1
        filename = ['task_specific_spiral_dist_leftHem_sub',num2str(subject),'.mat'];
        elseif hemisphere == 2
            filename = ['task_specific_spiral_dist_rightHem_sub',num2str(subject),'.mat'];
        end
        save([foldername,filename],'RH_block1_sub_dist','RH_block2_sub_dist','LH_block1_sub_dist','LH_block2_sub_dist' ,...
            'RF_block1_sub_dist','RF_block2_sub_dist','LF_block1_sub_dist','LF_block2_sub_dist','T_block2_sub_dist','T_block1_sub_dist');
    
    cd(main_folder)

end

end

