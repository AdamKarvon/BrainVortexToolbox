function [classification_accuracy_avg, classification_accuracy_stderr, classificaton_accuracy] = spiral_classifer_motor_task()
restoredefaultpath
cd '\\suphys.physics.usyd.edu.au\akar5239\BrainVortexToolbox-main'
addpath(genpath([pwd]))
main_folder = pwd;

flagSur = 0;%  0 for real data, 1 to generate surrogate data

% a range of parameters avaiable for different dataset, but for demonstration
% purpose, only use the parameters provided
No_of_Subject = 1; % number of subjects used for analysis, randomly selected from HCP database (S1200)
flagRest = 0; %  resting data => 1 , task data => 0

hemisphere = 2; % 1 for left hemisphere, 2 for right hemisphere, keep this as 1 for demonstration purpose

% flagTask:
% 1 = language task, original 100 subjects;
% 2 = language task, additional 100 subjects;
% 3 = working memory task;
% 4 = Motor Task
flagTask = 4;

% 0 = unsmoothed, raw data; 1 = temporally smoothed (bandpass filtered)
% data; 2 = spatiotemporally smoothed (bandpass filtered) data
flagSmooth = 1;

load('parcellation_template.mat')

%% load task label of each subject

foldername = [main_folder,'/Sample Data/Motor Task/TaskLabel'];
cd(foldername)
% name = dir(pwd);
load (['MotorTaskLabelAllSubject.mat']);
for isubject = 1:150
    fullTime_allsubject{isubject} = TaskLabel_AllSubject_motor{isubject};
end

foldername = [main_folder,'/Sample Data/Motor Task/Analysis'];
cd(foldername)
name = dir(pwd) ;
file_name2 = ['MotorTaskLabelAllSubject.mat'];
load (file_name2);
for isubject = 1:150
    fullTime_allsubject{isubject} = TaskLabel_AllSubject_motor{isubject};
end
%% load preprocessed fMRI data files
spiral_template_posi_nega_sparse = [];

for subject = 1:150 
    
    if hemisphere == 1
        foldername = [main_folder,'/Sample Data/Motor Task/Spiral Detected'];
        cd(foldername)
        filename = ['Spiral_detected_surfilt_motor_task_LEFT_sub',num2str(subject),'.mat'];
        load(filename,'spiral_filt_pos_real_centreONLY_95perc_extend','spiral_filt_nega_real_centreONLY_95perc_extend') 
    elseif hemisphere == 2
        foldername = [main_folder,'/Sample Data/Motor Task/Spiral Detected'];
        cd(foldername)
        filename = ['Spiral_detected_surfilt_motor_task_RIGHT_sub',num2str(subject),'.mat'];
        load(filename,'spiral_filt_pos_real_centreONLY_95perc_extend','spiral_filt_nega_real_centreONLY_95perc_extend')         
    end

% combine both clockwise and anticlockwise spirals in the same map across trials and subjects
    
    x = 1:251;
    y = 1:175;
    [x_grid, y_grid]= meshgrid(x,y);
    radius_threshold = 5;    % mark the 5 radius area surrounding spiral centre with 1 or -1
    for time = 1:size(spiral_filt_pos_real_centreONLY_95perc_extend,2)
        spiral_template_pos = zeros(175,251);
        spiral_template_nega = zeros(175,251);
        if time <= size(spiral_filt_pos_real_centreONLY_95perc_extend,2)

            % anticlockwise spirals
            for ipatt = 1:size(spiral_filt_pos_real_centreONLY_95perc_extend,1);
                temp1 = spiral_filt_pos_real_centreONLY_95perc_extend{ipatt,time};
                if nansum(temp1(:)) == 0; % skip empty cells
                    continue
                else
                    row = temp1(2); % find the row coordinate of spiral centre
                    col = temp1(1); % find the row coordinate of spiral centre
                    distance_matrix = sqrt((row - y_grid).^2 + (col - x_grid).^2);
                    distance_matrix(distance_matrix>radius_threshold) = nan;
                    idx_extend = find(~isnan(distance_matrix));
                    spiral_template_pos(idx_extend) = 1; % extend the region to 5x5 surrounding the random point and asign the same value as the centre point
                end
            end
        end

        % clockwise spirals
        if time <= size(spiral_filt_nega_real_centreONLY_95perc_extend,2)
            for ipatt = 1:size(spiral_filt_nega_real_centreONLY_95perc_extend,1);
                temp1 = spiral_filt_nega_real_centreONLY_95perc_extend{ipatt,time};


                if nansum(temp1(:)) == 0; % skip empty cells
                    continue
                else
                    row = temp1(2); % find the row coordinate of spiral centre
                    col = temp1(1); % find the row coordinate of spiral centre
                    distance_matrix = sqrt((row - y_grid).^2 + (col - x_grid).^2);
                    distance_matrix(distance_matrix>radius_threshold) = nan;
                    idx_extend = find(~isnan(distance_matrix));
                    spiral_template_nega(idx_extend) = -1; % extend the region to 5x5 surrounding the random point and asign the same value as the centre point
                end
            end
        end
        % calculate the combined distribution of both clockwise and
    % anticlockwise spirals for each subject at each time point
    spiral_template_posi_nega = spiral_template_nega + spiral_template_pos;
    spiral_template_posi_nega(isnan(spiral_template_posi_nega)) = 0; 
    spiral_template_posi_nega_sparse{time,subject} = spiral_template_posi_nega;
    
 end
end





%%

for subject = 1:150

session_duration = 17;   % (currently this is the number of frames 12 seconds/0.72 )

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
    %-> will also need to extract relevant spiral sizes DONE

    % Movement periods for the right hand
    if size(start_end_time_right_hand,1) == 2
        block_1 = start_end_time_right_hand(1,1):start_end_time_right_hand(1,1)+session_duration;
        block_2 = start_end_time_right_hand(2,1):start_end_time_right_hand(2,1)+session_duration;

        spiral_distribution_right_hand_block1{subject} = spiral_template_posi_nega_sparse(block_1,subject);
        spiral_distribution_right_hand_block2{subject} = spiral_template_posi_nega_sparse(block_2,subject);
    end



    % Movement periods for the right foot
    if size(start_end_time_right_foot,1) == 2
        block_1 = start_end_time_right_foot(1,1):start_end_time_right_foot(1,1)+session_duration;
        block_2 = start_end_time_right_foot(2,1):start_end_time_right_foot(2,1)+session_duration;

        spiral_distribution_right_foot_block1{subject} = spiral_template_posi_nega_sparse(block_1,subject);
        spiral_distribution_right_foot_block2{subject} = spiral_template_posi_nega_sparse(block_2,subject);
    end

    % Movement periods for the tongue
    if size(start_end_time_tongue,1) == 2
        block_1 = start_end_time_tongue(1,1):start_end_time_tongue(1,1)+session_duration;
        block_2 = start_end_time_tongue(2,1):start_end_time_tongue(2,1)+session_duration;

        spiral_distribution_tongue_block1{subject} = spiral_template_posi_nega_sparse(block_1,subject);

        spiral_distribution_tongue_block2{subject} = spiral_template_posi_nega_sparse(block_2,subject);
    end

    % extract task specifc spiral distribution data: "Left" tasks


    % Movement Periods for the Left Foot
    if size(start_end_time_left_foot,1) == 2
        block_1 = start_end_time_left_foot(1,1):start_end_time_left_foot(1,1)+session_duration;
        block_2 = start_end_time_left_foot(2,1):start_end_time_left_foot(2,1)+session_duration;

        spiral_distribution_left_foot_block1{subject} = spiral_template_posi_nega_sparse(block_1,subject);
        spiral_distribution_left_foot_block2{subject} = spiral_template_posi_nega_sparse(block_2,subject);
    end

    % Movement Periods for the Left Hand
    if size(start_end_time_left_hand,1) == 2
        block_1 = start_end_time_left_hand(1,1):start_end_time_left_hand(1,1)+session_duration;
        block_2 = start_end_time_left_hand(2,1):start_end_time_left_hand(2,1)+session_duration;

        spiral_distribution_left_hand_block1{subject} = spiral_template_posi_nega_sparse(block_1,subject);
        spiral_distribution_left_hand_block2{subject} = spiral_template_posi_nega_sparse(block_2,subject);
    end


end
%%
spiral_distribution_right_hand_matrix = zeros(175,251,18,300);
spiral_distribution_left_hand_matrix = zeros(175,251,18,300);
spiral_distribution_right_foot_matrix = zeros(175,251,18,300);
spiral_distribution_left_foot_matrix = zeros(175,251,18,300);
spiral_distribution_tongue_matrix = zeros(175,251,18,300);
for time = 1:18
for trial = 1:300
    temp1 = [cellfun(@(c) c{time}(1:175,1:251), spiral_distribution_right_hand_block1, 'uniformoutput',false), cellfun(@(c) c{time}(1:175,1:251), spiral_distribution_right_hand_block2, 'uniformoutput',false)];
    spiral_distribution_right_hand_matrix(:,:,time,trial) = full(temp1{trial});

    temp1 = [cellfun(@(c) c{time}(1:175,1:251), spiral_distribution_right_foot_block1, 'uniformoutput',false), cellfun(@(c) c{time}(1:175,1:251), spiral_distribution_right_foot_block2, 'uniformoutput',false)];
    spiral_distribution_right_foot_matrix(:,:,time,trial) = full(temp1{trial});

    temp1 = [cellfun(@(c) c{time}(1:175,1:251), spiral_distribution_tongue_block1, 'uniformoutput',false), cellfun(@(c) c{time}(1:175,1:251), spiral_distribution_tongue_block2, 'uniformoutput',false)];
    spiral_distribution_tongue_matrix(:,:,time,trial) = full(temp1{trial});

    temp1 = [cellfun(@(c) c{time}(1:175,1:251), spiral_distribution_left_foot_block1, 'uniformoutput',false), cellfun(@(c) c{time}(1:175,1:251), spiral_distribution_left_foot_block2, 'uniformoutput',false)];
    spiral_distribution_left_foot_matrix(:,:,time,trial) = full(temp1{trial});

    temp1 = [cellfun(@(c) c{time}(1:175,1:251), spiral_distribution_left_hand_block1, 'uniformoutput',false), cellfun(@(c) c{time}(1:175,1:251), spiral_distribution_left_hand_block2, 'uniformoutput',false)];
    spiral_distribution_left_hand_matrix(:,:,time,trial) = full(temp1{trial});

end
end

%%

    % find no. of trials for ALL tasks
    trials_right_hand = size(spiral_distribution_right_hand_matrix ,4);
    trials_right_foot = size(spiral_distribution_right_foot_matrix,4);

    trials_tongue = size(spiral_distribution_tongue_matrix,4);

    trials_left_hand = size(spiral_distribution_left_hand_matrix,4);
    trials_left_foot = size(spiral_distribution_left_foot_matrix,4);
    

    time_frame =  6:14;
    right_hand_classification_accuracy =[];
    right_foot_classification_accuracy = [];
    tongue_classification_accuracy = [];
    left_hand_classification_accuracy = [];
    left_foot_classification_accuracy = [];

        for repetition = 1:100

            % 80% trials for training data
            train_trial_right_hand = randperm(trials_right_hand,round(0.8*trials_right_hand));
            train_trial_right_foot = randperm(trials_right_foot,round(0.8*trials_right_foot));
            train_trial_tongue = randperm(trials_tongue,round(0.8*trials_tongue));
            train_trial_left_hand = randperm(trials_left_hand,round(0.8*trials_left_hand));
            train_trial_left_foot = randperm(trials_left_foot,round(0.8*trials_left_foot));

            % 20% remaining trials for testing
            test_trial_right_hand = setdiff([1:trials_right_hand],train_trial_right_hand);
            test_trial_right_foot = setdiff([1:trials_right_foot],train_trial_right_foot);
            test_trial_tongue = setdiff([1:trials_tongue],train_trial_tongue);
            test_trial_left_hand = setdiff([1:trials_left_hand],train_trial_left_hand);
            test_trial_left_foot = setdiff([1:trials_left_foot],train_trial_left_foot);


            % calculate the mean spiral distribution for training trials
            % for each task condition, as training template

            temp1_right_hand_train_template = nanmean(nanmean(spiral_distribution_right_hand_matrix(:,:,time_frame,train_trial_right_hand),4),3);
            temp1_right_foot_train_template = nanmean(nanmean(spiral_distribution_right_foot_matrix(:,:,time_frame,train_trial_right_foot),4),3);
            temp1_tongue_train_template = nanmean(nanmean(spiral_distribution_tongue_matrix(:,:,time_frame,train_trial_tongue),4),3);
            temp1_left_hand_train_template = nanmean(nanmean(spiral_distribution_left_hand_matrix(:,:,time_frame,train_trial_left_hand),4),3);
            temp1_left_foot_train_template = nanmean(nanmean(spiral_distribution_left_foot_matrix(:,:,time_frame,train_trial_left_foot),4),3);


            %select random test trial

            temp1_right_hand_test(:,:,:) = nanmean(spiral_distribution_right_hand_matrix(:,:,time_frame,test_trial_right_hand),3);
            temp1_right_foot_test(:,:,:) = nanmean(spiral_distribution_right_foot_matrix(:,:,time_frame,test_trial_right_foot),3);
            temp1_tongue_test(:,:,:) = nanmean(spiral_distribution_tongue_matrix(:,:,time_frame,test_trial_tongue),3);
            temp1_left_hand_test(:,:,:) = nanmean(spiral_distribution_left_hand_matrix(:,:,time_frame,test_trial_left_hand),3);
            temp1_left_foot_test(:,:,:) = nanmean(spiral_distribution_left_foot_matrix(:,:,time_frame,test_trial_left_foot),3);

            % Initialize arrays
            right_hand_classification_counts = zeros(1, 3); % [righthand, rightfoot, tongue, lefthand, leftfoot]
            right_foot_classification_counts = zeros(1, 3);
            tongue_classification_counts = zeros(1, 3);
            left_hand_classification_counts = zeros(1, 3);
            left_foot_classification_counts = zeros(1, 3);

            for trial = 1:60
                % Reset for each trial
                % % maxTimeFramesEachTrial contains the time frame of the maximum correlation for each trial
                % % maxCategoryEachTrial contains the category of the maximum correlation for each trial
                % % countMax contains the count of times each category had the overall max correlation across all trials
                % overallMaxCorrelation = 0;
                % overallMaxCategory = 0;
                % overallMaxTimeFrame = 0;


                correlations = zeros(1,3);


                if hemisphere == 1
                    righthand = corr2(temp1_right_hand_train_template, temp1_right_hand_test(:,:,trial));
                    rightfoot = corr2(temp1_right_foot_train_template, temp1_right_hand_test(:,:,trial));
                    tongue = corr2(temp1_tongue_train_template, temp1_right_hand_test(:,:,trial));

                    correlations = [righthand, rightfoot, tongue];

                    % Find overall max correlation and corresponding category and time frame
                    [~, maxIdx] = max(correlations);

                    % Increment count for the category with overall max correlation
                    right_hand_classification_counts(maxIdx) = right_hand_classification_counts(maxIdx) + 1;


                    righthand = corr2(temp1_right_hand_train_template, temp1_right_foot_test(:,:,trial));
                    rightfoot = corr2(temp1_right_foot_train_template, temp1_right_foot_test(:,:,trial));
                    tongue = corr2(temp1_tongue_train_template,temp1_right_foot_test(:,:,trial));


                    correlations = [righthand, rightfoot, tongue];
                    [~, maxIdx] = max(correlations);

                    right_foot_classification_counts(maxIdx) = right_foot_classification_counts(maxIdx) + 1;



                    righthand = corr2(temp1_right_hand_train_template, temp1_tongue_test(:,:,trial));
                    rightfoot = corr2(temp1_right_foot_train_template, temp1_tongue_test(:,:,trial));
                    tongue = corr2(temp1_tongue_train_template, temp1_tongue_test(:,:,trial));


                    correlations = [righthand, rightfoot, tongue];
                    [~, maxIdx] = max(correlations);

                    tongue_classification_counts(maxIdx) = tongue_classification_counts(maxIdx) + 1;


                elseif hemisphere == 2
                    tongue = corr2(temp1_tongue_train_template, temp1_tongue_test(:,:,trial));
                    lefthand = corr2(temp1_left_hand_train_template, temp1_tongue_test(:,:,trial));
                    leftfoot = corr2(temp1_left_foot_train_template, temp1_tongue_test(:,:,trial));

                    correlations = [lefthand, leftfoot, tongue];
                    [~, maxIdx] = max(correlations);

                    tongue_classification_counts(maxIdx) = tongue_classification_counts(maxIdx) + 1;

                    tongue = corr2(temp1_tongue_train_template, temp1_left_hand_test(:,:,trial));
                    lefthand = corr2(temp1_left_hand_train_template, temp1_left_hand_test(:,:,trial));
                    leftfoot = corr2(temp1_left_foot_train_template, temp1_left_hand_test(:,:,trial));

                    correlations = [lefthand, leftfoot, tongue];
                    [~, maxIdx] = max(correlations);

                    left_hand_classification_counts(maxIdx) = left_hand_classification_counts(maxIdx) + 1;


                    tongue = corr2(temp1_tongue_train_template, temp1_left_foot_test(:,:,trial));
                    lefthand = corr2(temp1_left_hand_train_template,temp1_left_foot_test(:,:,trial));
                    leftfoot = corr2(temp1_left_foot_train_template, temp1_left_foot_test(:,:,trial));

                    correlations = [lefthand, leftfoot, tongue];
                    [~, maxIdx] = max(correlations);

                    left_foot_classification_counts(maxIdx) = left_foot_classification_counts(maxIdx) + 1;

                end

            end
            
            if hemisphere == 1
                right_hand_classification_accuracy = [right_hand_classification_accuracy right_hand_classification_counts(1)/60];

                right_foot_classification_accuracy = [right_foot_classification_accuracy right_foot_classification_counts(2)/60];
            elseif hemisphere == 2
                left_hand_classification_accuracy = [left_hand_classification_accuracy left_hand_classification_counts(1)/60];
                %
                left_foot_classification_accuracy = [left_foot_classification_accuracy left_foot_classification_counts(2)/60];
            end

            tongue_classification_accuracy = [tongue_classification_accuracy tongue_classification_counts(3)/60];
        end

        %%
        classification_accuracy_avg = [mean(right_hand_classification_accuracy), mean(right_foot_classification_accuracy), mean(tongue_classification_accuracy)].*100;
        classification_accuracy_std = [std(right_hand_classification_accuracy.*100), std(right_foot_classification_accuracy.*100), std(tongue_classification_accuracy.*100)];
        classification_accuracy_stderr = classification_accuracy_std./sqrt(size(tongue_classification_accuracy,2));


% Data
means = classification_accuracy_avg;
std_errors = classification_accuracy_stderr;
tasks = {'Right Hand', 'Right Foot', 'Tongue'};

num_tasks = length(tasks);

% Create bar graph
figure;
b = bar(means, 'FaceColor', 'white');

% Add error bars
hold on;
% Assuming means is a 1D array, we'll directly use its length
numbars = length(means);
for i = 1:numbars
    x = b.XData(i); % Get the x-coordinate of the center of each bar
    errorbar(x, means(i), std_errors(i), 'k', 'linestyle', 'none');
end

% Overlay scatter plot
scatterData = {right_hand_classification_accuracy.*100, right_foot_classification_accuracy.*100, tongue_classification_accuracy.*100};
for i = 1:num_tasks
    scatterX = repmat(i, 1, length(scatterData{i})) + (rand(size(scatterData{i})) - 0.5) * 0.1; % Adjust 0.1 for more/less spread
    scatter(scatterX, scatterData{i}, 'filled', 'k', 'MarkerFaceAlpha',0.3);
end
% Add a red dotted line at 33% for chance level classification
yline(33, '--r', 'Chance Level');
hold off;

% Customize the graph
set(gca, 'xticklabel', tasks);
ylabel('Mean Classification Accuracy (%)');
xlabel('Task');
title('Left Hemisphere Classification Accuracy for Each Task');



        % right_hand_classification_accuracy_avg(time_frame) = mean(right_hand_classification_accuracy);
        % 
        % right_foot_classification_accuracy_avg(time_frame) = mean(right_foot_classification_accuracy);
        % 
        % tongue_classification_accuracy_avg(time_frame) = mean(tongue_classification_accuracy);
        % % 
        % % left_hand_classification_accuracy_avg(time_frame) = mean(left_hand_classification_accuracy);
        % % 
        % % left_foot_classification_accuracy_avg(time_frame) = mean(left_foot_classification_accuracy);
 end















