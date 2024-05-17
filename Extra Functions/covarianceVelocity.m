cd '\\suphys.physics.usyd.edu.au\akar5239\BrainVortexToolbox-main'
addpath(genpath([pwd]))
main_folder = pwd;

movement_block = 1;
task = 'RH';

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
task_times = fullTime_allsubject{isubject};
disp(['loading task data...'])
if movement_block == 1
    switch task
        case 'RH'
            task_onset = task_times(1,1);
            task_end = task_times(1,2);
            load ('Task_specific_spiralData_leftHEM.mat', 'RH_block1');
            task_data = RH_block1;
        case 'LF'
            task_onset = task_times(2,1);
            task_end = task_times(2,2);
            load ('Task_specific_spiralData_leftHEM.mat', 'LF_block1');
            task_data = LF_block1;
        case 'T'
            task_onset = task_times(3,1);
            task_end = task_times(3,2);
            load ('Task_specific_spiralData_leftHEM.mat', 'T_block1');
            task_data = T_block1;
        case 'RF'
            task_onset = task_times(4,1);
            task_end = task_times(4,2);
            load ('Task_specific_spiralData_leftHEM.mat', 'RF_block1');
            task_data = RF_block1;
        case 'LH'
            task_onset = task_times(5,1);
            task_end = task_times(5,2);
            load ('Task_specific_spiralData_leftHEM.mat', 'LH_block1');
            task_data = LH_block1;
    end

elseif movement_block == 2
    switch task

        case 'T'
            task_onset = task_times(6,1);
            task_end = task_times(6,2);
            load ('Task_specific_spiralData_leftHEM.mat', 'T_block2');
            task_data = T_block2;
        case 'LF'
            task_onset = task_times(7,1);
            task_end = task_times(7,2);
            load ('Task_specific_spiralData_leftHEM.mat', 'LF_block2');
            task_data = LF_block2;
        case 'RH'
            task_onset = task_times(8,1);
            task_end = task_times(8,2);
            load ('Task_specific_spiralData_leftHEM.mat', 'RH_block2');
            task_data = RH_block2;
        case 'LH'
            task_onset = task_times(9,1);
            task_end = task_times(9,2);
            load ('Task_specific_spiralData_leftHEM.mat', 'LH_block2');
            task_data = LH_block2;
        case 'RF'
            task_onset = task_times(10,1);
            task_end = task_times(10,2);
            load ('Task_specific_spiralData_leftHEM.mat', 'RF_block2');
            task_data = RF_block2;
    end

end

 %% Velocity Covariance matrix
% 
% Compare spirals at each timeslot for individual's spiral distribution (283 frames)
% Compute the Covariance Velocity_Cov{1:283} = Cov(A,B) 
% where: 
% A = subject i =[Velocity at each timestep x all spiral at frame [?]] 
% B = subject j= [Velocity at each timestep x all spiral at frame [?]] 
% 
% Velocity_Cov = spiral velocity Covariance matrix of every time frame for
% one subject

% Need to make 150 Velocity_Cov.

count_modulated_pos_modulated_pos = 0;
count_modulated_neg_modulated_neg = 0;
count_generated_neg_generated_neg = 0;
count_generated_pos_generated_pos = 0;
count_modulated_neg_generated_neg = 0;
count_modulated_pos_generated_pos = 0;
count_generated_pos_generated_neg = 0;
count_modulated_pos_modulated_neg = 0;
count_modulated_pos_generated_neg = 0;
count_generated_pos_modulated_neg = 0;





%%
% Assuming task_data.modulated.pos{subjectIndex} is the data structure
modulated_counter = 0 ;
modulated_pos_data = struct;


for subject = 1:150
% Note the order of the spiral types
subject_data = [task_data.modulated.pos{subject}; task_data.modulated.neg{subject};  task_data.generated.pos{subject}; task_data.generated.neg{subject}];  % Example for one subject

% Get sizes of each group
size_modulated_pos = size(task_data.modulated.pos{subject}, 1);
size_modulated_neg = size(task_data.modulated.neg{subject}, 1);
size_generated_pos = size(task_data.generated.pos{subject}, 1);
size_generated_neg = size(task_data.generated.neg{subject}, 1);

% Calculate cumulative sizes to determine row indices
cumulative_sizes = cumsum([0, size_modulated_pos, size_modulated_neg, size_generated_pos, size_generated_neg]);

% Store row indices for each group
indices_modulated_pos = (cumulative_sizes(1) + 1):cumulative_sizes(2);
indices_modulated_neg = (cumulative_sizes(2) + 1):cumulative_sizes(3);
indices_generated_pos = (cumulative_sizes(3) + 1):cumulative_sizes(4);
indices_generated_neg = (cumulative_sizes(4) + 1):cumulative_sizes(5);


%%

for modulated_spiral_id = indices_modulated_pos
    pos_modulated_velocity_data = subject_data{modulated_spiral_id,5};  
    x_vel = pos_modulated_velocity_data(1, :);
    y_vel = pos_modulated_velocity_data(2, :);
    pos_modulated_time_data = subject_data{modulated_spiral_id,3}; % Time frames for corresponding velocity measurements



    % Assuming you have a way to separate pre-task and post-task data, maybe based on time
    % For example, if task starts at time t = task_start_time
    
    % Adjust indices to exclude the last time frame for velocity calculation
    pre_task_indices = find(pos_modulated_time_data <= task_onset);
    post_task_indices = find(pos_modulated_time_data > task_onset);

    if isempty(post_task_indices) || size(post_task_indices,2) < 3

        continue

    else

    % Compute the lifetime of the object
    lifetime = length(pos_modulated_time_data);
    modulated_pos_data(subject).lifetime(modulated_spiral_id) = lifetime;

    
    % Ensure not to include the last time frame
    pre_task_indices = pre_task_indices(pre_task_indices < length(pos_modulated_time_data));
    post_task_indices = post_task_indices(post_task_indices < length(pos_modulated_time_data));
    
    
    % Pre-task data
    x_vel_pre = x_vel(pre_task_indices);
    y_vel_pre = y_vel(pre_task_indices);
    time_pre = pos_modulated_time_data(pre_task_indices);
    
    % Post-task data
    x_vel_post = x_vel(post_task_indices);  % Velocity data is not available for the last frame of movement
    y_vel_post = y_vel(post_task_indices);
    time_post = pos_modulated_time_data(post_task_indices);



    % Compute velocity magnitude and direction 
    vel_mag_pre = sqrt(x_vel_pre.^2 + y_vel_pre.^2);
    vel_dir_pre = atan2d(y_vel_pre, x_vel_pre);  % Direction in degrees
    
    vel_mag_post = sqrt(x_vel_post.^2 + y_vel_post.^2);
    vel_dir_post = atan2d(y_vel_post, x_vel_post);  % Direction in degrees

    % Compute mean velocity magnitude and direction before and after the task
    mean_vel_mag_pre = nanmean(vel_mag_pre);
    mean_vel_dir_pre = nanmean(vel_dir_pre);
    
    mean_vel_mag_post = nanmean(vel_mag_post);
    mean_vel_dir_post = nanmean(vel_dir_post);
    
    % Store the computed values for later analysis
    modulated_pos_data(subject).mean_vel_mag_pre(modulated_spiral_id) = mean_vel_mag_pre;
    modulated_pos_data(subject).mean_vel_dir_pre(modulated_spiral_id) = mean_vel_dir_pre;
    
    modulated_pos_data(subject).mean_vel_mag_post(modulated_spiral_id) = mean_vel_mag_post;
    modulated_pos_data(subject).mean_vel_dir_post(modulated_spiral_id) = mean_vel_dir_post;

    end
  

end




% Compute change in velocity magnitude
velocity_mag_change = modulated_pos_data(subject).mean_vel_mag_post - modulated_pos_data(subject).mean_vel_mag_pre;

modulated_pos_data(subject).velocity_mag_change = velocity_mag_change;


% Compute change in velocity direction
% Ensure the angle difference is within the range [-180, 180] degrees
direction_change = modulated_pos_data(subject).mean_vel_dir_post - modulated_pos_data(subject).mean_vel_dir_pre;
direction_change = mod(direction_change + 180, 360) - 180;  % Adjust angles to be within [-180, 180]

modulated_pos_data(subject).direction_change = direction_change;


%When interpreting the change in direction, it's important to note that a positive value indicates a counter-clockwise change, while a negative value indicates a clockwise change. The magnitude of the change tells you how much the direction has changed, not the actual final direction.



% Categorize and store the data based on object lifetime
for modulated_spiral_id = indices_modulated_pos
    lifetime = modulated_pos_data(subject).lifetime;
    
    short_lifetime{subject,:} = lifetime(:) <= 10;

    long_lifetime{subject,:} = lifetime(:) > 10;

end


modulated_pos_data(subject).longlived_velocitychange = modulated_pos_data(subject).velocity_mag_change(long_lifetime{subject});
modulated_pos_data(subject).shortlived_velocitychange = modulated_pos_data(subject).velocity_mag_change(short_lifetime{subject});

modulated_pos_data(subject).longlived_directionchange = modulated_pos_data(subject).direction_change(long_lifetime{subject});
modulated_pos_data(subject).shortlived_directionchange = modulated_pos_data(subject).direction_change(short_lifetime{subject});

% 
% %%
% total_spirals = size(subject_data, 1);
% 
% 
% 
% % Identify the common time frame range 
% % Step 1: Extract all time frames
% all_time_frames = [];
% for spiralIndex = 1:total_spirals
%     switch spiralIndex
%         case ismember(spiralIndex,indices_modulated_pos)
%             Spiral_type_ID = 1;
%         case ismember(spiralIndex,indices_modulated_neg)   
%             Spiral_type_ID = 2;
%         case ismember(spiralIndex,indices_generated_pos)
%             Spiral_type_ID = 3;
%         case ismember(spiralIndex,indices_generated_neg)
%             Spiral_type_ID = 4;
%     end
%     time_frames = [subject_data{spiralIndex, 3} Spiral_type_ID];  % Assuming time frames are in the 3rd column
%     all_time_frames = [all_time_frames; time_frames(:)];  % Flatten and append
% end
% 
% 
% % Step 2: Count occurrences of each time frame
% unique_time_frames = unique(all_time_frames(:,1));
% time_frame_counts = histcounts(all_time_frames,[unique_time_frames; max(unique_time_frames) + 1]);
% 
% % Step 3: Determine common start and end
% common_time_frames = unique_time_frames(time_frame_counts >= 2);
% common_start = min(common_time_frames);
% common_end = max(common_time_frames);
% 
% 
% % Initialize a matrix to store velocity covariance
% % The size of the matrix depends on the number of common time frames
% velocity_covariance = zeros(total_spirals, total_spirals);
% velocities_at_t = NaN(total_spirals,common_end - common_start + 1); 
% % Process each time frame
% for time = common_start:common_end
%      % To store velocities of all spirals at time t
% 
%     % Extract velocities for each spiral at time t
%     for spiralIndex = 1:total_spirals
%         spiral_velocity_data = subject_data{spiralIndex, 5};  % Assuming velocity data is in the 5th column
%         % Find the velocity at time t (you'll need a function to extract this)
%         t = find(subject_data{spiralIndex,3} == time);
%         if isempty(t) || t> size(subject_data{spiralIndex,5},2)
%             velocity_direction = NaN;
%         else 
%             velocity_x = subject_data{spiralIndex,5}(1,t);
%             velocity_y = subject_data{spiralIndex,5}(2,t);
% 
%             velocity_direction = atan2(velocity_y, velocity_x);
% 
% 
%         end
%         velocities_at_t(spiralIndex,time) = velocity_direction;
%     end
% 
% end
% 
% 
% 
% 
% for i = 1:size(velocities_at_t,1)
% for j = 1:size(velocities_at_t,1)
%     valid_indices_1 = ~isnan(velocities_at_t(i, :));
% 
%     valid_indices_2 = ~isnan(velocities_at_t(j, :));
% 
%     joint_valid_inices = valid_indices_1 + valid_indices_2 == 2;
% 
%     velocity1 = velocities_at_t(i,joint_valid_inices);
%     velocity2 = velocities_at_t(j,joint_valid_inices);
% 
%     if sum(joint_valid_inices) < 3
%         continue
%     else
%     Covariance = cov(velocity1,velocity2);
% 
%     velocity_covariance(i,j) = Covariance(1,2);
%     end
% 
% end
% end
% 
% % Calculate mean and standard deviation of off-diagonal elements
% off_diagonal_elements = velocity_covariance.*~eye(size(velocity_covariance));
% 
% mean_cov = mean(off_diagonal_elements,'all');
% std_dev_cov = std(off_diagonal_elements,0,'all');
% 
% % Define thresholds
% high_threshold = mean_cov + 4 * std_dev_cov; % Example: mean + 4*std deviation
% low_threshold = mean_cov - 4 * std_dev_cov; % Example: mean - 4*std deviation
% 
% 
% 
% plot_interaction = false;
% 
% [row,col] = find(tril(velocity_covariance,-1) > high_threshold);
% spiral_type1 = strings(length(row), 1);
% spiral_type2 = strings(length(row), 1);
% 
% for i = 1:size(row,1)
%     if ismember(row(i), indices_modulated_pos)
%         spiral_type1(i) = 'modulated_pos';
%     elseif ismember(row(i), indices_modulated_neg)
%         spiral_type1(i) = 'modulated_neg';
%     elseif ismember(row(i), indices_generated_pos)
%         spiral_type1(i) = 'generated_pos';
%     elseif ismember(row(i), indices_generated_neg)
%         spiral_type1(i) = 'generated_neg';
%     end
% 
%     if ismember(col(i), indices_modulated_pos)
%         spiral_type2(i) = 'modulated_pos';
%     elseif ismember(col(i), indices_modulated_neg)
%         spiral_type2(i) = 'modulated_neg';
%     elseif ismember(col(i), indices_generated_pos)
%         spiral_type2(i) = 'generated_pos';
%     elseif ismember(col(i), indices_generated_neg)
%         spiral_type2(i) = 'generated_neg';
%     end
% end
% 
% highCov_spiralPairs = [subject_data(row,:), subject_data(col,:), cellstr(spiral_type1), cellstr(spiral_type2)];
% 
% [row,col] = find(tril(velocity_covariance,-1) < low_threshold);
% spiral_type1 = strings(length(row), 1);
% spiral_type2 = strings(length(row), 1);
% for i = 1:size(row,1)
%     if ismember(row(i), indices_modulated_pos)
%         spiral_type1(i) = 'modulated_pos';
%     elseif ismember(row(i), indices_modulated_neg)
%         spiral_type1(i) = 'modulated_neg';
%     elseif ismember(row(i), indices_generated_pos)
%         spiral_type1(i) = 'generated_pos';
%     elseif ismember(row(i), indices_generated_neg)
%         spiral_type1(i) = 'generated_neg';
%     end
% 
%     if ismember(col(i), indices_modulated_pos)
%         spiral_type2(i) = 'modulated_pos';
%     elseif ismember(col(i), indices_modulated_neg)
%         spiral_type2(i) = 'modulated_neg';
%     elseif ismember(col(i), indices_generated_pos)
%         spiral_type2(i) = 'generated_pos';
%     elseif ismember(col(i), indices_generated_neg)
%         spiral_type2(i) = 'generated_neg';
%     end
% end
% 
% lowCov_spiralPairs = [subject_data(row,:) subject_data(col,:) cellstr(spiral_type1), cellstr(spiral_type2)];


% for spiralIndex = 1:size(highCov_spiralPairs,1)
%     % Use cellfun to extract x and y coordinates
%     spiral_path1 = highCov_spiralPairs{spiralIndex,2};
%     spiral_path2 = highCov_spiralPairs{spiralIndex,7};
%     minDistance = inf;
% 
%     % Loop through each point in the first cell array
%     for i = 1:length(spiral_path1)
%         point1 = spiral_path1{i};
%         % Loop through each point in the second cell array
%         for j = 1:length(spiral_path2)
%             point2 = spiral_path2{j};
% 
%             % Calculate the Euclidean distance between the points
%             distance = sqrt((point1(1) - point2(1))^2 + (point1(2) - point2(2))^2);
% 
%             % Update minimum distance if current distance is smaller
%             if distance < minDistance
%                 minDistance = distance;
%             end
%         end
%     end
% %%
%     if minDistance < 8
%         spiral_types = [highCov_spiralPairs{spiralIndex,11},'_',highCov_spiralPairs{spiralIndex,12}];
% 
%         switch spiral_types
%             case 'modulated_pos_modulated_pos'
%                 count_modulated_pos_modulated_pos = count_modulated_pos_modulated_pos + 1;
% 
%             case 'modulated_neg_modulated_neg'
%                 count_modulated_neg_modulated_neg = count_modulated_neg_modulated_neg + 1;
% 
%             case 'generated_neg_generated_neg'
%                 count_generated_neg_generated_neg = count_generated_neg_generated_neg + 1;
% 
%             case 'generated_pos_generated_pos'
%                 count_generated_pos_generated_pos = count_generated_pos_generated_pos + 1;
% 
%             case {'modulated_neg_generated_neg',  'generated_neg_modulated_neg'}
%                 count_modulated_neg_generated_neg = count_modulated_neg_generated_neg + 1;
% 
%             case {'modulated_pos_generated_pos',  'generated_pos_modulated_pos'}
%                 count_modulated_pos_generated_pos = count_modulated_pos_generated_pos + 1;
% 
%             case {'generated_pos_generated_neg',  'generated_neg_generated_pos'}
%                 count_generated_pos_generated_neg = count_generated_pos_generated_neg + 1;
% 
%             case {'modulated_pos_modulated_neg',  'modulated_neg_modulated_pos'}
%                 count_modulated_pos_modulated_neg = count_modulated_pos_modulated_neg + 1;
% 
%             case {'modulated_pos_generated_neg',  'generated_neg_modulated_pos'}
%                 count_modulated_pos_generated_neg = count_modulated_pos_generated_neg + 1;
% 
%             case {'generated_pos_modulated_neg',  'modulated_neg_generated_pos'}
%                 count_generated_pos_modulated_neg = count_generated_pos_modulated_neg + 1;
%         end
% 
% 
% 
% 
%         if plot_interaction == true
%             % Plot the line connecting all the points
%             figure(1); title('High Covariance')
%             hold on
% 
%             x_coords = cellfun(@(c) c(1), spiral_path1);
%             y_coords = cellfun(@(c) c(2), spiral_path1);
%             plot(x_coords, y_coords, 'k-');
%             plot(x_coords(1), y_coords(1),'o','MarkerFaceColor','b');
% 
%             x_coords = cellfun(@(c) c(1), spiral_path2);
%             y_coords = cellfun(@(c) c(2), spiral_path2);
%             plot(x_coords, y_coords, 'k-');
%             plot(x_coords(1), y_coords(1),'o','MarkerFaceColor','r');
% 
%             for parcellation_ID = 1:22
%                 parcellation_template_1par = parcellation_template;
%                 parcellation_template_1par(isnan(parcellation_template_1par)) = 0;
% 
%                 parcellation_template_1par(parcellation_template_1par~=parcellation_ID) = 0;
% 
%                 B = bwboundaries(parcellation_template_1par(:,:,1),'noholes');
%                 for k = 1:length(B)
%                     boundary = B{k};
%                     plot(boundary(:,2),boundary(:,1),'-','linewidth',2,'color',[0,0,0])
%                 end
%             end
% 
%             xlabel('X-coordinate');
%             ylabel('Y-coordinate');
%         end
%     end
% end

% for spiralIndex = 1:size(lowCov_spiralPairs,1)
%     % Use cellfun to extract x and y coordinates
%     spiral_path1 = lowCov_spiralPairs{spiralIndex,2};
%     spiral_path2 = lowCov_spiralPairs{spiralIndex,7};
%     minDistance = inf;
% 
% 
%     % Loop through each point in the first cell array
%     for i = 1:length(spiral_path1)
%         point1 = spiral_path1{i};
% 
%         % Loop through each point in the second cell array
%         for j = 1:length(spiral_path2)
%             point2 = spiral_path2{j};
% 
%             % Calculate the Euclidean distance between the points
%             distance = sqrt((point1(1) - point2(1))^2 + (point1(2) - point2(2))^2);
% 
%             % Update minimum distance if current distance is smaller
%             if distance < minDistance
%                 minDistance = distance;
%             end
%         end
%     end
% 
%     if minDistance < 8
%         spiral_types = [lowCov_spiralPairs{spiralIndex,11},'_',lowCov_spiralPairs{spiralIndex,12}];
% 
%         switch spiral_types
%             case 'modulated_pos_modulated_pos'
%                 count_modulated_pos_modulated_pos = count_modulated_pos_modulated_pos + 1;
% 
%             case 'modulated_neg_modulated_neg'
%                 count_modulated_neg_modulated_neg = count_modulated_neg_modulated_neg + 1;
% 
%             case 'generated_neg_generated_neg'
%                 count_generated_neg_generated_neg = count_generated_neg_generated_neg + 1;
% 
%             case 'generated_pos_generated_pos'
%                 count_generated_pos_generated_pos = count_generated_pos_generated_pos + 1;
% 
%             case {'modulated_neg_generated_neg',  'generated_neg_modulated_neg'}
%                 count_modulated_neg_generated_neg = count_modulated_neg_generated_neg + 1;
% 
%             case {'modulated_pos_generated_pos',  'generated_pos_modulated_pos'}
%                 count_modulated_pos_generated_pos = count_modulated_pos_generated_pos + 1;
% 
%             case {'generated_pos_generated_neg',  'generated_neg_generated_pos'}
%                 count_generated_pos_generated_neg = count_generated_pos_generated_neg + 1;
% 
%             case {'modulated_pos_modulated_neg',  'modulated_neg_modulated_pos'}
%                 count_modulated_pos_modulated_neg = count_modulated_pos_modulated_neg + 1;
% 
%             case {'modulated_pos_generated_neg',  'generated_neg_modulated_pos'}
%                 count_modulated_pos_generated_neg = count_modulated_pos_generated_neg + 1;
% 
%             case {'generated_pos_modulated_neg',  'modulated_neg_generated_pos'}
%                 count_generated_pos_modulated_neg = count_generated_pos_modulated_neg + 1;
%         end
% 
% 



%         if plot_interaction == true
%             clf
%             figure(2); title('Low Covariance')
%             hold on
%             x_coords = cellfun(@(c) c(1), spiral_path1);
%             y_coords = cellfun(@(c) c(2), spiral_path1);
%             plot(x_coords, y_coords, 'k-');
%             plot(x_coords(1), y_coords(1),'o','MarkerFaceColor','b');
% 
% 
%             x_coords = cellfun(@(c) c(1), spiral_path2);
%             y_coords = cellfun(@(c) c(2), spiral_path2);
%             plot(x_coords, y_coords, 'k-');
%             plot(x_coords(1), y_coords(1),'o','MarkerFaceColor','r');
% 
% 
%             for parcellation_ID = 1:22
%                 parcellation_template_1par = parcellation_template;
%                 parcellation_template_1par(isnan(parcellation_template_1par)) = 0;
% 
%                 parcellation_template_1par(parcellation_template_1par~=parcellation_ID) = 0;
% 
%                 B = bwboundaries(parcellation_template_1par(:,:,1),'noholes');
%                 for k = 1:length(B)
%                     boundary = B{k};
%                     plot(boundary(:,2),boundary(:,1),'-','linewidth',2,'color',[0,0,0])
%                 end
%             end
% 
%             xlabel('X-coordinate');
%             ylabel('Y-coordinate');
% 
%             pause(1)
%         end
%     end
% 
% 
% 
% end
end
%%

