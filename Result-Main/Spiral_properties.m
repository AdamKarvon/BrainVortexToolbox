restoredefaultpath
 cd '\\suphys.physics.usyd.edu.au\akar5239\BrainVortexToolbox-main'
addpath(genpath([pwd]))
main_folder = pwd;

% load task label (language task)
disp(['loading task label...'])

% load task label of each subject
foldername = [main_folder,'/Sample Data/Motor Task/TaskLabel'];
cd(foldername)
name = dir(pwd) ;
file_name2 = ['MotorTaskLabelAllSubject.mat'];
load (file_name2);
for isubject = 1:150
    fullTime_allsubject{isubject} = TaskLabel_AllSubject_motor{isubject};
end



RH1_gen_singularity = zeros(283,150);
LH1_gen_singularity = zeros(283,150);
RF1_gen_singularity = zeros(283,150);
LF1_gen_singularity = zeros(283,150);
T1_gen_singularity = zeros(283,150);
RH2_gen_singularity = zeros(283,150);
LH2_gen_singularity = zeros(283,150);
RF2_gen_singularity = zeros(283,150);
LF2_gen_singularity = zeros(283,150);
T2_gen_singularity = zeros(283,150);
non_task_gen_singularity = zeros(283,150);



for subject = 1:150
    instSpeed = [];
    lifetimes = [];

    for hemisphere = 2
    foldername = [main_folder,'/Sample Data/Motor Task/Analysis'];
    cd(foldername)
    if hemisphere == 1
        filename = ['spiral_detection_sub_left',num2str(subject),'.mat'];
    elseif hemisphere == 2
    filename = ['spiral_detection_sub_right',num2str(subject),'.mat'];
    end
    load(filename)
    detected_spirals = detected_spirals;


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




    for t = 1:283
        instSpeedTime = [];
        % Find graph at time frame
        fieldName = ['Frame' num2str(t)];
        G = detected_spirals.(fieldName);
        interaction = interaction_table.(fieldName);


    for nodeidx  = 1:numnodes(G)
        %Number of Singularities and when they were evoked
        generation_frame = G.Nodes.generationFrame(nodeidx);
        if ismember(generation_frame, start_end_time_right_hand(1,1):start_end_time_right_hand(1,2))
            RH1_gen_singularity(t,subject) =  RH1_gen_singularity(t,subject) + 1;
        elseif ismember(generation_frame, start_end_time_left_hand(1,1):start_end_time_left_hand(1,2))
            LH1_gen_singularity(t,subject) =  LH1_gen_singularity(t,subject) +1;
        elseif ismember(generation_frame, start_end_time_left_foot(1,1):start_end_time_left_foot(1,2))
            LF1_gen_singularity(t,subject) = LF1_gen_singularity(t,subject) + 1;
        elseif ismember(generation_frame, start_end_time_tongue(1,1):start_end_time_tongue(1,2))
            T1_gen_singularity(t,subject) = T1_gen_singularity(t,subject) +1;
        elseif ismember(generation_frame, start_end_time_right_foot(1,1):start_end_time_right_foot(1,2))
            RF1_gen_singularity(t,subject) = RF1_gen_singularity(t,subject) + 1;
        elseif ismember(generation_frame, start_end_time_tongue(2,1):start_end_time_tongue(2,2))
            T2_gen_singularity(t,subject) = T2_gen_singularity(t,subject) +1;
        elseif ismember(generation_frame, start_end_time_left_foot(2,1):start_end_time_left_foot(2,2))
            LF2_gen_singularity(t,subject) = LF2_gen_singularity(t,subject) + 1;
        elseif ismember(generation_frame, start_end_time_right_hand(2,1):start_end_time_right_hand(2,2))
            RH2_gen_singularity(t,subject) =  RH2_gen_singularity(t,subject) + 1;
        elseif ismember(generation_frame, start_end_time_left_hand(2,1):start_end_time_left_hand(2,2))
            LH2_gen_singularity(t,subject) =  LH2_gen_singularity(t,subject) +1;
        elseif ismember(generation_frame, start_end_time_right_foot(2,1):start_end_time_right_foot(2,2))
            RF2_gen_singularity(t,subject) = RF2_gen_singularity(t,subject) + 1;
        else
            non_task_gen_singularity(t,subject) = non_task_gen_singularity(t,subject) + 1;
        end



        % Distribution of instantaneous speeds

        instVelocity = G.Nodes.Velocity(nodeidx,:); 
        if ~isnan(instVelocity)
            instSpeed(end+1) = sqrt(instVelocity(1)^2 + instVelocity(2)^2);
            instSpeedTime(end+1) = sqrt(instVelocity(1)^2 + instVelocity(2)^2);
        end

        % Distribtuion of Spiral Lifetimes

        spiralDuration = G.Nodes.spiral_lifetime(nodeidx);
        if ~isnan(spiralDuration)
            lifetimes(end+1) = spiralDuration;
        end

    end %Node loop End

    averageSpeed(t,subject) = mean(instSpeedTime);
    end % time loop end

    lifetimeDistribution{subject,hemisphere} = lifetimes;
    speedDistribution{subject,hemisphere} = instSpeed;


    end % hemisphere loop end
end % subject loop end


%%

% PLOT NO. OF SINGULARITIES DIVIDED INTO THEIR GENRATION TIMES
figure(1)
    %bar([non_task_gen_singularity(:,1),LH1_gen_singularity(:,1),LF1_gen_singularity(:,1),RH1_gen_singularity(:,1),RF1_gen_singularity(:,1),T1_gen_singularity(:,1),LH2_gen_singularity(:,1),LF2_gen_singularity(:,1),RH2_gen_singularity(:,1),RF2_gen_singularity(:,1),T2_gen_singularity(:,1)],'stacked')
    b = bar([mean(RH1_gen_singularity(:,:),2),mean(LF1_gen_singularity(:,:),2),mean(T1_gen_singularity(:,:),2),mean(RF1_gen_singularity(:,:),2),mean(LH1_gen_singularity(:,:),2),mean(T2_gen_singularity(:,:),2),mean(LF2_gen_singularity(:,:),2),mean(RH2_gen_singularity(:,:),2),mean(LH2_gen_singularity(:,:),2),mean(RF2_gen_singularity(:,:),2),mean(non_task_gen_singularity(:,:),2)],'stacked');
% Define a list of unique colors for each stack
colors = [0 0 1; % Blue 
          1 0.5 0; % Orange
          0 1 1; % Cyan
          1 0 1; % Magenta
          1 0 0; %Red
          0 1 1; % Cyan
          1 0.5 0; % Orange   
          0 0 1; % Blue
          1 0 0; % Red
          1 0 1; %Magenta 
          0 0 0];% Black

% Assign colors to each stack
for i = 1:length(b)
    b(i).FaceColor = 'flat';
    b(i).CData = repmat(colors(i,:), size(b(i).CData,1), 1);
end
    leg = legend('Right Hand','Left Foot','Tongue','Right Foot','Left Hand','','','','','','Non-task evoked');
    leg.FontSize = 16;

    %title("Average Number of Vortices At each time frame")
    % xlabel('time')
    ylabel('Avergae number of Identified Singularites', 'FontSize',16)
    set(gca, 'FontSize', 16);  % Set font size of tick labels to 15
    xticksExtended = [start_end_time_right_hand(1,1), start_end_time_right_hand(1,2),...
                    start_end_time_left_foot(1,1), start_end_time_left_foot(1,2), ...
                    start_end_time_tongue(1,1), start_end_time_tongue(1,2),...
                    start_end_time_right_foot(1,1), start_end_time_right_foot(1,2),...
                    start_end_time_left_hand(1,1), start_end_time_left_hand(1,2)...
                  start_end_time_tongue(2,1),start_end_time_tongue(2,2),...
                  start_end_time_left_foot(2,1), start_end_time_left_foot(2,2),...
                  start_end_time_right_hand(2,1),start_end_time_right_hand(2,2),...
                  start_end_time_left_hand(2,1), start_end_time_left_hand(2,2), ...
                  start_end_time_right_foot(2,1), start_end_time_right_foot(2,2)];
    xticks(xticksExtended); % Include the task start and end times in the ticks
    tickLabelsExtended = {'RH Start','RH End',...
        'LF Start','LF End',...
        'T Start', 'T End',...
        'RF Start', 'RF End'...,
        'LH Start', 'LH End',...
        'T2 Start','T2 End',...
        'LF2 Start', 'LF2 End',...
        'RH2 Start','RH2 End',...
        'LH2 Start', 'LH2 End', ...
        'RF2 Start', 'RF2 End'};
    xticklabels(tickLabelsExtended)
    set(gca, 'FontSize', 16);  % Set font size of tick labels to 15
    xtickangle(90)


    %% Plot HIstograms of spiral lifetime and Normalized efect speed

figure(1)

% Define the number of subjects
numSubjects = 150;

% Concatenate all lifetime data into a single array
allLifetimes = [];
for i = 1:numSubjects
    allLifetimes = [allLifetimes, lifetimeDistribution{i, 2}];
end

% Create the histogram
[binCounts, binEdges] = histcounts(allLifetimes, 'Normalization', 'probability', 'BinWidth', 1);

% Convert bin edges to bin centers
binCenters = binEdges(1:end-1) + diff(binEdges) / 2;
binCentersSec = binCenters * 0.72; % Convert to seconds
% Calculate mean and SEM of lifetimes
meanLifetime = mean(allLifetimes);
semLifetime = std(allLifetimes) / sqrt(numSubjects);

% Create the line plot
subplot(1,2,1)
plot(binCentersSec, binCounts, 'b-', 'LineWidth', 3);
hold on;
fill([binCentersSec fliplr(binCentersSec)], [binCounts+SEM fliplr(binCounts-SEM)], 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

xlabel('Lifetime (s)', 'FontSize', 16);
ylabel('Probability density', 'FontSize', 16);
set(gca, 'FontSize', 15);  % Set font size of tick labels to 15

% Change the y-axis to a logarithmic scale
set(gca, 'YScale', 'log');
xtickFrames = 0:5:150; % Adjust the step based on your data range
xtickSeconds = xtickFrames * 0.72;
xticks(xtickFrames);
xticklabels(arrayfun(@(x) sprintf('%.1f', x), xtickSeconds, 'UniformOutput', false));

xlim([0,75])

% Concatenate all speed data into a single array
allspeeds = [];
for i = 1:numSubjects
    allspeeds = [allspeeds, speedDistribution{i, 2}];
end

% Create the histogram
[binCounts, binEdges] = histcounts(allspeeds, 'Normalization', 'probability', 'BinWidth', 1);

% Convert bin edges to bin centers
binCenters = binEdges(1:end-1) + diff(binEdges) / 2;

% Calculate mean and SEM of speeds
meanSpeed = mean(allspeeds);
semSpeed = std(allspeeds) / sqrt(numSubjects);


subplot(1,2,2)
% Create the line plot
plot(binCenters, binCounts, 'b-', 'LineWidth', 3);
hold on;
fill([binCenters fliplr(binCenters)], [binCounts+SEM fliplr(binCounts-SEM)], 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

xlabel('Propagation Speed (mm/s)', 'FontSize', 16);
ylabel('Probability density', 'FontSize', 16);
set(gca, 'FontSize', 15);  % Set font size of tick labels to 15

% Change the y-axis to a logarithmic scale
set(gca, 'YScale', 'log');

hold off;

disp(['Mean Lifetime: ', num2str(meanLifetime*0.72), ' s with SEM: ', num2str(semLifetime*0.72), ' s']);
disp(['Mean Propagation Speed: ', num2str(meanSpeed), ' mm/s with SEM: ', num2str(semSpeed), ' mm/s']);



    



