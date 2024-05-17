function [modulated_spiral_data, generated_spiral_data,modulated_clusters, generated_clusters] = task_clustering_centres(movement_block, task, task_frame,min_clusters,display_clusters,hemisphere) 
cd '\\suphys.physics.usyd.edu.au\akar5239\BrainVortexToolbox-main'
addpath(genpath([pwd]))
main_folder = pwd;

if hemisphere == 1
    load('parcellation_template.mat')
elseif hemisphere ==2
    load('parcellation_template22_RightBrain_subject1-100.mat')
    parcellation_template = parcellation_template22_RightBrain_100sub(:,:,1);
end
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
% name = dir(pwd) ;
task_times = TaskLabel_AllSubject_motor{1,1};
disp(['loading task data...'])
tic
if movement_block == 1
    switch task
        case 'RH'
            task_onset = task_times(1,1);
            task_end = task_times(1,2);
            if hemisphere == 1
                load ('Task_specific_spiralData_leftHEM.mat', 'RH_block1');
            elseif hemisphere == 2
                load ('Task_specific_spiralData_rightHEM.mat', 'RH_block1');
            end
            task_data = RH_block1;
        case 'LF'
            task_onset = task_times(2,1);
            task_end = task_times(2,2);
            if hemisphere == 1
                load ('Task_specific_spiralData_leftHEM.mat', 'LF_block1');
            elseif hemisphere == 2
                load ('Task_specific_spiralData_rightHEM.mat', 'LF_block1');
            end
            task_data = LF_block1;
        case 'T'
            task_onset = task_times(3,1);
            task_end = task_times(3,2);
            if hemisphere == 1
                load ('Task_specific_spiralData_leftHEM.mat', 'T_block1');
            elseif hemisphere == 2
                load ('Task_specific_spiralData_rightHEM.mat', 'T_block1');
            end
            task_data = T_block1;
        case 'RF'
            task_onset = task_times(4,1);
            task_end = task_times(4,2);
            if hemisphere == 1
                load ('Task_specific_spiralData_leftHEM.mat', 'RF_block1');
            elseif hemisphere == 2
                load ('Task_specific_spiralData_rightHEM.mat', 'RF_block1');
            end
            task_data = RF_block1;
        case 'LH'
            task_onset = task_times(5,1);
            task_end = task_times(5,2);
            if hemisphere == 1
                load ('Task_specific_spiralData_leftHEM.mat', 'LH_block1');
            elseif hemisphere == 2
                load ('Task_specific_spiralData_rightHEM.mat', 'LH_block1');
            end
            task_data = LH_block1;
    end

elseif movement_block == 2
    switch task

        case 'T'
            task_onset = task_times(6,1);
            task_end = task_times(6,2);
            if hemisphere == 1
                load ('Task_specific_spiralData_leftHEM.mat', 'T_block2');
            elseif hemisphere == 2
                load ('Task_specific_spiralData_rightHEM.mat', 'T_block2');
            end
            task_data = T_block2;
        case 'LF'
            task_onset = task_times(7,1);
            task_end = task_times(7,2);
            if hemisphere == 1
                load ('Task_specific_spiralData_leftHEM.mat', 'LF_block2');
            elseif hemisphere == 2
                load ('Task_specific_spiralData_rightHEM.mat', 'LF_block2');
            end
            task_data = LF_block2;
        case 'RH'
            task_onset = task_times(8,1);
            task_end = task_times(8,2);
            if hemisphere == 1
                load ('Task_specific_spiralData_leftHEM.mat', 'RH_block2');
            elseif hemisphere == 2
                load ('Task_specific_spiralData_rightHEM.mat', 'RH_block2');
            end
            task_data = RH_block2;
        case 'LH'
            task_onset = task_times(9,1);
            task_end = task_times(9,2);
            if hemisphere == 1
                load ('Task_specific_spiralData_leftHEM.mat', 'LH_block2');
            elseif hemisphere == 2
                load ('Task_specific_spiralData_rightHEM.mat', 'LH_block2');
            end
            task_data = LH_block2;
        case 'RF'
            task_onset = task_times(10,1);
            task_end = task_times(10,2);
            if hemisphere == 1
                load ('Task_specific_spiralData_leftHEM.mat', 'RF_block2');
            elseif hemisphere == 2
                load ('Task_specific_spiralData_rightHEM.mat', 'RF_block2');
            end
            task_data = RF_block2;
    end

end
toc

for time = 1:283
    for subject = 1:150
        if size(task_data.modulated.pos{subject},2) > 1
            pos_modulated_pattID = task_data.modulated.pos{subject}(:,1);
            pos_modulated_centres = task_data.modulated.pos{subject}(:,2);
            pos_modulated_time = task_data.modulated.pos{subject}(:,3);
            pos_modulated_velocity = task_data.modulated.pos{subject}(:,4);
        end
        if size(task_data.modulated.neg{subject},2) > 1
            neg_modulated_pattID = task_data.modulated.neg{subject}(:,1);
            neg_modulated_centres = task_data.modulated.neg{subject}(:,2);
            neg_modulated_time = task_data.modulated.neg{subject}(:,3);
            neg_modulated_velocity = task_data.modulated.neg{subject}(:,4);
        end
        if size(task_data.generated.pos{subject},2) > 1
            pos_generated_pattID = task_data.generated.pos{subject}(:,1);
            pos_generated_centres = task_data.generated.pos{subject}(:,2);
            pos_generated_time = task_data.generated.pos{subject}(:,3);
            pos_generated_velocity = task_data.generated.pos{subject}(:,4);
        end
        if size(task_data.generated.neg{subject},2) > 1
            neg_generated_pattID = task_data.generated.neg{subject}(:,1);
            neg_generated_centres = task_data.generated.neg{subject}(:,2);
            neg_generated_time = task_data.generated.neg{subject}(:,3);
            neg_generated_velocity = task_data.generated.neg{subject}(:,4);

        end

        for i = 1:size(pos_modulated_time)
            t = find(pos_modulated_time{i} == time);
            if isempty(t)
                pos_modulated{subject,i,time} = [];
            else
                if t > size(pos_modulated_velocity{i})
                    pos_modulated_velocity{i}(:,t) = [NaN];
                end
                vel = {[pos_modulated_velocity{i}(1,t) pos_modulated_velocity{i}(2,t)]};
                pos_modulated{subject,i,time} = cell2mat([pos_modulated_pattID{i} pos_modulated_centres{i}(t) vel]);
            end
        end

        for i = 1:size(pos_generated_time)
            t = find(pos_generated_time{i} == time);
            if isempty(t)
                pos_generated{subject,i,time} = [];
            else
                if t > size(pos_generated_velocity{i})
                    pos_generated_velocity{i}(:,t) = [NaN];
                end

                vel = {[pos_generated_velocity{i}(1,t) pos_generated_velocity{i}(2,t)]};

                pos_generated{subject,i,time} = cell2mat([pos_generated_pattID{i} pos_generated_centres{i}(t) vel]);
            end
        end

        for i = 1:size(neg_modulated_time)
            t = find(neg_modulated_time{i} == time);
            if isempty(t)
                neg_modulated{subject,i,time} = [];
            else
                if t > size(neg_modulated_velocity{i})
                    neg_modulated_velocity{i}(:,t) = [NaN];
                end

                vel = {[neg_modulated_velocity{i}(1,t) neg_modulated_velocity{i}(2,t)]};
                neg_modulated{subject,i,time} = cell2mat([neg_modulated_pattID{i} neg_modulated_centres{i}(t) vel]);
            end
        end

        for i = 1:size(neg_generated_time)
            t = find(neg_generated_time{i} == time);
            if isempty(t)
                neg_generated{subject,i,time} = [];
            else
                if t > size(neg_generated_velocity{i})
                    neg_generated_velocity{i}(:,t) = [NaN];
                end
                vel = {[neg_generated_velocity{i}(1,t) neg_generated_velocity{i}(2,t)]};
                neg_generated{subject,i,time} = cell2mat([neg_generated_pattID{i} neg_generated_centres{i}(t) vel]);
            end
        end

    end % Subject Loop end    
end %Time loop end

nSubjects = size(pos_modulated, 1);
nCoordinates = size(pos_modulated, 2);
nTime = size(pos_modulated, 3);


for time = 1:nTime
    % Initialize arrays to store the converted data
    subjectID = [];
    coordinateID = [];
    xCoordinates = [];
    yCoordinates = [];
    spiralVelocityX = [];
    spiralVelocityY = [];
    % Loop through the cell array and extract the data
    for subject = 1:nSubjects
        for coordinate = 1:nCoordinates
            % Extract subject, coordinate ID, x, and y values
            % Check if the cell is empty
            if ~isempty(pos_modulated{subject, coordinate, time})
                % Extract subject, coordinate ID, x, and y values
                subjectID = [subjectID; subject];
                coordinateID = [coordinateID; pos_modulated{subject, coordinate, time}(1)];
                xCoordinates = [xCoordinates; pos_modulated{subject, coordinate, time}(2)];
                yCoordinates = [yCoordinates; pos_modulated{subject, coordinate, time}(3)];
                spiralVelocityX = [spiralVelocityX; pos_modulated{subject, coordinate, time}(4)];
                spiralVelocityY = [spiralVelocityY; pos_modulated{subject, coordinate, time}(5)];
            end
        end
    end
    % Create the final double array by concatenating the columns
    if ~isempty([subjectID, coordinateID, xCoordinates, yCoordinates, spiralVelocityX, spiralVelocityY])
        task_pos_modulated{time} = [subjectID, coordinateID, xCoordinates, yCoordinates, spiralVelocityX, spiralVelocityY];
    end
end


nSubjects = size(pos_generated, 1);
nCoordinates = size(pos_generated, 2);
nTime = size(pos_generated, 3);


for time = 1:nTime
    % Initialize arrays to store the converted data
    subjectID = [];
    coordinateID = [];
    xCoordinates = [];
    yCoordinates = [];
    spiralVelocityX =[];
    spiralVelocityY = [];
    % Loop through the cell array and extract the data
    for subject = 1:nSubjects
        for coordinate = 1:nCoordinates
            % Extract subject, coordinate ID, x, and y values
            % Check if the cell is empty
            if ~isempty(pos_generated{subject, coordinate, time})
                % Extract subject, coordinate ID, x, and y values
                subjectID = [subjectID; subject];
                coordinateID = [coordinateID; pos_generated{subject, coordinate, time}(1)];
                xCoordinates = [xCoordinates; pos_generated{subject, coordinate, time}(2)];
                yCoordinates = [yCoordinates; pos_generated{subject, coordinate, time}(3)];
                spiralVelocityX = [spiralVelocityX; pos_generated{subject, coordinate, time}(4)];
                spiralVelocityY = [spiralVelocityY; pos_generated{subject, coordinate, time}(5)];
            end
        end
    end
    % Create the final double array by concatenating the columns
    if ~isempty( [subjectID, coordinateID, xCoordinates, yCoordinates, spiralVelocityX, spiralVelocityY])
        task_pos_generated{time} =  [subjectID, coordinateID, xCoordinates, yCoordinates, spiralVelocityX, spiralVelocityY];
    end
end

nSubjects = size(neg_modulated, 1);
nCoordinates = size(neg_modulated, 2);
nTime = size(neg_modulated, 3);


for time = 1:nTime
    % Initialize arrays to store the converted data
    subjectID = [];
    coordinateID = [];
    xCoordinates = [];
    yCoordinates = [];
    spiralVelocityX =[];
    spiralVelocityY = [];
    % Loop through the cell array and extract the data
    for subject = 1:nSubjects
        for coordinate = 1:nCoordinates
            % Extract subject, coordinate ID, x, and y values
            % Check if the cell is empty
            if ~isempty(neg_modulated{subject, coordinate, time})
                % Extract subject, coordinate ID, x, and y values
                subjectID = [subjectID; subject];
                coordinateID = [coordinateID; neg_modulated{subject, coordinate, time}(1)];
                xCoordinates = [xCoordinates; neg_modulated{subject, coordinate, time}(2)];
                yCoordinates = [yCoordinates; neg_modulated{subject, coordinate, time}(3)];
                spiralVelocityX = [spiralVelocityX; neg_modulated{subject, coordinate, time}(4)];
                spiralVelocityY = [spiralVelocityY; neg_modulated{subject, coordinate, time}(5)];
            end
        end
    end
    % Create the final double array by concatenating the columns
    if ~isempty( [subjectID, coordinateID, xCoordinates, yCoordinates, spiralVelocityX, spiralVelocityY])
        task_neg_modulated{time} =  [subjectID, coordinateID, xCoordinates, yCoordinates, spiralVelocityX, spiralVelocityY];
    end
end

nSubjects = size(neg_generated, 1);
nCoordinates = size(neg_generated, 2);
nTime = size(neg_generated, 3);


for time = 1:nTime
    % Initialize arrays to store the converted data
    subjectID = [];
    coordinateID = [];
    xCoordinates = [];
    yCoordinates = [];
    spiralSize = [];
    spiralVelocityX =[];
    spiralVelocityY = [];
    % Loop through the cell array and extract the data
    for subject = 1:nSubjects
        for coordinate = 1:nCoordinates
            % Extract subject, coordinate ID, x, and y values
            % Check if the cell is empty
            if ~isempty(neg_generated{subject, coordinate, time})
                % Extract subject, coordinate ID, x, and y values
                subjectID = [subjectID; subject];
                coordinateID = [coordinateID; neg_generated{subject, coordinate, time}(1)];
                xCoordinates = [xCoordinates; neg_generated{subject, coordinate, time}(2)];
                yCoordinates = [yCoordinates; neg_generated{subject, coordinate, time}(3)];
                spiralVelocityX = [spiralVelocityX; neg_generated{subject, coordinate, time}(4)];
                spiralVelocityY = [spiralVelocityY; neg_generated{subject, coordinate, time}(5)];
            end
        end
    end
    % Create the final double array by concatenating the columns
    if ~isempty( [subjectID, coordinateID, xCoordinates, yCoordinates, spiralVelocityX, spiralVelocityY])
        task_neg_generated{time} =  [subjectID, coordinateID, xCoordinates, yCoordinates, spiralVelocityX, spiralVelocityY];
    end
end


%%
for l = task_onset+task_frame
      
    eps = 5;       % Search radius around core points
    for i = 1:2
        for iNeighbours = 100:-1:15
            min_no_neigbours = iNeighbours;  % cluster must be surrounded by min_neighbours to be considered a cluster
            
            if i == 1
                
                Z = task_neg_generated{l}(:,3:4);
                distances = influenceDistance(Z);
                [neg_generated_clusters, ~] = dbscan(distances,eps,min_no_neigbours,'Distance','precomputed');   % HyperParamaters neighbourhood search radius, minimum number of neighbors for core point
                if max(neg_generated_clusters) >= min_clusters 
                    disp(min_no_neigbours)
                    break
                end

            elseif i == 2
                
                Z = task_pos_generated{l}(:,3:4);
                distances = influenceDistance(Z);
                [pos_generated_clusters, ~] = dbscan(distances,eps,min_no_neigbours,'Distance','precomputed');   % HyperParamaters neighbourhood search radius, minimum number of neighbors for core point

                if max(pos_generated_clusters) >= min_clusters
                    disp(min_no_neigbours)
                    break
                end
            end

        end
    end

    for i = 1:2
        for iNeighbours = 100:-1:15
            min_no_neigbours = iNeighbours;  % cluster must be surrounded by min_neighbours to be considered a cluster
            
            if i == 1
                Z = task_neg_modulated{l}(:,3:4);
                distances = influenceDistance(Z);
                [neg_modulated_clusters, ~] = dbscan(distances,eps,min_no_neigbours,'Distance','precomputed');   % HyperParamaters neighbourhood search radius, minimum number of neighbors for core point
                if max(neg_modulated_clusters) >= min_clusters 
                    disp(min_no_neigbours)
                    break
                end

            elseif i == 2
                Z = task_pos_modulated{l}(:,3:4);
                distances = influenceDistance(Z);
                [pos_modulated_clusters, ~] = dbscan(distances,eps,min_no_neigbours,'Distance','precomputed');   % HyperParamaters neighbourhood search radius, minimum number of neighbors for core point

                if max(pos_modulated_clusters) >= min_clusters
                    disp(min_no_neigbours)
                    break
                end
            end

        end
    end
   

    %% Cluster Analysis
    modulated_clusters = struct;
    modulated_clusters.pos = cell(max(pos_modulated_clusters),1);
    modulated_clusters.neg = cell(max(neg_modulated_clusters),1);
    clusters = transpose(unique(pos_modulated_clusters));    
    for c = clusters
        mask = pos_modulated_clusters == c;
        if c == -1
            modulated_clusters.posNoise = task_pos_modulated{l}(mask,:);
        else
            modulated_clusters.pos{c} = task_pos_modulated{l}(mask,:);
        end
    end
    clusters = transpose(unique(neg_modulated_clusters));
    for c = clusters
        mask = neg_modulated_clusters == c;
        if c == -1
            modulated_clusters.negNoise = task_neg_modulated{l}(mask,:);
        else
            modulated_clusters.neg{c} = task_neg_modulated{l}(mask,:);
        end
    end


    generated_clusters = struct;
    generated_clusters.pos = cell(max(pos_generated_clusters),1);
    generated_clusters.neg = cell(max(neg_generated_clusters),1);
    clusters = transpose(unique(pos_generated_clusters));
    for c = clusters
        mask = pos_generated_clusters == c;
        if c == -1
            generated_clusters.posNoise = task_pos_generated{l}(mask,:);
        else
            generated_clusters.pos{c} = task_pos_generated{l}(mask,:);
        end
    end
    clusters = transpose(unique(neg_generated_clusters));
    for c = clusters
        mask = neg_generated_clusters == c;
        if c == -1
            generated_clusters.negNoise = task_neg_generated{l}(mask,:);
        else
            generated_clusters.neg{c} = task_neg_generated{l}(mask,:);
        end
    end

% Generated Spiral data

generated_spiral_data = struct;

% AntiClockwise and Clockwise clustered sprial data
for cluster = 1:size(generated_clusters.pos,1)
    cluster_data = cell(size(generated_clusters.pos{cluster},1),5);
    for clustered_spiral = 1:size(generated_clusters.pos{cluster})
        subject = generated_clusters.pos{cluster}(clustered_spiral,1);
        pattID = generated_clusters.pos{cluster}(clustered_spiral,2);
        spiral_ID = cell2mat(task_data.generated.pos{subject}(:,1)) == pattID;
        BUGFIX = task_data.generated.pos{subject}(spiral_ID,:);               % Duplicate PattIDs occur when a spiral disappears and then reappears
        cluster_data(clustered_spiral,:) = [subject BUGFIX(1,:)];
    end
    generated_spiral_data.pos{cluster} = cluster_data;
end
for cluster = 1:size(generated_clusters.neg,1)
    cluster_data = cell(size(generated_clusters.neg{cluster},1),5);
    for clustered_spiral = 1:size(generated_clusters.neg{cluster})
        subject = generated_clusters.neg{cluster}(clustered_spiral,1);
        pattID = generated_clusters.neg{cluster}(clustered_spiral,2);
        spiral_ID = cell2mat(task_data.generated.neg{subject}(:,1)) == pattID;
        BUGFIX = task_data.generated.neg{subject}(spiral_ID,:);      
        cluster_data(clustered_spiral,:) = [subject BUGFIX(1,:)];
    end
    generated_spiral_data.neg{cluster} = cluster_data;
end

% Unclustered spiral data

%Anticlockwise
unclustered_data = cell(size(generated_clusters.posNoise,1),5);
for unclustered_spiral = 1:size(generated_clusters.posNoise,1)
    subject = generated_clusters.posNoise(unclustered_spiral,1);
    pattID = generated_clusters.posNoise(unclustered_spiral,2);
    spiral_ID = cell2mat(task_data.generated.pos{subject}(:,1)) == pattID;
    BUGFIX = task_data.generated.pos{subject}(spiral_ID,:);
    unclustered_data(unclustered_spiral,:) = [subject BUGFIX(1,:)];
end
generated_spiral_data.posNoise = unclustered_data;

%CLockwise
unclustered_data = cell(size(generated_clusters.negNoise,1),5);
for unclustered_spiral = 1:size(generated_clusters.negNoise,1)
    subject = generated_clusters.negNoise(unclustered_spiral,1);
    pattID = generated_clusters.negNoise(unclustered_spiral,2);
    spiral_ID = cell2mat(task_data.generated.neg{subject}(:,1)) == pattID;
    BUGFIX = task_data.generated.neg{subject}(spiral_ID,:);
    unclustered_data(unclustered_spiral,:) = [subject BUGFIX(1,:)];
end
generated_spiral_data.negNoise = unclustered_data;


% Modulated Spiral Data 
modulated_spiral_data = struct;

for cluster = 1:size(modulated_clusters.pos,1)
    cluster_data = cell(size(modulated_clusters.pos{cluster},1),5);
    for clustered_spiral = 1:size(modulated_clusters.pos{cluster})
        subject = modulated_clusters.pos{cluster}(clustered_spiral,1);    % Retireve subject no. and subject pattID from clustered data
        pattID = modulated_clusters.pos{cluster}(clustered_spiral,2);
        spiral_ID = cell2mat(task_data.modulated.pos{subject}(:,1)) == pattID;  % Find the pattern within subjects unprocessed data
        cluster_data(clustered_spiral,:) = [subject task_data.modulated.pos{subject}(spiral_ID,:)]; % store the full the spiral infromation 
    end
    modulated_spiral_data.pos{cluster} = cluster_data;
end
for cluster = 1:size(modulated_clusters.neg,1)
    cluster_data = cell(size(modulated_clusters.neg{cluster},1),5);
    for clustered_spiral = 1:size(modulated_clusters.neg{cluster})
        subject = modulated_clusters.neg{cluster}(clustered_spiral,1);
        pattID = modulated_clusters.neg{cluster}(clustered_spiral,2);
        spiral_ID = cell2mat(task_data.modulated.neg{subject}(:,1)) == pattID;
        cluster_data(clustered_spiral,:) = [subject task_data.modulated.neg{subject}(spiral_ID,:)];
    end
    modulated_spiral_data.neg{cluster} = cluster_data;
end

% Unclustered Spiral Data
%Anticlockwise
unclustered_data = cell(size(modulated_clusters.posNoise,1),5);
for unclustered_spiral = 1:size(modulated_clusters.posNoise,1)
    subject = modulated_clusters.posNoise(unclustered_spiral,1);
    pattID = modulated_clusters.posNoise(unclustered_spiral,2);
    spiral_ID = cell2mat(task_data.modulated.pos{subject}(:,1)) == pattID;
    BUGFIX = task_data.modulated.pos{subject}(spiral_ID,:);
    unclustered_data(unclustered_spiral,:) = [subject BUGFIX(1,:)];
end
modulated_spiral_data.posNoise = unclustered_data;

%clockwise
unclustered_data = cell(size(modulated_clusters.negNoise,1),5);
for unclustered_spiral = 1:size(modulated_clusters.negNoise,1)
    subject = modulated_clusters.negNoise(unclustered_spiral,1);
    pattID = modulated_clusters.negNoise(unclustered_spiral,2);
    spiral_ID = cell2mat(task_data.modulated.neg{subject}(:,1)) == pattID;
    BUGFIX = task_data.modulated.neg{subject}(spiral_ID,:);
    unclustered_data(unclustered_spiral,:) = [subject BUGFIX(1,:)];
end
modulated_spiral_data.negNoise = unclustered_data;


    %% Order clusters within struct to be descending in size
    [~,id] = sort(cellfun(@(c) size(c,1), modulated_clusters.pos ),'descend');
    modulated_clusters.pos = modulated_clusters.pos(id);

    [~,id] = sort(cellfun(@(c) size(c,1), modulated_clusters.neg ),'descend');
    modulated_clusters.neg = modulated_clusters.neg(id);


    [~,id] = sort(cellfun(@(c) size(c,1), generated_clusters.pos ),'descend');
    generated_clusters.pos = generated_clusters.pos(id);

    [~,id] = sort(cellfun(@(c) size(c,1), generated_clusters.neg ),'descend');
    generated_clusters.neg = generated_clusters.neg(id);

    [~,id] = sort(cellfun(@(c) size(c,1), modulated_spiral_data.pos ),'descend');
    modulated_spiral_data.pos = modulated_spiral_data.pos(id);

    [~,id] = sort(cellfun(@(c) size(c,1), modulated_spiral_data.neg ),'descend');
    modulated_spiral_data.neg = modulated_spiral_data.neg(id);


    [~,id] = sort(cellfun(@(c) size(c,1), generated_spiral_data.pos ),'descend');
    generated_spiral_data.pos = generated_spiral_data.pos(id);

    [~,id] = sort(cellfun(@(c) size(c,1), generated_spiral_data.neg ),'descend');
    generated_spiral_data.neg = generated_spiral_data.neg(id);






 

end

end

function resultVector= influenceDistance(Z)
numPoints = size(Z, 1);
distanceMatrix = zeros(numPoints);

for i = 1:numPoints
    for j = i+1:numPoints
        d = norm(Z(i, :) - Z(j, :)); % Euclidean distance between points i and j

        distanceMatrix(i, j) = d;

        % Since the distance matrix is symmetric, we can copy the result to the other side
        distanceMatrix(j, i) = distanceMatrix(i, j);
    end
end
resultVector = squareform(distanceMatrix);
end

