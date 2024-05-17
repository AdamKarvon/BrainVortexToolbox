restoredefaultpath
cd '\\suphys.physics.usyd.edu.au\akar5239\BrainVortexToolbox-main'
addpath(genpath([pwd]))
main_folder = pwd;

%%

    image1 = zeros(175,251);

    image2 = zeros(175,251);

for hemisphere =1:2% 1 for left hemisphere, 2 for right hemisphere, keep this as 1 for demonstration purpose

% flagTask:
% 1 = language task, original 100 subjects;
% 2 = language task, additional 100 subjects;
% 3 = working memory task;
% 4 = Motor Task
flagTask = 4;



% load task label (language task)
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

if hemisphere == 1
    load('parcellation_template.mat')
    load('Task_specific_spiralDistribution_leftHEM.mat')
    task_ID= [1,8];
elseif hemisphere ==2
    load('parcellation_template22_RightBrain_subject1-100.mat')
    parcellation_template_right = parcellation_template22_RightBrain_100sub(:,:,1) ;
    load('Task_specific_spiralDistribution_rightHEM.mat')
    task_ID= [5,9];
end



for task_ID = task_ID
    task_start = fullTime_allsubject{1,1}(task_ID,1);
    task_end = fullTime_allsubject{1,1}(task_ID,2);
    
    switch task_ID
        case 1
            distribution = RH_block1_subAvg_dist;
            task = 'Right Finger Movement Block 1';
        case 2
            distribution = LF_block1_subAvg_dist;
            task = 'Left Toes Movement Block 1';
        case 3
            distribution = T_block1_subAvg_dist;
            task = 'Tongue Movement Block 1';
        case 4
            distribution = RF_block1_subAvg_dist;
            task = 'Right Toes Movement Block 1';
        case 5
            distribution = LH_block1_subAvg_dist;
            task = 'Left Finger Movement Block 1';
        case 6
            distribution = T_block2_subAvg_dist;
            task = 'Tongue Movement Block 2';
        case 7
            distribution = LF_block2_subAvg_dist;
            task = 'Left Toes Movement Block 2';
        case 8
            distribution = RH_block2_subAvg_dist;
            task = 'Right Finger Movement Block 2';
        case 9
            distribution = LH_block2_subAvg_dist;
            task = 'Left Finger Movement Block 2';
        case 10
            distribution = RF_block2_subAvg_dist;
            task = 'Right Toes Movement Block 2';
    end


    time = task_end-3;
    

    if hemisphere == 1
        image1 = image1 + (distribution.modulated.neg(:,:,time) + distribution.modulated.pos(:,:,time) + distribution.generated.neg(:,:,time) + distribution.generated.pos(:,:,time));

        clear distribution

    elseif hemisphere == 2
        image2 = image2 + (distribution.modulated.neg(:,:,time)+ distribution.modulated.pos(:,:,time)+ distribution.generated.neg(:,:,time) + distribution.generated.pos(:,:,time));


        clear distribution
    end
end
end


%%

% First subplot for Pre-exisiting
%subplot(1,2,1); % This selects the second grid cell
figure;
imagesc(image1./300);
title('Right Hand','FontSize',16);
hold on;
for parcellation_ID = 1:22
    parcellation_template_1par = parcellation_template>0;
    parcellation_template_1par(isnan(parcellation_template_1par)) = 0;
    parcellation_template_1par(parcellation_template_1par~=parcellation_ID) = 0;
    B = bwboundaries(parcellation_template_1par(:,:,1), 'noholes');
    for k = 1:length(B)
        boundary = B{k};
        plot(boundary(:,2), boundary(:,1), '-', 'linewidth', 2, 'color', [0,0,0]);
    end
end
set(gca, "ydir", 'normal');
c= colorbar('Fontsize',20);
xticklabels([""]);
yticklabels([""]);
colormap jet;
axis equal
axis tight

% Second subplot for Task Evoked
figure;
imagesc(image2./300);
title('Left Hand','FontSize',16);
hold on 
for parcellation_ID = 1:22
    parcellation_template_1par = parcellation_template_right>0;
    parcellation_template_1par(isnan(parcellation_template_1par)) = 0;
    parcellation_template_1par(parcellation_template_1par~=parcellation_ID) = 0;
    B = bwboundaries(parcellation_template_1par(:,:,1), 'noholes');
    for k = 1:length(B)
        boundary = B{k};
        plot(boundary(:,2), boundary(:,1), '-', 'linewidth', 2, 'color', [0,0,0]);
    end
end
set(gca, "ydir", 'normal');
c = colorbar('Fontsize',20);
colormap jet;
axis equal
axis tight
xticklabels([""]);
yticklabels([""]);
c.Label.String = 'Trial Averaged Spiral Density';

% Coordinates for scale bar placement
x_start = 200; % X-coordinate for the start of the scale bar
y_position = 20; % Y-coordinate for both start and end of the scale bar
x_end = x_start + 25; % 5 pixels from x_start

% Draw the scale bar line
line([x_start, x_end], [y_position, y_position], 'Color', 'k', 'LineWidth', 4);

% Add text label for the scale bar
% Adjust the position as needed to not overlap with the scale bar itself
text(x_start, y_position - 5, '50 mm', 'Color', 'k', 'FontSize', 15, 'VerticalAlignment', 'top');


