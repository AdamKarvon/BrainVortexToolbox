% add all subfolders of the current directory to the MATLAB path
addpath(genpath([pwd]));
main_folder = pwd;

% Set the number of subjects and task flag
No_of_Subject = 1; % Number of subjects used for analysis, randomly selected from HCP database (S1200)


% Task-specific phase gradient field for Motor Task
disp('loading task label...');

% Load task label for each subject for Motor Task
foldername = [main_folder, '/Sample Data/Motor Task/TaskLabel'];
cd(foldername);
name = dir(pwd);
file_name2 = 'MotorTaskLabelAllSubject.mat';
load(file_name2);
for isubject = 1:150
    fullTime_allsubject{isubject} = TaskLabel_AllSubject_motor{isubject};
end


cd(main_folder);
load('parcellation_template22_RightBrain_subject1-100.mat');

% Loop over each hemisphere (1: left, 2: right)
for hemisphere = 1

    % Initialize complex vector field for all subjects
    complex_vec_field_ALLsubject = nan(175, 251, 284, 150);
    for subject = 1:No_of_Subject
        % Load spatiotemporal bandpass filtered fMRI signal file
        foldername = [main_folder, '/Sample Data/Motor Task/Preprocessed Data'];
        if hemisphere == 1
            filename = ['Preprocessed_spatiotemporalbandpass_data_motor_task_LEFT_sub', num2str(subject), '.mat'];
        elseif hemisphere == 2
            filename = ['Preprocessed_spatiotemporalbandpass_data_motor_task_RIGHT_sub', num2str(subject), '.mat'];
        end
        load(filename);
        sigBPass = permute(DataOut, [1, 2, 3, 4]);

        % Compute the phase field
        phaseSig = nan(size(sigBPass));
        for irow = 1:size(sigBPass, 1)
            for icol = 1:size(sigBPass, 2)
                temp1 = sigBPass(irow, icol, :);
                phaseSig(irow, icol, :) = angle(hilbert(temp1(:)));
            end
        end

        % Compute the phase gradient (vector) field
        Vx_flowmap_norm_phase = [];
        Vy_flowmap_norm_phase = [];
        vPhaseX = zeros(size(phaseSig));
        vPhaseY = zeros(size(phaseSig));
        cd(main_folder);
        for iTime = 1:size(phaseSig, 3)
            for iX = 1:size(phaseSig, 1)
                vPhaseX(iX, 2:end-1, iTime) = (anglesubtract(phaseSig(iX, 3:end, iTime), phaseSig(iX, 1:end-2, iTime))) / 2;
            end
            for iY = 1:size(phaseSig, 2)
                vPhaseY(2:end-1, iY, iTime) = (anglesubtract(phaseSig(3:end, iY, iTime), phaseSig(1:end-2, iY, iTime))) / 2;
            end
        end
        complex_vec_field(:, :, :) = (-vPhaseX(:, :, :)) + 1i .* (-vPhaseY(:, :, :));

        % Clear temporary variables
        clearvars vPhaseX vPhaseY phaseSig;

        % Store the task-specific phase vector field for each subject
        complex_vec_field_ALLsubject(:, :, :, subject) = complex_vec_field;
    end

    % Compute the average complex vector field for all subjects in the hemisphere
    if hemisphere == 1
        left_hem_complex_vec_field_avg = nanmean(complex_vec_field_ALLsubject, 4);
    elseif hemisphere == 2
        right_hem_complex_vec_field_avg = nanmean(complex_vec_field_ALLsubject, 4);
        parcellation_template_right = parcellation_template22_RightBrain_100sub(:, :, 1);
    end
end

% Plotting the results
% Task IDs for different limb movements
task_ID_right_limb = [1, 8, 4, 10];
task_ID_left_limb = [5, 9, 2, 7];

figure(1);
pre_end = 3;

% Loop over the number of subplots to create
for i = 1:2
    if i == 1
        % Compute the vector fields for right hand and right foot tasks
        RH_task_vec_fields = (left_hem_complex_vec_field_avg(:, :, fullTime_allsubject{1, 1}(task_ID_right_limb(1), 2) - pre_end) + left_hem_complex_vec_field_avg(:, :, fullTime_allsubject{1, 1}(task_ID_right_limb(2), 2) - pre_end)) / 2;
        RF_task_vec_fields = (left_hem_complex_vec_field_avg(:, :, fullTime_allsubject{1, 1}(task_ID_right_limb(3), 2) - pre_end) + left_hem_complex_vec_field_avg(:, :, fullTime_allsubject{1, 1}(task_ID_right_limb(4), 2) - pre_end)) / 2;

        directions_left_hem(:) = [squeeze(angle(complex_vec_field_ALLsubject(117, 128, fullTime_allsubject{1, 1}(task_ID_right_limb(1), 2) - pre_end, :))); squeeze(angle(complex_vec_field_ALLsubject(117, 128, fullTime_allsubject{1, 1}(task_ID_right_limb(2), 2) - pre_end, :)))];
        complex_vec_field_avg = RH_task_vec_fields;
    elseif i == 2
        % Compute the vector fields for left hand and left foot tasks
        % Uncomment the following lines if needed
        % LH_task_vec_fields = (right_hem_complex_vec_field_avg(:, :, fullTime_allsubject{1, 1}(task_ID_left_limb(1), 2) - pre_end) + right_hem_complex_vec_field_avg(:, :, fullTime_allsubject{1, 1}(task_ID_left_limb(2), 2) - pre_end)) / 2;
        % LF_task_vec_fields = (right_hem_complex_vec_field_avg(:, :, fullTime_allsubject{1, 1}(task_ID_left_limb(3), 2) - pre_end) + right_hem_complex_vec_field_avg(:, :, fullTime_allsubject{1, 1}(task_ID_left_limb(4), 2) - pre_end)) / 2;
        % complex_vec_field_avg = LH_task_vec_fields;

        % directions_right_hem(:) = [squeeze(angle(complex_vec_field_ALLsubject(130, 136, fullTime_allsubject{1, 1}(task_ID_left_limb(1), 2) - pre_end, :))); squeeze(angle(complex_vec_field_ALLsubject(130, 136, fullTime_allsubject{1, 1}(task_ID_left_limb(2), 2) - pre_end, :)))];
    end

    subplot(1, 2, i); % Create subplot

    % Calculate direction (angle) and magnitude of the vector field
    directions = angle(complex_vec_field_avg(:, :)); % Angle in radians
    magnitudes = abs(complex_vec_field_avg(:, :)); % Magnitude

    meanMagnitude = nanmean(magnitudes(:));
    stdMagnitude = nanstd(magnitudes(:));
    zScores(:, :) = (magnitudes - meanMagnitude) ./ stdMagnitude;

    if i == 1
        hold on;
        imagesc(zScores);
        h = streamslice([1:251], [1:175], cos(directions), sin(directions), 6, 'w');
        set(h, 'Color', [0 0 0]);
        set(h, 'LineWidth', 1.5);
    end

    if i == 2
        p = polarhistogram(directions_left_hem, 12);
        title('Propagation direction distribution', 'FontSize', 16);
        ax = gca;
        ax.FontSize = 16;
        p.FaceColor = [0.5, 0.5, 0.5];
        p.LineWidth = 1;
    else
        for parcellation_ID = [6, 8]
            parcellation_template_1par = parcellation_template;
            parcellation_template_1par(isnan(parcellation_template_1par)) = 0;
            parcellation_template_1par(parcellation_template_1par ~= parcellation_ID) = 0;
            B = bwboundaries(parcellation_template_1par, 'noholes');
            for k = 1:length(B)
                boundary = B{k};
                plot(boundary(:, 2), boundary(:, 1), '-', 'LineWidth', 2, 'Color', [1, 1, 1]);
            end
        end
    end

    switch i
        case 1
            xlabel('X Coordinate');
            xlim([105, 165]);
            ylim([65, 160]);
            ylabel('Y Coordinate');
            scatter(128, 117, 1000, "pentagram", 'filled', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'LineWidth', 2);
            % Coordinates for scale bar placement
            x_start = 110; % X-coordinate for the start of the scale bar
            y_position = 68; % Y-coordinate for both start and end of the scale bar
            x_end = x_start + 25; % 5 pixels from x_start
            line([x_start, x_end], [y_position, y_position], 'Color', 'w', 'LineWidth', 2);
            text(x_start, y_position - 2, '50 mm', 'Color', 'w', 'FontSize', 15, 'VerticalAlignment', 'top');
            xticks([]);
            yticks([]);
            xticklabels([]);
            yticklabels([]);
            c = colorbar('FontSize', 15); % Create the colorbar
            c.Label.String = 'Vector Magnitude Z-score'; % Set the label text
    end

    % Turn off hold for the current subplot
    hold off;
    clim([0, 6]);
    colormap(summer(6));
    clear zScores;
end




