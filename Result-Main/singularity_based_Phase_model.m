% Add all subfolders of the current directory to the MATLAB path
addpath(genpath([pwd]));

% Set the main folder to the current directory
main_folder = pwd;

% Set the hemisphere and subject identifiers
hemisphere = 1;
subject = 2;

% Load spatiotemporal bandpass filtered fMRI signal file
foldername = [main_folder, '/Sample Data/Motor Task/Preprocessed Data'];

% Determine the filename based on the hemisphere and subject
if hemisphere == 1
    filename = ['Preprocessed_spatiotemporalbandpass_data_motor_task_LEFT_sub', num2str(subject), '.mat'];
elseif hemisphere == 2
    filename = ['Preprocessed_spatiotemporalbandpass_data_motor_task_RIGHT_sub', num2str(subject), '.mat'];
end

% Load the data file
load(filename);
sigBPass = permute(DataOut, [1, 2, 3, 4]); % Order the data correctly

% Compute the phase field
phaseSig = nan(size(sigBPass));
for irow = 1:size(sigBPass, 1)
    for icol = 1:size(sigBPass, 2)
        temp1 = sigBPass(irow, icol, :);
        phaseSig(irow, icol, :) = angle(hilbert(temp1(:)));
    end
end

% Compute the phase gradient (vector) field
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

% Load the appropriate parcellation template based on the hemisphere
if hemisphere == 1
    load('parcellation_template.mat');
elseif hemisphere == 2
    load('parcellation_template22_RightBrain_subject1-100.mat');
    parcellation_template = parcellation_template22_RightBrain_100sub(:, :, 1);
end

% Select a specific time point for analysis
time = 90;

% Detect spirals in the phase field at the specified time point
detected_spirals = singularity_tracking(-vPhaseX(:, :, time), -vPhaseY(:, :, time), phaseSig(:, :, time));
Frame = detected_spirals.Frame1;

% Extract the smoothed phase map at the specified time point
smooth_phase_map = phaseSig(:, :, time);

% Define kernels for computing gradients
kx = [-0.5 0 0.5; -1 0 1; -0.5 0 0.5]; % Kernel for x-gradient
ky = [0.5 1 0.5; 0 0 0; -0.5 -1 -0.5]; % Kernel for y-gradient

% Extract positive and negative vortex centroids
WCentroids_pos = Frame.Nodes.PeakVorticityLoc(Frame.Nodes.Rotation == "Positive", :);
WCentroids_nega = Frame.Nodes.PeakVorticityLoc(Frame.Nodes.Rotation == "Negative", :);

% Measure spiral rotation
[omega_pos, omega_nega] = measure_spiral_rotation(WCentroids_pos, WCentroids_nega, phaseSig, time);

% Set the number of iterations for optimization
iterations = 1; % This optimization is tedious above 10 iterations.

% Optimize the radial influence
[sigma_pos, sigma_nega, costHistory] = optimize_radial_influence(WCentroids_pos, WCentroids_nega, omega_pos, omega_nega, parcellation_template, smooth_phase_map, iterations);

% Generate the model phase using the optimized parameters
[~, ~, Model_Phase] = vortex_based_phase_vector_field_model(WCentroids_pos, WCentroids_nega, omega_pos, omega_nega, sigma_pos, sigma_nega, parcellation_template);

% Plot the model phase with centroids
figure;
imagesc(Model_Phase(1:175, 1:251)); % Specific data range
set(gca, 'YDir', 'normal');
colormap(hsv); % Set colormap to 'hsv'
hold on;
scatter(WCentroids_pos(:, 1), WCentroids_pos(:, 2), 'wo', 'LineWidth', 3); % White 'o' for positive centroids
scatter(WCentroids_nega(:, 1), WCentroids_nega(:, 2), 'ko', 'LineWidth', 3); % Black 'o' for negative centroids
hold off;
axis tight;
axis equal;
yticklabels([]);
xticklabels([]);

% Plot boundaries of the parcellation template
hold on;
for parcellation_ID = 1:22
    parcellation_template_1par = parcellation_template > 0;
    parcellation_template_1par(isnan(parcellation_template_1par)) = 0;
    parcellation_template_1par(parcellation_template_1par ~= parcellation_ID) = 0;
    B = bwboundaries(parcellation_template_1par(:, :, 1), 'noholes');
    for k = 1:length(B)
        boundary = B{k};
        plot(boundary(:, 2), boundary(:, 1), '-', 'LineWidth', 2, 'Color', [0, 0, 0]);
    end
end
col = colorbar;
col.Label.String = 'Phase (rad)';
col.Ticks = [-pi + 0.1, 0, pi - 0.1];
col.TickLabels = {'-\pi', '0', '\pi'};
set(gca, 'FontSize', 18);

% Plot the real phase map with centroids
figure;
imagesc(phaseSig(:, :, time)); % Specific time slice
set(gca, 'YDir', 'normal');
colormap(hsv); % Ensure consistent colormap
hold on;
scatter(WCentroids_pos(:, 1), WCentroids_pos(:, 2), 'wo', 'LineWidth', 3); % White 'o' for positive centroids
scatter(WCentroids_nega(:, 1), WCentroids_nega(:, 2), 'ko', 'LineWidth', 3); % Black 'o' for negative centroids
hold off;
axis tight;
axis equal;
yticklabels([]);
xticklabels([]);

% Plot boundaries of the parcellation template
hold on;
for parcellation_ID = 1:22
    parcellation_template_1par = parcellation_template > 0;
    parcellation_template_1par(isnan(parcellation_template_1par)) = 0;
    parcellation_template_1par(parcellation_template_1par ~= parcellation_ID) = 0;
    B = bwboundaries(parcellation_template_1par(:, :, 1), 'noholes');
    for k = 1:length(B)
        boundary = B{k};
        plot(boundary(:, 2), boundary(:, 1), '-', 'LineWidth', 2, 'Color', [0, 0, 0]);
    end
end
col = colorbar;
col.Label.String = 'Phase (rad)';
col.Ticks = [-pi + 0.1, 0, pi - 0.1];
col.TickLabels = {'-\pi', '0', '\pi'};
set(gca, 'FontSize', 18);
