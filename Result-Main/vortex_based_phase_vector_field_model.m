function [Model_vx, Model_vy, Model_phase] = vortex_based_phase_vector_field_model(WCentroids_pos, WCentroids_nega, omega_pos, omega_nega, sigma_pos, sigma_nega, parcellation_template)
% This function models the phase field based on vortex centroids
% Inputs:
% - WCentroids_pos: Centroids of counterclockwise (positive) vortices
% - WCentroids_nega: Centroids of clockwise (negative) vortices
% - omega_pos: Angular offset for positive vortices
% - omega_nega: Angular offset for negative vortices
% - sigma_pos: Gaussian scaling factor for positive vortices
% - sigma_nega: Gaussian scaling factor for negative vortices
% - parcellation_template: Template for region of interest

% Initialize grid and phase matrix
x = linspace(1, 251, 251);
y = linspace(1, 251, 251);
[X, Y] = meshgrid(x, y);
S = zeros(251, 251);

% Impose boundary conditions
% Create padding and concatenate with the template
padding = zeros(76, 251);
parcellation_template_padded = [parcellation_template; padding];

% Filter the parcellation template and resize
parcellation_template_filt = parcellation_template_padded ./ parcellation_template_padded;
parcellation_template_filt = logical(~isnan(parcellation_template_filt));
parcellation_template_filt = imresize(parcellation_template_filt, [251 251], 'nearest');

% Set region of interest to 1
S(parcellation_template_filt) = 1;

% Find boundaries of the region
B = bwboundaries(parcellation_template_filt(:, :), 'noholes');
boundaryY = B{1}(:, 1);
boundaryX = B{1}(:, 2);

% Determine points inside the boundary
[in, on] = inpolygon(X, Y, boundaryX, boundaryY);
idx = in | on;

% Set points outside the boundary to NaN
X(~idx) = NaN;
Y(~idx) = NaN;

% Initialize accumulators for vector components
X_component = zeros(size(S));
Y_component = zeros(size(S));

% Process clockwise (negative) vortices
for CW_spiral = 1:size(WCentroids_nega, 1)
    if isnan(WCentroids_nega(CW_spiral))
        continue;
    else
        % Define spiral location and charge
        x0 = WCentroids_nega(CW_spiral, 1);
        y0 = WCentroids_nega(CW_spiral, 2);
        m = -1;

        % Calculate radial distance and relative coordinates
        r = sqrt((X - x0).^2 + (Y - y0).^2);
        Y_rel = Y - y0;
        X_rel = X - x0;

        % Calculate angle and Gaussian scaling factor
        theta = mod(atan2(Y_rel, X_rel) - omega_nega(CW_spiral), 2 * pi) - pi;
        gaussian_scaling = exp(-(r.^2) / (2 * sigma_nega(CW_spiral)^2));

        % Calculate phase change and update vector components
        phase_change = m * theta;
        X_component = X_component + gaussian_scaling .* cos(phase_change);
        Y_component = Y_component + gaussian_scaling .* sin(phase_change);
    end
end

% Process counterclockwise (positive) vortices
for ACW_spiral = 1:size(WCentroids_pos, 1)
    if isnan(WCentroids_pos(ACW_spiral))
        continue;
    else
        % Define spiral location and charge
        x0 = WCentroids_pos(ACW_spiral, 1);
        y0 = WCentroids_pos(ACW_spiral, 2);
        m = 1;

        % Calculate radial distance and relative coordinates
        r = sqrt((X - x0).^2 + (Y - y0).^2);
        Y_rel = Y - y0;
        X_rel = X - x0;

        % Calculate angle and Gaussian scaling factor
        theta = mod(atan2(Y_rel, X_rel) - omega_pos(ACW_spiral), 2 * pi) - pi;
        gaussian_scaling = exp(-(r.^2) / (2 * sigma_pos(ACW_spiral)^2));

        % Calculate phase change and update vector components
        phase_change = m * theta;
        X_component = X_component + gaussian_scaling .* cos(phase_change);
        Y_component = Y_component + gaussian_scaling .* sin(phase_change);
    end
end

% Calculate phase field from vector components
S_phase = atan2(Y_component, X_component);
S = exp(1i * S_phase);
Model_phase = angle(S);

% Compute the spatial gradient of the phase
PhaseX = zeros(size(Model_phase));
PhaseY = zeros(size(Model_phase));
for iTime = 1
    for iX = 1:size(Model_phase, 1)
        PhaseX(iX, 2:end-1, iTime) = (anglesubtract(Model_phase(iX, 3:end, iTime), Model_phase(iX, 1:end-2, iTime))) / 2;
    end
    for iY = 1:size(Model_phase, 2)
        PhaseY(2:end-1, iY, iTime) = (anglesubtract(Model_phase(3:end, iY, iTime), Model_phase(1:end-2, iY, iTime))) / 2;
    end
end

% Standardize vector length
PhaseX = PhaseX ./ sqrt(PhaseX.^2 + PhaseY.^2);
PhaseY = PhaseY ./ sqrt(PhaseX.^2 + PhaseY.^2);

% Final output of the modeled phase vector field
Model_vx = PhaseX(1:175, 1:251);
Model_vy = PhaseY(1:175, 1:251);

end

