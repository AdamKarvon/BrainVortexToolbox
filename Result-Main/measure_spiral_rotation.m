function [omega_pos, omega_nega] = measure_spiral_rotation(WCentroids_pos, WCentroids_nega, phaseSig, time)
% This function measures the rotation of spiral waves in the phase data.
% Inputs:
% - WCentroids_pos: Centroids of counterclockwise (positive) vortices
% - WCentroids_nega: Centroids of clockwise (negative) vortices
% - phaseSig: 3D matrix of phase signals (x, y, time)
% - time: The time slice to analyze
% Outputs:
% - omega_pos: Rotation angles for positive vortices
% - omega_nega: Rotation angles for negative vortices

% Initialize rotation angle arrays
omega_nega = zeros(size(WCentroids_nega, 1), 1);
omega_pos = zeros(size(WCentroids_pos, 1), 1);

% Edge detection kernels
kx = [-1, 0, 1; -2, 0, 2; -1, 0, 1]; % Sobel kernel for X direction
ky = [-1, -2, -1; 0, 0, 0; 1, 2, 1]; % Sobel kernel for Y direction

% Identify rotation of real phase data
EdgeX = conv2(phaseSig(:,:,time), kx, 'same');
EdgeY = conv2(phaseSig(:,:,time), ky, 'same');
waveFront = sqrt(EdgeX.^2 + EdgeY.^2);
wavefront_mask = waveFront > 2*pi;

% Find the coordinates of points to overplot
[y, x] = find(wavefront_mask);

% Measure rotation for positive vortices
desiredRadius = 3;
for idx = 1:size(WCentroids_pos, 1)
    center = WCentroids_pos(idx, :);
    
    % Convert points to polar coordinates relative to the current center
    % and filter points outside the desired radius
    [theta, r, xWithinRadius, yWithinRadius] = convertToPolarRelative(center(1), center(2), x, y, desiredRadius);

    % Check if there are enough points within the radius
    if length(theta) <= 2
        WCentroids_pos(idx, :) = nan(1, 2);
        continue;
    else
        desiredRadius = 5;
        [theta, r, xWithinRadius, yWithinRadius] = convertToPolarRelative(center(1), center(2), x, y, desiredRadius);
        
        % Calculate median angle
        sinSum = sum(sin(theta));
        cosSum = sum(cos(theta));
        medianTheta = atan2(sinSum, cosSum);

        omega_pos(idx) = medianTheta;
    end
end

% Measure rotation for negative vortices
desiredRadius = 3;
for idx = 1:size(WCentroids_nega, 1)
    center = WCentroids_nega(idx, :);
    
    % Convert points to polar coordinates relative to the current center
    % and filter points outside the desired radius
    [theta, r, xWithinRadius, yWithinRadius] = convertToPolarRelative(center(1), center(2), x, y, desiredRadius);

    % Check if there are enough points within the radius
    if length(theta) <= 2
        WCentroids_nega(idx, :) = nan(1, 2);
        continue;
    else
        desiredRadius = 5;
        [theta, r, xWithinRadius, yWithinRadius] = convertToPolarRelative(center(1), center(2), x, y, desiredRadius);

        % Calculate median angle
        sinSum = sum(sin(theta));
        cosSum = sum(cos(theta));
        medianTheta = atan2(sinSum, cosSum);

        omega_nega(idx) = medianTheta;
    end
end

end

function [theta, r, xWithinRadius, yWithinRadius] = convertToPolarRelative(centerX, centerY, x, y, desiredRadius)
% Convert (x, y) points to polar coordinates relative to a center (centerX, centerY)
% and filter points within a desired radius
% Outputs:
% - theta: Angles of points within the radius
% - r: Radii of points within the radius
% - xWithinRadius: x-coordinates of points within the radius
% - yWithinRadius: y-coordinates of points within the radius

% Calculate relative coordinates
xRel = x - centerX;
yRel = y - centerY;

% Convert to polar coordinates
theta = atan2(yRel, xRel);
r = sqrt(xRel.^2 + yRel.^2);

% Filter points within the desired radius
withinRadius = r <= desiredRadius;
theta = theta(withinRadius);
r = r(withinRadius);
xWithinRadius = x(withinRadius);
yWithinRadius = y(withinRadius);
end
