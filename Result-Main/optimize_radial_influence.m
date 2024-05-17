function [sigma_pos, sigma_nega, costHistory, variance_captured_x] = optimize_radial_influence(WCentroids_pos, WCentroids_nega, omega_pos, omega_nega, parcellation_template, smooth_phase_map, iterations)
% This function optimizes the radial influence (sigma) for vortex-based phase vector field model using gradient descent.
% Inputs:
% - WCentroids_pos: Centroids of counterclockwise (positive) vortices
% - WCentroids_nega: Centroids of clockwise (negative) vortices
% - omega_pos: Rotation angles for positive vortices
% - omega_nega: Rotation angles for negative vortices
% - parcellation_template: Template for region of interest
% - smooth_phase_map: The phase map used for cost computation
% Outputs:
% - sigma_pos: Optimized radial influence for positive vortices
% - sigma_nega: Optimized radial influence for negative vortices
% - costHistory: History of cost values during optimization

% Initialize sigma values
sigma_pos = 10 .* ones(size(WCentroids_pos, 1), 1);
sigma_nega = 10 .* ones(size(WCentroids_nega, 1), 1);

% Gradient descent settings
alpha = 10; % Learning rate
delta = 1; % Perturbation for finite difference approximation

% Initialize array to store cost at each iteration
costHistory = zeros(iterations, 1);

tic;
for i = 1:iterations
    % Calculate the current Model_Phase
    [~, ~, Model_Phase] = vortex_based_phase_vector_field_model(WCentroids_pos, WCentroids_nega, omega_pos, omega_nega, sigma_pos, sigma_nega, parcellation_template);
    current_cost = computeCost(smooth_phase_map, Model_Phase);

    % Initialize gradient vectors   
    grad_sigma_pos = zeros(size(sigma_pos));
    grad_sigma_nega = zeros(size(sigma_nega));

    % Store the current cost in the history array
    costHistory(i) = current_cost;

    % Approximate gradient for sigma_pos
    for j = 1:length(sigma_pos)
        sigma_pos_temp = sigma_pos;
        sigma_pos_temp(j) = sigma_pos_temp(j) - delta;
        [~, ~, Model_Phase_temp] = vortex_based_phase_vector_field_model(WCentroids_pos, WCentroids_nega, omega_pos, omega_nega, sigma_pos_temp, sigma_nega, parcellation_template);
        cost_temp = computeCost(smooth_phase_map, Model_Phase_temp);
        grad_sigma_pos(j) = (cost_temp - current_cost) / delta;
    end

    % Approximate gradient for sigma_nega
    for j = 1:length(sigma_nega)
        sigma_nega_temp = sigma_nega;
        sigma_nega_temp(j) = sigma_nega_temp(j) - delta;
        [~, ~, Model_Phase_temp] = vortex_based_phase_vector_field_model(WCentroids_pos, WCentroids_nega, omega_pos, omega_nega, sigma_pos, sigma_nega_temp, parcellation_template);
        cost_temp = computeCost(smooth_phase_map, Model_Phase_temp);
        grad_sigma_nega(j) = (cost_temp - current_cost) / delta;
    end

    % Update parameters
    sigma_pos = sigma_pos + alpha * grad_sigma_pos;
    sigma_nega = sigma_nega + alpha * grad_sigma_nega;
end

% Plot the cost function history to visualize convergence
figure;
plot(1:iterations, costHistory, '-o');
xlabel('Iteration');
ylabel('Cost');
set(gca, 'FontSize', 18);
axis tight;

variance_captured_x = 1 - computeCost(smooth_phase_map, Model_Phase);

end

function cost = computeCost(smooth_phase_map, Model_Phase)  

        Model_Phase= Model_Phase(1:175,1:251);

        a = smooth_phase_map - nanmean(smooth_phase_map);
        b = Model_Phase - nanmean(Model_Phase);

        r = nansum(nansum(a.*b))/sqrt(nansum(nansum(a.*a))*nansum(nansum(b.*b)));

        cost = 1-r.^2;
end
