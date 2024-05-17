function [outerTable,interaction_table]= pattDetection_v5(vPhaseX,vPhaseY, phaseSig)
%% Function Description
% The function uses vorticity contours to identify vortex patches, labeling them based on generation time, rotation, and size.
N = size(vPhaseX,3); %number of time points
dt = 0.72;
addpath(genpath([pwd]))
main_folder = pwd;

% Computes the Q-criterion from the 2D cartesian velocity
% field Vx, Vy and the grid spacing grid_h

Vx_flowmap_norm_phase = -vPhaseX./sqrt(vPhaseX.^2+vPhaseY.^2) ; % standardize vector length as 1
Vy_flowmap_norm_phase = -vPhaseY./sqrt(vPhaseX.^2+vPhaseY.^2) ; % standardize vector length as 1


% Vorticity values for vortex identification
vorticity = [];

for time = 1:size(vPhaseX,3)
    temp1_vx = Vx_flowmap_norm_phase(:,:,time);
    temp1_vy = Vy_flowmap_norm_phase(:,:,time);
    [x,y] = meshgrid(1:size(temp1_vx,2),1:size(temp1_vx,1));
    [vorticity(:,:,time), ~]  = curl(x,y,temp1_vx,temp1_vy);  % curl value of phase vector field
end


% Prepare output variable: A table with a new column for each timestep
% containing the spiral distribution as a graph within the first row of the
% columns

% Convert numeric range to string array of valid variable names
varNames = "Frame" + (1:N); % Prefixing with 'Var' to ensure names are valid identifiers

% Create an empty cell array to hold the cell arrays for each variable
cellArrays = cell(1, N);

% Create a table by converting cell arrays to a table directly
outerTable = cell2table(cellArrays, 'VariableNames', varNames);
interaction_table = cell2table(cellArrays, 'VariableNames', varNames);

% Inititialize graph to recursively update
main = graph;

% Initialize unique node identifier
nodeID = 0;


% Loop over each time step to detect and analyze vortices.
for t = 1:N%:timeSteps
    edge_contour_list = cell(1);
    % vortex_contour_list = cell(1);

    % Extract the slice of matrix corresponding to the current time step.

    flow_magnitude = sqrt(vPhaseX(:,:,time).^2+vPhaseY(:,:,time).^2);

    edge_contours = contourc(flow_magnitude, [0.25 0.25]);

    % Initialize indices for tracking contours and vertices within each vortex.
    contour_index = 1;
    edge_index = 0;
    vertex_index = 1;

    % Iterate through the contour data to identify the relevant contours 
    for i = 1:size(edge_contours,2)
        if i == contour_index
            contour_level = edge_contours(1,contour_index);
            if contour_level == 0.25 % Check if we are at the correct contour level
                contour_no_points = edge_contours(2,contour_index); % Number of points in this contour
                vertex_index = 1; % Reset vertex index for a new vortex
                edge_index = edge_index + 1; % Move to the next vortex
                contour_index = contour_index + 1 + contour_no_points; % Skip to the next contour start
            else
                disp('ERROR: Misaligned contour index') % Error handling for unexpected contour data
            end
        else
            % Store the contour vertices for each vortex
            edge_contour_list{edge_index}(vertex_index,:) = [edge_contours(1,i) , edge_contours(2,i)];
            vertex_index = vertex_index + 1; % Move to the next vertex within the current vortex
        end
    end


    % Determine Node Contours

    % Extract the slice of matrix corresponding to the current time step.
    currentSlice = vorticity(:,:,t);

    vorticity_filt(:,:) = currentSlice > 1 | currentSlice < -1; 

    % BEGIN SPIRAL IDENTIFICATION
    % -- Selecting regions with vorticity contours as spirals
    % -- CC will contain all pixels within the contours for Frame i


    CC = bwconncomp(vorticity_filt,8) ;
    S = regionprops(CC,vorticity_filt,{'Centroid','WeightedCentroid'} );

    B1 = regionprops(CC,'BoundingBox');        % bounding box of patterns
    boundPatt = cat(1, B1.BoundingBox);

    % CONSTRUCT GRAPH NETWORK
    current_timestep = graph;

    % NODES
    % Add Nodes For each identified vorticity contour
    for iPatt = 1: size(CC.PixelIdxList,2)
        velocity = [nan, nan];
        spiral_lifetime = nan;
        edge_contour = cell(1);

        potentialNode = CC.PixelIdxList{iPatt};

        % Create a mask
        maskImage = false(175,251);
        maskImage(potentialNode) = true;

        % Isolate Vorticity within the vorticity Contour
        nodeVorticity = vorticity(:,:,t);
        nodeVorticity(~maskImage) = NaN;
       

        % Identify the Peak Vorticity  (Note: this max may be positive or negative)
        [~, peak_vorticity_loc] = max(abs(nodeVorticity),[],'all','linear');
        peak_vorticity = nodeVorticity(peak_vorticity_loc);

        [row, col]  = ind2sub([175 251],peak_vorticity_loc);

        peak_vorticity_loc = S(iPatt).WeightedCentroid; % x and y

        % Identfiy the Rotation of the Vortex
        if peak_vorticity > 0
            rotation = "Positive";
        elseif peak_vorticity < 0
            rotation = "Negative";
        end


        % Determine Whether Multiple spirals are contained within the vorticity contour
        if sum(nodeVorticity > 0 , 'all') > 0 && sum( nodeVorticity < 0 , 'all') > 0
            interacting = true;
        else
            interacting = false;
        end

        bounding_box{1} = boundPatt(iPatt,:);
        
        if interacting == true
            continue
        else
            nodeID = nodeID + 1;
            nodeProperties = table(nodeID, t ,rotation, peak_vorticity_loc, velocity, spiral_lifetime,bounding_box, interacting, edge_contour,...
                'VariableNames', ["NodeID","generationFrame" ,"Rotation", "PeakVorticityLoc", "Velocity" ,"spiral_lifetime","BoundingBox", "Interacting","EdgeContour"]);


            current_timestep = addnode(current_timestep, nodeProperties);
        end
    end

    % EDGES
    

    % for contour_index = 1:size(edge_contour_list,2)
    % 
    %     edge_contour = edge_contour_list{contour_index};
    % 
    %     node_list = [];
    %     for node_index = 1:size(current_timestep.Nodes,1)
    % 
    %         centre = current_timestep.Nodes.PeakVorticityLoc(node_index,:);
    % 
    %         % Circle parameters
    %         xc = centre(1); % x-coordinate of the circle's center
    %         yc = centre(2); % y-coordinate of the circle's center
    %         r = 5; % Radius of the circle
    % 
    %         % Number of points along the circumference
    %         numPoints = 10;
    % 
    %         % Angles for generating points along the circle's circumference
    %         theta = linspace(0, 2*pi, numPoints);
    % 
    %         % Coordinates of points along the circumference
    %         verticesX = xc + r * cos(theta);
    %         verticesY = yc + r * sin(theta);
    % 
    %         within_contour = inpolygon(verticesX,verticesY,edge_contour(:,1),edge_contour(:,2));
    % 
    %         if sum(within_contour,'all')>0
    %             current_timestep.Nodes.EdgeContour(node_index) = {edge_contour};
    % 
    %             node_list = [node_list current_timestep.Nodes.NodeID(node_index)];
    % 
    %         end
    %     end
    % 
    %     if length(node_list)>1
    %         pairs = nchoosek(node_list, 2); % Generate all combinations of two nodes
    % 
    % 
    %         for i = 1:size(pairs, 1)
    % 
    %             node1_index = find(current_timestep.Nodes.NodeID == pairs(i, 1));
    %             node2_index = find(current_timestep.Nodes.NodeID == pairs(i, 2));
    % 
    % 
    %             node1Loc = current_timestep.Nodes.PeakVorticityLoc(node1_index,:);
    %             node2Loc = current_timestep.Nodes.PeakVorticityLoc(node2_index,:);
    %             weight = norm(node1Loc - node2Loc);
    % 
    %             if (current_timestep.Nodes.Rotation(node1_index) ~= current_timestep.Nodes.Rotation(node2_index) && weight <= 30)
    %                 current_timestep = addedge(current_timestep, node1_index, node2_index, weight);
    % 
    %             end
    %         end
    %     end
    % end


    %current_timestep = minspantree(current_timestep,'Type','forest');
    


    
    % Assuming current_timestep.Nodes.PeakVorticityLoc contains the node locations
    nodeLocations = current_timestep.Nodes.PeakVorticityLoc;

    % Number of nodes
    numNodes = size(nodeLocations, 1);

    % Initialize a matrix to store distances between each pair of nodes
    distances = zeros(numNodes, numNodes);

    % Calculate distances between all pairs of nodes
    for i = 1:numNodes
        for j = 1:numNodes
            if i ~= j && (current_timestep.Nodes.Rotation(i) ~= current_timestep.Nodes.Rotation(j))
                distances(i,j) = sqrt((nodeLocations(i,1) - nodeLocations(j,1))^2 + (nodeLocations(i,2) - nodeLocations(j,2))^2);
            else
                distances(i,j) = inf; % To ensure a node doesn't consider itself as a neighbour
            end
        end
    end

    % For each node, find the indices of the 4 closest neighbours
    closestNeighbours = zeros(numNodes, 4);
    for i = 1:numNodes
        [weight, sortedIndices] = sort(distances(i,:));
        closestNeighbours(i,:) = sortedIndices(1:4);



        current_timestep = addedge(current_timestep, i, closestNeighbours(i,1), weight(1));
        current_timestep = addedge(current_timestep, i, closestNeighbours(i,2), weight(2));
        current_timestep = addedge(current_timestep, i, closestNeighbours(i,3), weight(3));
        current_timestep = addedge(current_timestep, i, closestNeighbours(i,4), weight(4));
    end
    % 
    current_timestep = simplify(current_timestep);

    % update main graph

    if t == 1
        main = current_timestep;
        interactions = [];
    else

        % Extract peak vorticity locations for the current and previous timesteps
        currentPeakLocs = current_timestep.Nodes.PeakVorticityLoc;
        previousPeakLocs = main.Nodes.PeakVorticityLoc;

        % Extract rotation data for the current and previous timesteps
        currentRotations = current_timestep.Nodes.Rotation;
        previousRotations = main.Nodes.Rotation;

        % Calculate Euclidean distances between nodes' peak vorticity locations across timesteps


        distances = zeros(size(previousPeakLocs,1),size(currentPeakLocs,1));


        for prevNode_index = 1:size(previousPeakLocs,1)
            for currNode_index = 1:size(currentPeakLocs,1)
                if currentRotations(currNode_index) ~= previousRotations(prevNode_index)
                    distances(prevNode_index,currNode_index) = inf;
                else
                    distances(prevNode_index,currNode_index) = sqrt((currentPeakLocs(currNode_index,1) - previousPeakLocs(prevNode_index,1))^2 + (currentPeakLocs(currNode_index,2) - previousPeakLocs(prevNode_index,2))^2);
                end
            end
        end

        % distances = [];
        % 
        % for currNode_index = 1:size(currentPeakLocs,1)
        %     for prevNode_index = 1:size(previousPeakLocs,1)
        %         distances(currNode_index,prevNode_index) =  sqrt((currentPeakLocs(currNode_index,1)-previousPeakLocs(prevNode_index,1)).^2 + (currentPeakLocs(currNode_index,2)-previousPeakLocs(prevNode_index,2)).^2);
        %         if currentRotations(currNode_index) ~= previousRotations(prevNode_index)
        %             distances(currNode_index,prevNode_index) = inf;
        %         end
        %     end
        % end

        costofnonassignment = {5*ones(numnodes(main),1),5*ones(numnodes(current_timestep),1)};


        [assignments,annihilatedNodes,generatedNodes] = assignkbest(distances,costofnonassignment);


        % These describe the list of what nodes got assigned where,
        % may be used to construct a full graph of the nodes over time
        finalPreviousNodes = [assignments{1}(:,1)];
        finalCurrentNodes = [assignments{1}(:,2)];


        % Assuming assignments have been made
        uniqueAssignments = unique(assignments{1}(:,2), 'stable');
        if length(uniqueAssignments) < size(assignments{1}, 1)
            error('Duplicate NodeID assignments detected. Review matching criteria.');
        else
            % Proceed with NodeID assignment
            for i = 1:size(assignments{1}, 1)
                prevNodeID = main.Nodes.NodeID(assignments{1}(i, 1));
                currentNodeIndex = assignments{1}(i, 2);
                current_timestep.Nodes.NodeID(currentNodeIndex) = prevNodeID;
                current_timestep.Nodes.generationFrame(currentNodeIndex) = main.Nodes.generationFrame(assignments{1}(i, 1));

            end
        end

        for annihilated = 1:size(annihilatedNodes{1},1)
            spiral_lifetime = t - main.Nodes.generationFrame(annihilatedNodes{1}(annihilated));

            fieldName = ['Frame' num2str(t-1)];
            outerTable.(fieldName).Nodes.spiral_lifetime(annihilatedNodes{1}(annihilated))=spiral_lifetime;
        end

        pair = ismember(annihilatedNodes{1},  main.Edges.EndNodes);


        interactingNodeIndex = annihilatedNodes{1}(pair);
        if ~isempty(interactingNodeIndex) && ~isempty(annihilatedNodes{1})

            for i =  1:size(interactingNodeIndex,1)
                [row,~] = find(main.Edges.EndNodes == interactingNodeIndex(i,1));
                [~,idx] = min(main.Edges.Weight(row));

                interactingNodeIndex(i,2) = row(idx);

            end

            annihilated_pair_rows = unique(interactingNodeIndex(:,2));

            interaction_type = cell(1);
            interaction_location = cell(1);
            interactingNodeIDs = cell(1);
            for irow = 1:size(annihilated_pair_rows,1)

                interacting_pair = main.Edges.EndNodes(annihilated_pair_rows(irow),:);

                if sum(interactingNodeIndex(:,2) == annihilated_pair_rows(irow)) == 2
                    if main.Nodes.Rotation(interacting_pair(1)) ~= main.Nodes.Rotation(interacting_pair(2))
                        interaction_type{1}(end+1,1) = "FA";
                    elseif main.Nodes.Rotation(interacting_pair(1)) == main.Nodes.Rotation(interacting_pair(2))
                        interaction_type{1}(end+1,1) = "PA";
                    end
                elseif sum(interactingNodeIndex(:,2) == annihilated_pair_rows(irow)) == 1
                    if any(ismember(main.Nodes.NodeID(interacting_pair),current_timestep.Nodes.NodeID)) % 1 Node continues
                        if main.Nodes.Rotation(interacting_pair(1)) ~= main.Nodes.Rotation(interacting_pair(2))
                            interaction_type{1}(end+1,1) = "PM";
                        elseif main.Nodes.Rotation(interacting_pair(1)) == main.Nodes.Rotation(interacting_pair(2))
                            interaction_type{1}(end+1,1) = "FM";
                        end
                    end
                end


                interaction_location{1}(end+1,:) = ( main.Nodes.PeakVorticityLoc(interacting_pair(1),:) + main.Nodes.PeakVorticityLoc(interacting_pair(2),:) )./2;

                interactingNodeIDs{1}(end+1,:) = main.Nodes.NodeID(interacting_pair);


            end

            interactions = table(interaction_location', interaction_type', interactingNodeIDs', ...
                'VariableNames', {'interaction_location', 'interaction_type', 'interacting_NodeIDs'});

        elseif isempty(interactingNodeIndex) && ~isempty(annihilatedNodes{1})

            interactions = [num2str(size(annihilatedNodes{1},1)), ' unexplained annihilations at frame ', num2str(t)];

        end





        % After Unique identification of continued vortices we may compare previous
        % locations to estimate velocity


        for inode = 1:current_timestep.numnodes
            continued_node = find(main.Nodes.NodeID == current_timestep.Nodes.NodeID(inode));

            if ~isempty(continued_node)
                dx = current_timestep.Nodes.PeakVorticityLoc(inode,1) - main.Nodes.PeakVorticityLoc(continued_node,1);
                dy = current_timestep.Nodes.PeakVorticityLoc(inode,2) - main.Nodes.PeakVorticityLoc(continued_node,2);

                vx = 2*dx/dt;   % 2mm cifti ord
                vy = 2*dy/dt;

                velocity = [vx,vy];


                current_timestep.Nodes.Velocity(inode,:) = velocity;
            end
        end



        main = current_timestep;

    end

    fieldName = ['Frame' num2str(t)];
    outerTable.(fieldName) = main;

    interaction_table.(fieldName)=  {interactions};
end
end


function Q = computeQCriterion2D(Vx_flowmap_norm_phase, Vy_flowmap_norm_phase, grid_h)
    % Computes the Q-criterion from the 2D cartesian velocity 
    % field Vx, Vy and the grid spacing grid_h

    % Calculate spatial derivatives of velocity components
    [dxVx, dyVx] = gradient(Vx_flowmap_norm_phase, grid_h);
    [dxVy, dyVy] = gradient(Vy_flowmap_norm_phase, grid_h);

    % Components of the rate of strain tensor for 2D
    S_xx = dxVx; % Normal strain in x-direction
    S_yy = dyVy; % Normal strain in y-direction
    S_xy = 0.5 * (dyVx + dxVy); % Shear strain

    % Rotation rate tensor (vorticity for 2D) component
    Omega_z = 0.5 * (dyVx - dxVy); % Only one component in 2D, representing out-of-plane rotation

    % Compute Q-criterion for 2D
    % In 2D, the Q-criterion simplifies to half the difference between the square of vorticity 
    % (rotation rate) and the square of the strain rate components.
    Q = 0.5 * (Omega_z.^2 - (S_xx.^2 + S_yy.^2 + 2*S_xy.^2));

    % Note: The Q-criterion definition used here is adapted for 2D flow analysis, focusing on the balance between
    % rotational and strain energies in the plane of the flow.
end



