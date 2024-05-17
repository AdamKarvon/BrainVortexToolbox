restoredefaultpath
 cd '\\suphys.physics.usyd.edu.au\akar5239\BrainVortexToolbox-main'
addpath(genpath([pwd]))
main_folder = pwd;


interacting_spirals_all = zeros(150,283);
plotting = 1;

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


Num_spirals = zeros(283,150);
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


GeneratedGenerated = zeros(283,150);
GeneratedModulated = zeros(283,150);
ModulatedModulated = zeros(283,150);




for subject = 1
    no_interacting_spirals = [];

    nodesToRemove = [];

    % load spatiotemporal bandpass filtered fMRI signal file
    foldername = [main_folder,'/Sample Data/Motor Task/Preprocessed Data'];
    %cd(foldername)
    if hemisphere == 1            % spatio
        filename = ['Preprocessed_spatiotemporalbandpass_data_motor_task_LEFT_sub',num2str(subject),'.mat'];
    elseif hemisphere == 2
        filename = ['Preprocessed_spatiotemporalbandpass_data_motor_task_RIGHT_sub',num2str(subject),'.mat'];
    end
    load(filename)

    sigBPass = permute(DataOut,[1,2,3,4]); % Order the data correctly

    % phase field
    phaseSig = nan(size(sigBPass));
    for irow = 1:size(sigBPass,1)
        for icol = 1:size(sigBPass,2)
            temp1 = sigBPass(irow,icol,:);
            phaseSig(irow,icol,:) = angle(hilbert(temp1(:)));
        end
    end

    % phase gradient (vector) field
    Vx_flowmap_norm_phase = [];
    Vy_flowmap_norm_phase = [];
    vPhaseX = zeros(size(phaseSig)) ;
    vPhaseY = zeros(size(phaseSig)) ;
    cd(main_folder)
    for iTime = 1:size(phaseSig,3)
        for iX = 1:size(phaseSig,1)
            vPhaseX(iX,2:end-1,iTime) = (anglesubtract(phaseSig(iX,3:end,iTime),phaseSig(iX,1:end-2,iTime)))/2 ;
        end
        for iY = 1:size(phaseSig,2)
            vPhaseY(2:end-1,iY,iTime) = (anglesubtract(phaseSig(3:end,iY,iTime),phaseSig(1:end-2,iY,iTime)))/2 ;
        end
    end
    Vx_flowmap_norm_phase = vPhaseX;%./sqrt(vPhaseX.^2+vPhaseY.^2) ; % standardize vector length as 1
    Vy_flowmap_norm_phase = vPhaseY;%./sqrt(vPhaseX.^2+vPhaseY.^2) ; % standardize vector length as 1


    for time = 1:size(Vx_flowmap_norm_phase,3)
        temp1_vx = Vx_flowmap_norm_phase(:,:,time);
        temp1_vy = Vy_flowmap_norm_phase(:,:,time);
        [x,y] = meshgrid(1:size(temp1_vx,2),1:size(temp1_vx,1));
        [vorticity(:,:,time), ~]  = curl(x,y,temp1_vx,temp1_vy);  % curl value of phase vector field
    end





%% Identify Task Times

    foldername = [main_folder,'/Sample Data/Motor Task/Analysis'];
    cd(foldername)
    filename = ['spiral_detection_sub',num2str(subject),'.mat'];
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




    for t = 195
% Find graph at time frame
    fieldName = ['Frame' num2str(t)];
    G = detected_spirals.(fieldName);


    GeneratedGenerated_edge = [];

    ModulatedGenerated_edge = [];

    ModulatedModulated_edge = [];

    
    %%Vorticity Location Plottig code
    if plotting == 1
    figure(1)
    
    imagesc(sigBPass(:,:,t))
    hold on
    %contour(Q(:,:,t), [0 0],'color','g')
    contour(Q(:,:,t), [-0.05 -0.05],'color', 'k', 'linewidth', 1)
    %contour(vorticity(:,:,t), [-1 1],'color', 'k', 'linewidth', 1)

    quiver(Vx_flowmap_norm_phase(:,:,t),Vy_flowmap_norm_phase(:,:,t))
    h = plot(G,'XData', G.Nodes.PeakVorticityLoc(:,1), 'YData', G.Nodes.PeakVorticityLoc(:,2));
    h.NodeLabel = G.Nodes.NodeID;

    nodeColors = zeros(size(G.Nodes, 1), 3); % Default color white

    % Loop through each node
    for i = 1:size(G.Nodes, 1)
        if G.Nodes.Rotation(i) == "Positive"
            nodeColors(i,:) = [1,0,0]; % red for Positive
        elseif G.Nodes.Rotation(i) == "Negative"
            nodeColors(i,:) = [0,0,1]; % blue for Negative
        end
    end



    %h.NodeLabel = [];
    h.NodeLabel = [repmat("Gen Frame ", [numnodes(G),1])+num2str(G.Nodes.generationFrame)];
    
    h.NodeLabelColor = 'k';
    h.NodeColor = nodeColors;
    h.EdgeColor = 'k';
    hold off

    hold on
    for parcellation_ID = 1:22
        parcellation_template_1par = parcellation_template;
        parcellation_template_1par(isnan(parcellation_template_1par)) = 0;

        parcellation_template_1par(parcellation_template_1par~=parcellation_ID) = 0;

        B = bwboundaries(parcellation_template_1par(:,:,1),'noholes');
        for k = 1:length(B)
            boundary = B{k};
            plot(boundary(:,2),boundary(:,1),'-','linewidth',2,'color',[0,0,0])
        end
    end
    hold off
    set(gca,'ydir','normal')
    colorbar()
    colormap jet
    title(['Phase singularity and rotational contours over Spatially Filtered Amplitude, Frame',num2str(t)])
    %title(['Vortex peak Voritcity locations Frame' num2str(t)])
    % 
    xlim([0 250])
    ylim([0 180])
     % xlim([110 160])
     % ylim([100 150])
    pause(1)
    end

 

    %%
        %getcurrenttask()

        % Initialize currentTask to empty
        currentTask = [];

        if ismember(t, start_end_time_right_hand(1,1):start_end_time_right_hand(1,2))
            currentTask = [start_end_time_right_hand(1,1):start_end_time_right_hand(1,2)];
        elseif ismember(t, start_end_time_left_hand(1,1):start_end_time_left_hand(1,2))
            currentTask = [start_end_time_left_hand(1,1): start_end_time_left_hand(1,2)];
        elseif ismember(t, start_end_time_left_foot(1,1):start_end_time_left_foot(1,2))
            currentTask = [start_end_time_left_foot(1,1): start_end_time_left_foot(1,2)];
        elseif ismember(t, start_end_time_tongue(1,1):start_end_time_tongue(1,2))
            currentTask = [start_end_time_tongue(1,1): start_end_time_tongue(1,2)];
        elseif ismember(t, start_end_time_right_foot(1,1):start_end_time_right_foot(1,2))
            currentTask = [start_end_time_right_foot(1,1): start_end_time_right_foot(1,2)];

        elseif ismember(t, start_end_time_tongue(2,1):start_end_time_tongue(2,2))
            currentTask = [start_end_time_tongue(2,1): start_end_time_tongue(2,2)];
        elseif ismember(t, start_end_time_left_foot(2,1):start_end_time_left_foot(2,2))
            currentTask = [start_end_time_left_foot(2,1): start_end_time_left_foot(2,2)];
        elseif ismember(t, start_end_time_right_hand(2,1):start_end_time_right_hand(2,2))
            currentTask = [start_end_time_right_hand(2,1): start_end_time_right_hand(2,2)];
        elseif ismember(t, start_end_time_left_hand(2,1):start_end_time_left_hand(2,2))
            currentTask = [start_end_time_left_hand(2,1): start_end_time_left_hand(2,2)];
        elseif ismember(t, start_end_time_right_foot(2,1):start_end_time_right_foot(2,2))
            currentTask = [start_end_time_right_foot(2,1): start_end_time_right_foot(2,2)];
        end



        % Loop through edges to find the generation time of interactions
        for interacting_nodes = 1:numedges(G)

          %  if G.Edges.Weight(interacting_nodes) <= 2
                firstNodeId = G.Edges.EndNodes(interacting_nodes,1);
                secondNodeId = G.Edges.EndNodes(interacting_nodes,2);

           % if (G.Nodes.Rotation(firstNodeId) == "Positive" && G.Nodes.Rotation(secondNodeId) == "Negative") || (G.Nodes.Rotation(secondNodeId) == "Positive" && G.Nodes.Rotation(firstNodeId) == "Negative" )

                firstGenerationtime = G.Nodes.generationFrame(firstNodeId);
                secondGenerationtime = G.Nodes.generationFrame(secondNodeId);

                if ismember(firstGenerationtime,currentTask) && ismember(secondGenerationtime,currentTask)
                    % Generated / Generated Interaction
                    GeneratedGenerated_edge(end+1) = interacting_nodes;
                elseif  (~ismember(firstGenerationtime,currentTask) && ismember(secondGenerationtime,currentTask)) ||  (ismember(firstGenerationtime,currentTask) && ~ismember(secondGenerationtime,currentTask))
                    % Generated / Modulated Interaction
                    ModulatedGenerated_edge(end+1) = interacting_nodes;
                elseif ~ismember(firstGenerationtime,currentTask) && ~ismember(secondGenerationtime,currentTask)
                    % Modulated / Modulated Interaction
                    ModulatedModulated_edge(end+1) = interacting_nodes;
                end
           % end
            %end
        end



        Context_dependent_interactions = rmedge(G, [GeneratedGenerated_edge ModulatedModulated_edge]);

        nodeDegrees = degree(Context_dependent_interactions);
        nodesWithEdges = find(nodeDegrees > 0);
        % % Extract the subgraph that contains only nodes with edges
        Context_dependent_interactions = subgraph(Context_dependent_interactions, nodesWithEdges);



        task_generated_interactions = rmedge(G, [ModulatedGenerated_edge ModulatedModulated_edge]);
        nodeDegrees = degree(task_generated_interactions);
        nodesWithEdges = find(nodeDegrees > 0);
        % % Extract the subgraph that contains only nodes with edges
        task_generated_interactions = subgraph(task_generated_interactions, nodesWithEdges);

        if plotting == 1
            figure(2)

            g = plot(Context_dependent_interactions,'XData', Context_dependent_interactions.Nodes.PeakVorticityLoc(:,1), 'YData', Context_dependent_interactions.Nodes.PeakVorticityLoc(:,2));

            nodeColors = zeros(size(Context_dependent_interactions.Nodes, 1), 3); % Default color white

            % Loop through each node
            for i = 1:size(Context_dependent_interactions.Nodes, 1)
                if Context_dependent_interactions.Nodes.Rotation(i) == "Positive"
                    nodeColors(i,:) = [1,0,0]; % red for Positive
                elseif Context_dependent_interactions.Nodes.Rotation(i) == "Negative"
                    nodeColors(i,:) = [0,0,1]; % blue for Negative
                end
            end



            %g.NodeLabel = [];
            g.NodeLabel = [repmat("Gen Frame ", [numnodes(Context_dependent_interactions),1])+num2str(Context_dependent_interactions.Nodes.generationFrame)];
            g.NodeLabelColor = 'k';
            g.NodeColor = nodeColors;
            g.EdgeColor = 'k';

            hold on
            for parcellation_ID = 1:22
                parcellation_template_1par = parcellation_template;
                parcellation_template_1par(isnan(parcellation_template_1par)) = 0;

                parcellation_template_1par(parcellation_template_1par~=parcellation_ID) = 0;

                B = bwboundaries(parcellation_template_1par(:,:,1),'noholes');
                for k = 1:length(B)
                    boundary = B{k};
                    plot(boundary(:,2),boundary(:,1),'-','linewidth',2,'color',[0,0,0])
                end
            end
            hold off

            set(gca,'ydir','normal')
    title(['Context Dependent Interactions ' ,num2str(t) ])

    xlim([0 250])
    ylim([0 180])

    figure(3)

    p = plot(task_generated_interactions,'XData', task_generated_interactions.Nodes.PeakVorticityLoc(:,1), 'YData', task_generated_interactions.Nodes.PeakVorticityLoc(:,2));

    nodeColors = zeros(size(task_generated_interactions.Nodes, 1), 3); % Default color white

    % Loop through each node
    for i = 1:size(task_generated_interactions.Nodes, 1)
        if task_generated_interactions.Nodes.Rotation(i) == "Positive"
            nodeColors(i,:) = [1,0,0]; % red for Positive
        elseif task_generated_interactions.Nodes.Rotation(i) == "Negative"
            nodeColors(i,:) = [0,0,1]; % blue for Negative
        end
    end



    %p.NodeLabel = [];
    p.NodeLabel = [repmat("Gen Frame ", [numnodes(task_generated_interactions),1])+num2str(task_generated_interactions.Nodes.generationFrame)];
    p.NodeLabelColor = 'k';
    p.NodeColor = nodeColors;
    p.EdgeColor = 'k';

    hold on
    for parcellation_ID = 1:22
        parcellation_template_1par = parcellation_template;
        parcellation_template_1par(isnan(parcellation_template_1par)) = 0;

        parcellation_template_1par(parcellation_template_1par~=parcellation_ID) = 0;

        B = bwboundaries(parcellation_template_1par(:,:,1),'noholes');
        for k = 1:length(B)
            boundary = B{k};
            plot(boundary(:,2),boundary(:,1),'-','linewidth',2,'color',[0,0,0])
        end
    end
    hold off
    set(gca,'ydir','normal')
    title(['Task Generated Interactions ',num2str(t)])



    xlim([0 250])
    ylim([0 180])

end












        % Loop through edges to find the generation time of interactions
        for interacting_nodes = 1:numedges(G)

            %if G.Edges.Weight(interacting_nodes) <= 3
                firstNodeId = G.Edges.EndNodes(interacting_nodes,1);
                secondNodeId = G.Edges.EndNodes(interacting_nodes,2);

                if (G.Nodes.Rotation(firstNodeId) == "Positive" && G.Nodes.Rotation(secondNodeId) == "Negative") || (G.Nodes.Rotation(secondNodeId) == "Positive" && G.Nodes.Rotation(firstNodeId) == "Negative" )

                    firstGenerationtime = G.Nodes.generationFrame(firstNodeId);
                    secondGenerationtime = G.Nodes.generationFrame(secondNodeId);

                    if ismember(firstGenerationtime,currentTask) && ismember(secondGenerationtime,currentTask)
                        % Generated / Generated Interaction
                        GeneratedGenerated(t,subject) = GeneratedGenerated(t,subject) + 1;
                    elseif  (~ismember(firstGenerationtime,currentTask) && ismember(secondGenerationtime,currentTask)) ||  (ismember(firstGenerationtime,currentTask) && ~ismember(secondGenerationtime,currentTask))
                        % Generated / Modulated Interaction
                        GeneratedModulated(t,subject) = GeneratedModulated(t,subject) + 1;
                    elseif ~ismember(firstGenerationtime,currentTask) && ~ismember(secondGenerationtime,currentTask)
                        % Modulated / Modulated Interaction
                        ModulatedModulated(t,subject) = ModulatedModulated(t,subject) + 1;
                    end
                end
           % end
        end




%Number of Singularities and when they were evoked
    for nodeID  = 1:numnodes(G)
        generation_frame = G.Nodes.generationFrame(nodeID);


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
    end




   
    no_interacting_spirals(end+1) = size(Gsub.Edges.EndNodes,1);


    end % time loop end
end % subject loop end


%%
%plot(1:283,mode(interacting_spirals_all,1))

% PLOT INTERACTION TYPES OVER TIME 
figure(5)
%bar([GeneratedGenerated(1:283,1),GeneratedModulated(1:283,1),ModulatedModulated(1:283,1)],'stacked')
bar([mean(GeneratedGenerated(:,:),2),mean(GeneratedModulated(:,:),2),mean(ModulatedModulated(:,:),2)],'stacked')
title("Average Vortex Interactions At each time frame")
xlabel('time')
ylabel('Number of Interactions')
legend('GeneratedGenerated','GeneratedModulated','ModulatedModulated')
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

% PLOT NO. OF SINGULARITIES DIVIDED INTO THEIR GENRATION TIMES
figure(6)
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
    legend('RH1','LF1','T1','RF1','LH1','T2','LF2','RH2','LH2','RF2','Non-task evoked')
    title("Average Number of Vortices At each time frame")
    xlabel('time')
    ylabel('Number of Identified Singularites')
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
