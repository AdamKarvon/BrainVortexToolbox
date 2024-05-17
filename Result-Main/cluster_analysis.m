
restoredefaultpath
cd '\\suphys.physics.usyd.edu.au\akar5239\BrainVortexToolbox-main'
addpath(genpath([pwd]))
main_folder = pwd;

task = 'RF';

hemisphere = 1;

task_frame = 10;
movement_block = 2;
display_clusters = 0;
min_clusters = 2;
[block1_modulated_spiral_data, block1_generated_spiral_data,block1_modulated_clusters, block1_generated_clusters] = task_clustering_centres(movement_block, task, task_frame,min_clusters,display_clusters,hemisphere);

%%

% task_frame = 12;
% 
% block1_modulated_spiral_data = all_modulated_spiral_data.taskframe(task_frame);
% block1_generated_spiral_data = all_generated_spiral_data.taskframe(task_frame);
% block1_modulated_clusters  =  all_modulated_clusters.taskframe(task_frame);
% block1_generated_clusters = all_generated_clusters.taskframe(task_frame);

%% Spiral Duration within task specific clusters



for cluster = 1:3
    clear spiral_duration
    for pos_spiral = 1:size(block1_modulated_spiral_data.pos{cluster},1)
        spiral_duration(pos_spiral) = length(block1_modulated_spiral_data.pos{cluster}{pos_spiral,4});
    end 
    switch cluster
        case 1
            x1 = spiral_duration;
        case 2
            x2 = spiral_duration;
        case 3
            x3 = spiral_duration;
    end
end
for pos_spiral = 1:size(block1_modulated_spiral_data.posNoise,1)
    x4(pos_spiral) = length(block1_modulated_spiral_data.posNoise{pos_spiral,4});
end

g1 = repmat({'Cluster 1'},size(x1));
g2 = repmat({'Cluster 2'},size(x2));
g3 = repmat({'Cluster 3'},size(x3));
g4 = repmat({'Unclustered'},size(x4));

x = transpose([x1 x2 x3 x4]);
g = transpose([g1 g2 g3 g4]);

figure(3)
boxplot(x,g)

[xbar,s2,grp]  = grpstats(x,g,["mean","var","gname"]);

[h,p,ci,stats] = ttest2(x3,x4,0.05,'right');

%% Location of longest duration spirals

locs = block1_generated_spiral_data.negNoise(x4>50,3);

hold on
for i = 1:size(locs,1)
    for j = 1:size(locs{i},2)
        x = locs{i,1}{1,j}(1);
        y = locs{i,1}{1,j}(2);
        scatter(x,y,'filled','b')
    end
end

%% Clustered spiral dipoles
% stuck looking at static frames so for frame 12 well see if a subject that
% has a spiral in the largest pos cluster has any spirals in the other two
% neg clusters
clear Neighbouring_spirals
clear cluster_of_interest

spiral_type = 'generated';  % 'generated' or 'modulated'
focus_cluster = 'pos';  % look at pairs connected to the largest 'neg' or 'pos' cluster
 
cluster_order = 1;
clear spiral_pairs



switch spiral_type
      case 'modulated'
        if strcmp(focus_cluster, 'neg')
            clusters = block1_modulated_spiral_data.neg;
            selected_cluster = clusters{cluster_order};
            neighbouring_clusters = block1_modulated_spiral_data.pos;
            neighbouring_noise_pos = block1_modulated_spiral_data.posNoise;
            neighbouring_noise_neg = block1_modulated_spiral_data.negNoise;
            focus_color = 'b';
            neighbour_color = 'r';

        elseif strcmp(focus_cluster, 'pos')
            clusters = block1_modulated_spiral_data.pos;
            selected_cluster = clusters{cluster_order};
            neighbouring_clusters = block1_modulated_spiral_data.neg;
            neighbouring_noise_pos = block1_modulated_spiral_data.posNoise;
            neighbouring_noise_neg = block1_modulated_spiral_data.negNoise;
            focus_color = 'r';
            neighbour_color = 'b';

        end
    case 'generated'
        if strcmp(focus_cluster, 'neg')
            clusters = block1_generated_spiral_data.neg;
            selected_cluster = clusters{cluster_order}; 
            neighbouring_clusters = block1_generated_spiral_data.pos;
            neighbouring_noise_pos = block1_generated_spiral_data.posNoise;
            neighbouring_noise_neg = block1_generated_spiral_data.posNoise;
            focus_color = 'b';
            neighbour_color = 'r';

        elseif strcmp(focus_cluster, 'pos')
            clusters = block1_generated_spiral_data.pos;
            selected_cluster = clusters{cluster_order};
            neighbouring_clusters = block1_generated_spiral_data.neg;
            neighbouring_noise_pos = block1_generated_spiral_data.posNoise;
            neighbouring_noise_neg = block1_generated_spiral_data.negNoise;
            focus_color = 'r';
            neighbour_color = 'b';

        end
end


for focus_spiral =  1:size(selected_cluster,1)
    % Re-form data of selected cluster
    neighbours_of_interest = cell(1,1);

    spiral_time = selected_cluster{focus_spiral,4};
    spiral_loc = selected_cluster{focus_spiral,3};
    subject_id = selected_cluster{focus_spiral,1};
    spiral_no = 0;
    cluster_of_interest(subject_id,spiral_time) = spiral_loc;

    % Find surrounding spirals
    for icluster = 1:size(neighbouring_clusters,2) %along rows        
        spiral_time = neighbouring_clusters{icluster}(cell2mat(cellfun(@(c) c == subject_id , neighbouring_clusters{icluster}(:,1),'Uniform', false)),4);  
        spiral_loc = neighbouring_clusters{icluster}(cell2mat(cellfun(@(c) c == subject_id , neighbouring_clusters{icluster}(:,1),'Uniform', false)),3);
        
        for matching_subject_id = 1:size(spiral_time,1)
            spiral_no = spiral_no+1;
            neighbours_of_interest(spiral_no,spiral_time{matching_subject_id}) = spiral_loc{matching_subject_id};
        end      
    end
   
% NOISE/ Unclustered spirals
switch focus_cluster
    case 'pos'
        spiral_time = neighbouring_noise_neg(cell2mat(cellfun(@(c) c == subject_id , neighbouring_noise_neg(:,1),'Uniform', false)),4);
        spiral_loc = neighbouring_noise_neg(cell2mat(cellfun(@(c) c == subject_id , neighbouring_noise_neg(:,1),'Uniform', false)),3);
    
        for matching_subject_id = 1:size(spiral_time,1)
            spiral_no = spiral_no+1;
            neighbours_of_interest(spiral_no,spiral_time{matching_subject_id}) = spiral_loc{matching_subject_id};
        end

    case 'neg'
        spiral_time = neighbouring_noise_pos(cell2mat(cellfun(@(c) c == subject_id , neighbouring_noise_pos(:,1),'Uniform', false)),4);
        spiral_loc = neighbouring_noise_pos(cell2mat(cellfun(@(c) c == subject_id , neighbouring_noise_pos(:,1),'Uniform', false)),3);
    
        for matching_subject_id = 1:size(spiral_time,1)
            spiral_no = spiral_no+1;
            neighbours_of_interest(spiral_no,spiral_time{matching_subject_id}) = spiral_loc{matching_subject_id};
        end
end
    


    Neighbouring_spirals{subject_id} = neighbours_of_interest;  
end
%%


 for isubject = 1:150
     clf
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


 itime = 15:27;

 z = cluster_of_interest(isubject,itime);

 for i = 1:size(z,1)
     mask = ~cellfun(@isempty, z(i,:));
     nonEmptyTimes = z(i,mask);
     x = cellfun(@(c) c(1), nonEmptyTimes);
     y = cellfun(@(c) c(2), nonEmptyTimes);
     if ~isempty(x) && ~isempty(y)
         plot(x, y,'k-')
         end_point = plot(x(end), y(end),'s','MarkerFaceColor',focus_color);
     end
 end
 
if ~isempty(Neighbouring_spirals{isubject})
 z = Neighbouring_spirals{isubject}(:,itime);
 
     for i = 1:size(z,1)
         mask = ~cellfun(@isempty, z(i,:));
         nonEmptyTimes = z(i,mask);
         x = cellfun(@(c) c(1), nonEmptyTimes);   % The indexing here controls whether we look at all times for 1 subject or all subject for 1 time
         y = cellfun(@(c) c(2), nonEmptyTimes);
         if ~isempty(x) && ~isempty(y)
             plot(x, y,'k-')
             end_point = plot(x(end), y(end),'o','MarkerFaceColor',neighbour_color);
         end
     end
end
pause()
 end %subject for loop end

%%

for pairs = 1: size(spiral_pairs,1)

    EmptyCells = cellfun('isempty', spiral_pairs(pairs,:));


    if sum(EmptyCells) <= no_of_neighbouring_clusters %spiral missing in most prominent location
        focus_path =  cellfun(@(C) C(1,3), spiral_pairs(pairs,1));
        non_empty_clusters = spiral_pairs(pairs,~cellfun('isempty', spiral_pairs(pairs,:)));
        neighbouring_path = cellfun(@(C) C(1,3), non_empty_clusters(1,2:end));


        figure(2)
        hold on
        title('Subject has spiral in all clusters')
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

        x_coords = cellfun(@(c) c(1), focus_path{1});
        y_coords = cellfun(@(c) c(2), focus_path{1});
        plot(x_coords, y_coords, 'k-');
        end_point = plot(x_coords(end), y_coords(end),'o','MarkerFaceColor',focus_color);

        for path = 1:size(neighbouring_path,2)
            x_coords = cellfun(@(c) c(1), neighbouring_path{path});
            y_coords = cellfun(@(c) c(2), neighbouring_path{path});

            plot(x_coords, y_coords, 'k--');
            plot(x_coords(end), y_coords(end),'o','MarkerFaceColor',neighbour_color);
        end
    end
    if EmptyCells(3) == 1 && sum(EmptyCells) <= 2 % spiral missing in 2nd most prominent location
        focus_path =  cellfun(@(C) C(1,3), spiral_pairs(pairs,1));
        non_empty_clusters = spiral_pairs(pairs,~cellfun('isempty', spiral_pairs(pairs,:)));
        neighbouring_path = cellfun(@(C) C(1,3), non_empty_clusters(1,2:end));

        figure(2)
        hold on
        title('Subject has spiral pair between 3 clusters')
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

        x_coords = cellfun(@(c) c(1), focus_path{1});
        y_coords = cellfun(@(c) c(2), focus_path{1});
        plot(x_coords, y_coords, 'k-');
        end_point = plot(x_coords(end), y_coords(end),'o','MarkerFaceColor',focus_color);
        for path = 1:size(neighbouring_path,2)
            x_coords = cellfun(@(c) c(1), neighbouring_path{path});
            y_coords = cellfun(@(c) c(2), neighbouring_path{path});

            plot(x_coords, y_coords, 'k--');
            plot(x_coords(end), y_coords(end),'o','MarkerFaceColor',neighbour_color);
        end
    end
    if EmptyCells(4) == 1 && sum(EmptyCells) <= 2 % spiral missing in 3rd most prominent location
        focus_path =  cellfun(@(C) C(1,3), spiral_pairs(pairs,1));
        non_empty_clusters = spiral_pairs(pairs,~cellfun('isempty', spiral_pairs(pairs,:)));
        neighbouring_path = cellfun(@(C) C(1,3), non_empty_clusters(1,2:end));


        figure(3)
        hold on
        title('Subject has spiral pair between 2 clusters')
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

        x_coords = cellfun(@(c) c(1), focus_path{1});
        y_coords = cellfun(@(c) c(2), focus_path{1});
        plot(x_coords, y_coords, 'k-');
        end_point = plot(x_coords(end), y_coords(end),'o','MarkerFaceColor',focus_color);
        for path = 1:size(neighbouring_path,2)
            x_coords = cellfun(@(c) c(1), neighbouring_path{path});
            y_coords = cellfun(@(c) c(2), neighbouring_path{path});

            plot(x_coords, y_coords, 'k--');
            plot(x_coords(end), y_coords(end),'o','MarkerFaceColor',neighbour_color);
        end
    end
    if EmptyCells(5) == 1  % spiral missing in 4th most prominent location
        focus_path =  cellfun(@(C) C(1,3), spiral_pairs(pairs,1));
        non_empty_clusters = spiral_pairs(pairs,~cellfun('isempty', spiral_pairs(pairs,:)));
        neighbouring_path = cellfun(@(C) C(1,3), non_empty_clusters(1,2:end));

        figure(4)
        hold on
        title('Subject has spiral pair in 1 clusters')
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

        x_coords = cellfun(@(c) c(1), focus_path{1});
        y_coords = cellfun(@(c) c(2), focus_path{1});
        plot(x_coords, y_coords, 'k-');
        end_point = plot(x_coords(end), y_coords(end),'o','MarkerFaceColor',focus_color);
        for path = 1:size(neighbouring_path,2)
            x_coords = cellfun(@(c) c(1), neighbouring_path{path});
            y_coords = cellfun(@(c) c(2), neighbouring_path{path});

            plot(x_coords, y_coords, 'k--');
            plot(x_coords(end), y_coords(end),'o','MarkerFaceColor',neighbour_color);
        end
    end
    if sum(EmptyCells)>4  % subject has no obvious task specific spiral pairs
        continue
    end
end
%% Spiral Duration as function of frame specific Spiral size
spiralSize = block1_generated_clusters.neg{2}(:,5);



%meanSpiralSize = cellfun(@size , RH1_generated_spiral_data.neg{2}{:,5});



figure(3)
hold on
scatter(spiralSize, x2);
f = fit(spiralSize, transpose(x2),'b*x^m'); 

plot(f,'k')


%% Spiral Duration as a function of distance from cluster centre

X = block1_modulated_clusters.neg{2}(:,3);
Y = block1_modulated_clusters.neg{2}(:,4);
cluster_centre = mean([X,Y]);
distances = sqrt((X - cluster_centre(1)).^2 + (Y - cluster_centre(2)).^2);



figure(3)
hold on
scatter(distances,x2)
f = fit(distances, transpose(x2),'b*x^m'); 

plot(f,'k')

xlabel('Distance to cluster centre')
ylabel('Spiral DUration')

%% Spiral displacement from initial location
for cluster = 1:3
distance = zeros(size(block1_modulated_clusters.neg{cluster},1),1);
for pos_spiral = 1:size(block1_modulated_clusters.neg{cluster},1)
    spiral_locs = block1_modulated_spiral_data.neg{cluster}{pos_spiral,3};
    
    X = block1_modulated_clusters.neg{cluster}(pos_spiral,3);
    Y = block1_modulated_clusters.neg{cluster}(pos_spiral,4);
    clustered_loc = [X Y];
    initial_loc = cell2mat(spiral_locs(1));
    distance(pos_spiral) = pdist([clustered_loc; initial_loc],'euclidean');
end

switch cluster
    case 1
        d1 = distance;
    case 2
        d2 = distance;
    case 3
        d3 = distance;
end

end
g1 = repmat({'Cluster 1'},size(d1));
g2 = repmat({'Cluster 2'},size(d2));
g3 = repmat({'Cluster 3'},size(d3));


D = [d1; d2; d3];
G = [g1; g2; g3];

figure(4)
boxplot(D,G)


%% Cluster Trajectories prior to cluster
cluster = 1;

for spiral = 1:size(block1_modulated_clusters.neg{cluster})
    clustered_loc = [block1_modulated_clusters.neg{cluster}(spiral,3) block1_modulated_clusters.neg{cluster}(spiral,4)];


    spiral_path = block1_modulated_spiral_data.neg{cluster}{spiral,3};
    figure(2)
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
    
    % Use cellfun to extract x and y coordinates
    x_coords = cellfun(@(c) c(1), spiral_path);
    y_coords = cellfun(@(c) c(2), spiral_path);

    % Plot the line connecting all the points
    plot(x_coords, y_coords, 'k-');
    plot(x_coords(1), y_coords(1),'o','MarkerFaceColor','b');
    xlabel('X-coordinate');
    ylabel('Y-coordinate');
    title('Line Connecting Points');

end


% relating spiral pairs
% 
% for spiral
% 
%   if ismember(subject,) && ismember(subject)



%% Spiral Size as a function of distance from spiral centre
X = block2_generated_clusters.pos{2}(:,3);
Y = block2_generated_clusters.pos{2}(:,4);
cluster_centre = mean([X,Y]);

distances = sqrt((X - cluster_centre(1)).^2 + (Y - cluster_centre(2)).^2);
figure(3)
hold on
scatter(distances,block2_generated_clusters.pos{2}(:,5))
f = fit(distances,block2_generated_clusters.pos{2}(:,5),'b*x+c'); 

plot(f,'k')

xlabel('Distance to cluster centre')
ylabel('Spiral Size')

% 
% figure(1)
% hold on
%  histogram(RH1_modulated_clusters.pos{1}(:,5),10)
%  xline(mean(RH1_modulated_clusters.pos{1}(:,5)))
% figure(2)
%  hold on
%  histogram([RH1_generated_clusters.pos{1}(:,5)],10)
%  xline(mean(RH1_generated_clusters.pos{1}(:,5)))

 % xline(mean(RH1_generated_clusters.pos{2}(:,5))+sqrt(var(RH1_generated_clusters.pos{2}(:,5))))
 % xline(mean(RH1_generated_clusters.pos{2}(:,5))-sqrt(var(RH1_generated_clusters.pos{2}(:,5))))






%% Number of modulated and generated spirals as a fucntion of time


%% Plotting Code
for num = task_frame
% block1_modulated_clusters = all_modulated_clusters.taskframe(num);
% block1_generated_clusters = all_generated_clusters.taskframe(num);
plot_noise = 1;
modulated_clusters = block1_modulated_clusters;
generated_clusters = block1_generated_clusters;

[~,id] = sort( cellfun(@(c) size(c,1), modulated_clusters.pos ),'descend');
modulated_clusters.pos = modulated_clusters.pos(id);
if hemisphere == 1
    load('parcellation_template.mat')
elseif hemisphere ==2
    load('parcellation_template22_RightBrain_subject1-100.mat')
    parcellation_template = parcellation_template22_RightBrain_100sub(:,:,1);
end



disp (['Pos Modulated clusters' ])
disp (modulated_clusters.pos)
disp (['Neg  Modulated clusters' ])
disp(modulated_clusters.neg)


disp (['Pos Generated clusters' ])
disp (generated_clusters.pos)
disp (['Neg Generated clusters' ])
disp(generated_clusters.neg)

figure(1)
clf
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
x = [];
y = [];
i = [];
j = [];
G = [];
H = [];


for cluster_id = 1:size(modulated_clusters.neg,1)
    x = [x; modulated_clusters.neg{cluster_id}(:,3)];
    y = [y; modulated_clusters.neg{cluster_id}(:,4)];
    G = [G; size(modulated_clusters.neg{cluster_id},1)*ones(size(modulated_clusters.neg{cluster_id}(:,3),1),1)];
end
for cluster_id = 1:size(modulated_clusters.pos,1)
    i = [i; modulated_clusters.pos{cluster_id}(:,3)];
    j = [j; modulated_clusters.pos{cluster_id}(:,4)];
    H = [H; size(modulated_clusters.pos{cluster_id},1)*ones(size(modulated_clusters.pos{cluster_id}(:,3),1),1)];
end
if plot_noise == 1
    s3 = scatter(modulated_clusters.posNoise(:,3),modulated_clusters.posNoise(:,4),'rx','MarkerEdgeAlpha',.2);
    s4 = scatter(modulated_clusters.negNoise(:,3),modulated_clusters.negNoise(:,4),'bx','MarkerEdgeAlpha',.2);
end
xticklabels([""]);
yticklabels([""]);
gs1 = gscatter(x,y,G,'b','+*sd^v><ph',8);
gs2 = gscatter(i,j,H,'r','+*sd^v><ph',8);
axis equal
a1=gca;
% title(a1,['Modulated spiral clusters ', num2str(num) ,' frames After onset'])
title(a1,['Task Modulated'],'FontSize',20)
a2=axes('position',get(gca,'position'),'visible','off');
leg1 = legend(a1, gs1);
fontsize(leg1,20,'points')
leg2 = legend(a2, gs2);
title(leg1,'CW Cluster Sizes');
title(leg2,'ACW Cluster Sizes'); 
fontsize(leg2,20,'points')
axis equal


figure(2)
clf
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
x = [];
y = [];
i = [];
j = [];
G = [];
H = [];
for cluster_id = 1:size(generated_clusters.neg,1)
    x = [x; generated_clusters.neg{cluster_id}(:,3)];
    y = [y; generated_clusters.neg{cluster_id}(:,4)];
    G = [G; size(generated_clusters.neg{cluster_id},1)*ones(size(generated_clusters.neg{cluster_id}(:,3),1),1)];
end
for cluster_id = 1:size(generated_clusters.pos,1)
    i = [i; generated_clusters.pos{cluster_id}(:,3)];
    j = [j; generated_clusters.pos{cluster_id}(:,4)];
    H = [H; size(generated_clusters.pos{cluster_id},1)*ones(size(generated_clusters.pos{cluster_id}(:,3),1),1)];
end
xticklabels([""]);
yticklabels([""]);
if plot_noise == 1
    s3 = scatter(generated_clusters.posNoise(:,3),generated_clusters.posNoise(:,4),'rx','MarkerEdgeAlpha',.2,'DisplayName','unclustered');
    s4 = scatter(generated_clusters.negNoise(:,3),generated_clusters.negNoise(:,4),'bx','MarkerEdgeAlpha',.2,'DisplayName','unclustered');
end
gs1 = gscatter(x,y,G,'b','+*sd^v><ph',8);
gs2 = gscatter(i,j,H,'r','+*sd^v><ph',8);

axis equal
a1=gca;
a2=axes('position',get(gca,'position'),'visible','off');
% title(a1,['Generated spiral clusters ', num2str(num) ,' frames After onset'])
title(a1,['Task Generated'],'FontSize',20)
leg1 = legend(a1, gs1);
fontsize(leg1,20,'points')

leg2 = legend(a2, gs2);
title(leg1,'CW Cluster Sizes');
title(leg2,'ACW Cluster Sizes'); 
fontsize(leg2,20,'points')

hold off
set(gca,'ydir','normal')
axis equal

end