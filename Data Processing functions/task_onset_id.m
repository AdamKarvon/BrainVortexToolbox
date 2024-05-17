function [spontaneuos_vortices, modulated_vortices] = task_onset_id(spiral_filt, centres, sizes , blocks)

%%%%%%%%%%%%%%%%%
%Identifies Vortices that exist prior to the task onset, and Vortices that appear after the task onset

% records:
% -Centre location
% -Vortex lifetime (modulation/generation time)
% -Vortex size (as estimated is spiral_detection_surfilt.m)

%NOTE: CW and A-CW vortices are considered seperately for the sake of clarity
%%%%%%%%%%%%%%%%%
dt= 0.72;
% Pre-allocating variables
vortexIndex_block1 = 0;
vortexIndex_block2 = 0;
modulatedIndex_block1 = 0;
modulatedIndex_block2 = 0;
possible_modulation = 0;


% Generated Vortices in first block of motion 
for ipatt = 1:size(spiral_filt,1)
    for time = blocks(1,:)
% If a vortex appears in the detected vortices when there wasn't a
% vortex in the previous time step check if it has been generated
% during the movement block
        if isempty(spiral_filt{ipatt,time-1}) && ~isempty(spiral_filt{ipatt,time})
            possible_generation = ipatt;
            generation_time = find(~cellfun('isempty', spiral_filt(possible_generation,:)));

            if size(generation_time,2) > 3 && ismember(generation_time(1),blocks(1,:))
                %vortexIndex keeps track of the number of generated vortices
                vortexIndex_block1 = vortexIndex_block1 + 1;
                block1_generated_vortices{vortexIndex_block1,1} = possible_generation;
                block1_generated_vortices{vortexIndex_block1,2} = centres(ipatt,generation_time);
                block1_generated_vortices{vortexIndex_block1,3} = generation_time;
                block1_generated_vortices{vortexIndex_block1,4} = sizes(ipatt,generation_time);
                 
                dx = diff(cellfun(@(c)c(:,1) ,centres(ipatt,generation_time)));
                dy = diff(cellfun(@(c)c(:,2) ,centres(ipatt,generation_time)));
                vx = dx/dt;
                vy = dy/dt;
                velocities = [vx; vy];
                block1_generated_vortices{vortexIndex_block1,5} = velocities;
            end

        elseif ~isempty(spiral_filt{ipatt,blocks(1,1)-1}) && ~isempty(spiral_filt{ipatt,blocks(1,1)}) && ipatt ~= possible_modulation
            possible_modulation = ipatt;
            modulation_time = find(~cellfun('isempty', spiral_filt(possible_modulation,:)));
% Now we check what vortices existed before the task block and continued to
% exist for some time during the task period
            if size(modulation_time,2) > 3 && modulation_time(1) < blocks(1,1)      %Considering strictly less than task onset 
                modulatedIndex_block1 = modulatedIndex_block1 + 1;
                block1_modulated_vortices{modulatedIndex_block1,1} = possible_modulation;
                block1_modulated_vortices{modulatedIndex_block1,2} = centres(ipatt,modulation_time);

                block1_modulated_vortices{modulatedIndex_block1,3} = modulation_time;
                block1_modulated_vortices{modulatedIndex_block1,4} = sizes(ipatt,modulation_time);
                
                dx = diff(cellfun(@(c)c(:,1) ,centres(ipatt,modulation_time)));
                dy = diff(cellfun(@(c)c(:,2) ,centres(ipatt,modulation_time)));
                vx = dx/dt;
                vy = dy/dt;
                velocities = [vx; vy];
                block1_modulated_vortices{modulatedIndex_block1,5} = velocities;
            end
        end
    end

    for time = blocks(2,:)
        if isempty(spiral_filt{ipatt,time-1}) && ~isempty(spiral_filt{ipatt,time})
            possible_generation = ipatt;
            generation_time = find(~cellfun('isempty', spiral_filt(possible_generation,:)));
            if size(generation_time,2) > 3 && ismember(generation_time(1),blocks(2,:))
                vortexIndex_block2 = vortexIndex_block2 + 1;
                block2_generated_vortices{vortexIndex_block2,1} = possible_generation;
                block2_generated_vortices{vortexIndex_block2,2} = centres(ipatt,generation_time);

                block2_generated_vortices{vortexIndex_block2,3} = generation_time;
                block2_generated_vortices{vortexIndex_block2,4} = sizes(ipatt,generation_time);
                
                
                dx = diff(cellfun(@(c)c(:,1) ,centres(ipatt,generation_time)));
                dy = diff(cellfun(@(c)c(:,2) ,centres(ipatt,generation_time)));
                vx = dx/dt;
                vy = dy/dt;
                velocities = [vx; vy];
                block2_generated_vortices{vortexIndex_block2,5} = velocities;
            end
        elseif ~isempty(spiral_filt{ipatt,blocks(2,1)-1}) && ~isempty(spiral_filt{ipatt,blocks(2,1)}) && ipatt ~= possible_modulation
            possible_modulation = ipatt;
            modulation_time = find(~cellfun('isempty', spiral_filt(possible_modulation,:)));
            if size(modulation_time,2) > 3 && modulation_time(1) < blocks(2,1)   %Considering strictly less than task onset 
                modulatedIndex_block2 = modulatedIndex_block2 +1;
                block2_modulated_vortices{modulatedIndex_block2,1} = possible_modulation;
                block2_modulated_vortices{modulatedIndex_block2,2} = centres(ipatt,modulation_time);
                block2_modulated_vortices{modulatedIndex_block2,3} = modulation_time;
                block2_modulated_vortices{modulatedIndex_block2,4} = sizes(ipatt,modulation_time);
                
                dx = diff(cellfun(@(c)c(:,1) ,centres(ipatt,modulation_time)));
                dy = diff(cellfun(@(c)c(:,2) ,centres(ipatt,modulation_time)));
                vx = dx/dt;
                vy = dy/dt;
                velocities = [vx; vy];
                block2_modulated_vortices{modulatedIndex_block2,5} = velocities;

            end
        end
    end
end

spontaneuos_vortices{1} = block1_generated_vortices;
spontaneuos_vortices{2} = block2_generated_vortices;
modulated_vortices{1} = block1_modulated_vortices ;
modulated_vortices{2} = block2_modulated_vortices;

end