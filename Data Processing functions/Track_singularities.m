function spiral_detection = sprial_detect(subject)
%%
for hemisphere = 2
    disp('done')
for subject = 1:150
restoredefaultpath
cd '\\suphys.physics.usyd.edu.au\akar5239\BrainVortexToolbox-main'

addpath(genpath([pwd]))
main_folder = pwd;


 % 1 for left hemisphere, 2 for right hemisphere, keep this as 1 for demonstration purpose

% flagTask:
% 1 = language task, original 100 subjects;
% 2 = language task, additional 100 subjects;
% 3 = working memory task;
% 4 = Motor Task
flagTask = 4;


cd(main_folder)


% load spatiotemporal bandpass filtered fMRI signal file
foldername = [main_folder,'/Sample Data/Motor Task/Preprocessed Data'];
cd(foldername)
if hemisphere == 1            % spatio
    filename = ['Preprocessed_spatiotemporalbandpass_data_motor_task_LEFT_sub',num2str(subject),'.mat'];
elseif hemisphere == 2
    filename = ['Preprocessed_spatiotemporalbandpass_data_motor_task_RIGHT_sub',num2str(subject),'.mat'];
end
load(filename)
disp(subject)

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


[detected_spirals,interaction_table] = singularity_tracking(vPhaseX,vPhaseY,phaseSig);  


% Save File
    folder_name = [main_folder,'/Sample Data/Motor Task/Analysis'];
    cd(folder_name) 
if hemisphere == 1
    filename = ['spiral_detection_sub_Left', num2str(subject), '.mat'];
elseif hemisphere == 2 
    filename = ['spiral_detection_sub_Right', num2str(subject), '.mat'];
end
    save([folder_name,'/',filename],'detected_spirals','interaction_table')



end
end

