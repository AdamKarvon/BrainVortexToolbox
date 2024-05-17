# BrainVortexToolbox
MATLAB toolbox to automatically detect and analyse spiral wave patterns in the fMRI data, developed by Dr. Pulin Gong's group at University of Sydney. 

## Instructions for use
System dependencies: N/A <br />
Software version: MATLAB 2016b and above (has been tested on MATLAB 2016b and 2022b) <br />
Hardware requirement: N/A

Data format: fMRI data in standard CIFTI grayordinate space, comprising 32K cortical vertices <br />
Data tested: 150 subjects from the Human Connectome Project (HCP)<br /> [https://db.humanconnectome.org/data/projects/HCP_1200]

### Launch: <br />
Sample data is not avaialble in GitHub due to size limitation, please download from the HCP site link above. No subject ID is provided as all subjects should be selected randomly from a corhort of 1200 subjects.

Please download all folders from this repository and allocate the raw fMRI data files (i.e., 'tfMRI_LANGUAGE_LR_Atlas.dtseries.nii') and structural data files (i.e.,'L.flat.32k_fs_LR.surf.gii') downloaded from HCP database into subfolders named by the subject ID under 'Raw Data' and 'Data Pos' folders, repectively. For task state data, please also download the task label files for each subject (witihn the 'EVs' subfolder next to the raw fMRI data file).

For example, the data under subject ID 100206 recorded during a language task should be allocated in the following folders: <br />
Raw fMRI data: '/main_folder/Sample Data/Language Task Original 100 sub/Raw Data/100206/tfMRI_LANGUAGE_LR_Atlas.dtseries.nii'; <br />
Structural data: '/main_folder/Sample Data/Language Task Original 100 sub/Data Pos/100206/L.flat.32k_fs_LR.surf.gii'; <br />
Task label files:  '/main_folder/Sample Data/Language Task Original 100 sub/Raw Data/100206/EVs/present_math.txt'; <br />




## Authors
Brain Vortex Toolbox functions -
* **Yiben Xu** - yixu4976@uni.sydney.edu.au
* **Xian Long** - [Xian Long](https://github.com/longxian319)
* **Pulin Gong** - pulin.gong@sydney.edu.au
Results, Extra functions, Data Processing -
* **Adam Karvon** - akar5239@uni.sydney.edu.au







