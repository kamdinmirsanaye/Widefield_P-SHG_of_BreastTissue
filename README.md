# Widefield P-SHG of BreastTissue
The enclosed software was created to carryout analysis of widefield polarimetric second-harmonic generation images of human breast tissue.
## System Requirements
The standalone software was created for Windows 10 64-bit, however it has also been tested successfully on a Mac platform.  
The program is written in MATLAB 2021a. Version 2017a is required for displaying figures, however, the remainder of code should run on older versions of MATLAB as well.  
## Installation
#### MATLAB Access
Various versions of MATLAB are available on [MATLAB's website](https://www.mathworks.com/products/matlab.html). Software free trials, as well as full licenses are available on the same website. Alternatively, the program can be executed on [MATLAB Online](https://www.mathworks.com/products/matlab-online.html) with an appropriate license.
#### Software Access
To successfully execute the program, download and unzip all files present in this [repository](https://github.com/kamdinmirsanaye/Widefield_P-SHG_of_BreastTissue) on a local drive.
## Demo
#### File/Folder Structure
- The main executible program is "WidefieldPolarimetricSHGofBreastTissue.m", which requires "DivideImage.m", "GLCMFeat.m", and "LogisticRegCode.m", as well as "HeatmapMaker.m" when display graphics option is selected.  
- In addition to the required MATLAB code files, 5 folders titled: "Background", "Data", "DataSubset", "ImageResults", and "ExpectedResults" can be found in the downloaded repository.  
- The "Background" folder contains 2 sets of background images, each corresponding to approximately half of the data. These images will be subtracted from the data in the program.  
- The "Data" folder includes 4 subfolders corresponding to normal and 3 tumor groups. Each subfolder contains a complete set of analyzed images of each tissue type.  
- In order to test the code on a smaller subset of the data, please use the "DataSubset" folder.  
- To avoid generating numerous images of the data, pregenerated polarimetric parameter images are provided in the "ImageResults" folder. 
- The expected results of the complete analysis are provided in folder "ExpectedResults".
- Ensure all above files and folders are present under a unified folder 
#### Running the program
- Open "WidefieldPolarimetricSHGofBreastTissue.m" in MATLAB
- Select appropriate options under section entitled: "Tunable Options". It is important to note that current options are set to replicate the results of the submitted manuscript.
- If you wish to display polarimetric parameter images, set "DisplayImages" option to "1". If you wish to display the corresponding histograms, set "DisplayHistograms" option to "1".
- If using a Mac computer, direction of the slashes will need to be reversed. To do so:
  - Press Command+F
  - Type " '\' " for "Find what"
  - Type " '/' " for "Replace with"
  - Click on "Replace all"
- Run the program
- The first window prompt will ask you to select the background images. Please select the "Background" folder.
- The second window prompt will ask you to select the data for analysis. Please select either the full dataset ("Data" folder) or the subset of the data for testing the code ("SubsetData"). Analysis of the full dataset (withtout displaying the graphics) takes less than 12 minutes, while the analysis of the subset take less than 2 minutes.
- At this stage, the program will read the background and data, compute the polarimetric and texture parameters, carryout classification of normal and tumor tissue images, perform statistical multiple comparisons testing, and display the boxplots of all parameters. 
- In addition, 4 results tables will be displayed in the command prompt, showing the repeated cross-validation mean, repeated cross-validation standard deviation, prediction results, as well as, multiple comparisons results. Screnshots of the expected results can be found in the "ExpectedResults" folder.
