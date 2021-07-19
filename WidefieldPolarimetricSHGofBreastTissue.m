%%%%%%%%%%%%%%%% Widefield Polarimetry Analysis Program %%%%%%%%%%%%%%%%%%%

% Name:     WidefieldPolarimetricSHGofBreastTissue.m
% Authors:  Kamdin Mirsanaye & Yasmeen Kamaliddin
% Date:     July, 2021

%   This code computes polarimetric paramters of images gathered from the
%   widefield nonlinear microscope in Barzda lab at University of Toronto.
%   It then automatically calculates statistics and five textural features 
%   (contrast, correlation, entropy, ASM, and IDM) of an entire data set of
%   P-SHG images.

% NECESSARY FUNCTIONS: HeatmapMaker.m, DivideImage.m, GLCMFeat.m, LogisticRegCode.m

% Copyright (C) Kamdin Mirsanaye & Yasmeen Kamaliddin 2021-

%% %%%%%%%%%%%%%%%%%%%%%%%%%%% Tunable Options %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Signal to Noise Ratio Cutoff
SNRcutoff = 3;

% Minimum and Maximum Paramter Ranges
minR=1.4; maxR=2.4;
minLD=-1.5; maxLD=1.5;
minCD=-0.8; maxCD=0.8;
minDCP=0.5; maxDCP=1;
MinSHGInt=20; MaxSHGInt=400;

% Texture Analysis parameters
d = 1;
numOfSubImagesVec  = [64]; % Include any of the following subdivision levels in the array: 
                           % 1,4,16,64,256,1024,4096,16384,65536,262144,104856
                           % Please note that subdivision above 256 are unnecessary and result in
                           % significantly longer computation times

% Display Polarimetric Parameter Images
DisplayImages = 0; % set to 1 to enable

% Display Polarimetric Parameter Histograms
DisplayHistograms = 0; % set to 1 to enable

% Perform Multiple Comparisons
MultipleComparisons=1; % set to 1 to enable

% Perform Texture Analysis
TextureAnalysis=1; % set to 1 to enable

% Perform Classification
Classification=1; % set to 1 to enable

% Display Parameter Boxplots
Boxplots=1; % set to 1 to enable

%% %%%%%%%%%%%%%%%%% Turn off unnecessary warnings %%%%%%%%%%%%%%%%%%%%%%%%

addpath(fileparts(mfilename('fullpath')))
cd(fileparts(mfilename('fullpath')))
warning('off','all')
warning

%% %%%%%%%%%%%%%%%%%%%%% Directory Specifications %%%%%%%%%%%%%%%%%%%%%%%%%

% Choose Background Directory
BGNDPath = uigetdir(pwd, 'Select Background Folder');
BGNDPath1 = dir(BGNDPath);
BGNDPath2 = BGNDPath1([BGNDPath1(:).isdir]); 
BGNDPath2 = BGNDPath2(~ismember({BGNDPath2(:).name},{'.','..'}));
BGNDPaths = cell(1,length(BGNDPath2));
BGNDNumber = zeros(1,length(BGNDPath2));
for i=1:length(BGNDPath2)
    BGNDName = BGNDPath2(i).name;
    BGNDNumber(i) = str2num(BGNDName(1));
    BGNDPaths{i}=[BGNDPath '\' BGNDName];
end

% Choose Group Directories
GroupPath = uigetdir(pwd, 'Select Groups Folder');
GroupPath1 = dir(GroupPath);
GroupPath2 = GroupPath1([GroupPath1(:).isdir]); 
GroupPath2 = GroupPath2(~ismember({GroupPath2(:).name},{'.','..'}));
GroupPaths = cell(1,length(GroupPath2));
GroupNumber = zeros(1,length(GroupPath2));
for i=1:length(GroupPath2)
    GroupName = GroupPath2(i).name;
    GroupNumber(i) = str2num(GroupName(1));
    GroupPaths{i}=[GroupPath '\' GroupName];
end

%% %%%%%%%%%%%% Load Polarization-Dependent Background Data %%%%%%%%%%%%%%%

fprintf('Loading data ...\n');

bgndMatC1C1=[];bgndMatC2C1=[];bgndMatL1C1=[];bgndMatL2C1=[];
bgndMatC1C2=[];bgndMatC2C2=[];bgndMatL1C2=[];bgndMatL2C2=[];
bgndMatC1L2=[];bgndMatC2L2=[];bgndMatL1L2=[];bgndMatL2L2=[];
bgndMatC1L1=[];bgndMatC2L1=[];bgndMatL1L1=[];bgndMatL2L1=[];

for bg = 1:length(BGNDPaths)
    cd(BGNDPaths{bg})    
    % Variable Initialization (Notation: PSG-PSA)

    bgndC1C1=[];bgndC2C1=[];bgndL1C1=[];bgndL2C1=[];
    bgndC1C2=[];bgndC2C2=[];bgndL1C2=[];bgndL2C2=[];
    bgndC1L2=[];bgndC2L2=[];bgndL1L2=[];bgndL2L2=[];
    bgndC1L1=[];bgndC2L1=[];bgndL1L1=[];bgndL2L1=[];

    F = dir('*.tif');

    for iii = 1:length(F)
        Fname = F(iii).name;
        ind1 = strfind(Fname,'_');
        ind2 = strfind(Fname,'.');
        ImageNumber = str2double(Fname(ind1(end)+1:ind2-1));
        if (ImageNumber<=length(F)/2 && mod(ImageNumber,8)==1)
            FileData = imread(Fname, 'tif');
            FileData = reshape(FileData,[],1);
            FileData = double(FileData);
            bgndC1C1 = [bgndC1C1 FileData];
        elseif (ImageNumber<=length(F)/2 && mod(ImageNumber,8)==2)
            FileData = imread(Fname, 'tif');
            FileData = reshape(FileData,[],1);
            FileData = double(FileData);
            bgndC1C2 = [bgndC1C2 FileData];
        elseif (ImageNumber<=length(F)/2 && mod(ImageNumber,8)==3)
            FileData = imread(Fname, 'tif');
            FileData = reshape(FileData,[],1);
            FileData = double(FileData);
            bgndC1L2 = [bgndC1L2 FileData];
        elseif (ImageNumber<=length(F)/2 && mod(ImageNumber,8)==4)
            FileData = imread(Fname, 'tif');
            FileData = reshape(FileData,[],1);
            FileData = double(FileData);
            bgndC1L1 = [bgndC1L1 FileData];
        elseif (ImageNumber<=length(F)/2 && mod(ImageNumber,8)==5)
            FileData = imread(Fname, 'tif');
            FileData = reshape(FileData,[],1);
            FileData = double(FileData);
            bgndC2C1 = [bgndC2C1 FileData];
        elseif (ImageNumber<=length(F)/2 && mod(ImageNumber,8)==6)
            FileData = imread(Fname, 'tif');
            FileData = reshape(FileData,[],1);
            FileData = double(FileData);
            bgndC2C2 = [bgndC2C2 FileData];
        elseif (ImageNumber<=length(F)/2 && mod(ImageNumber,8)==7)
            FileData = imread(Fname, 'tif');
            FileData = reshape(FileData,[],1);
            FileData = double(FileData);
            bgndC2L2 = [bgndC2L2 FileData];
        elseif (ImageNumber<=length(F)/2 && mod(ImageNumber,8)==0)
            FileData = imread(Fname, 'tif');
            FileData = reshape(FileData,[],1);
            FileData = double(FileData);
            bgndC2L1 = [bgndC2L1 FileData];   
        elseif (ImageNumber>length(F)/2 && mod(ImageNumber,8)==1)
            FileData = imread(Fname, 'tif');
            FileData = reshape(FileData,[],1);
            FileData = double(FileData);
            bgndL1C1 = [bgndL1C1 FileData];
        elseif (ImageNumber>length(F)/2 && mod(ImageNumber,8)==2)
            FileData = imread(Fname, 'tif');
            FileData = reshape(FileData,[],1);
            FileData = double(FileData);
            bgndL1C2 = [bgndL1C2 FileData];
        elseif (ImageNumber>length(F)/2 && mod(ImageNumber,8)==3)
            FileData = imread(Fname, 'tif');
            FileData = reshape(FileData,[],1);
            FileData = double(FileData);
            bgndL1L2 = [bgndL1L2 FileData];
        elseif (ImageNumber>length(F)/2 && mod(ImageNumber,8)==4)
            FileData = imread(Fname, 'tif');
            FileData = reshape(FileData,[],1);
            FileData = double(FileData);
            bgndL1L1 = [bgndL1L1 FileData];
        elseif (ImageNumber>length(F)/2 && mod(ImageNumber,8)==5)
            FileData = imread(Fname, 'tif');
            FileData = reshape(FileData,[],1);
            FileData = double(FileData);
            bgndL2C1 = [bgndL2C1 FileData];
        elseif (ImageNumber>length(F)/2 && mod(ImageNumber,8)==6)
            FileData = imread(Fname, 'tif');
            FileData = reshape(FileData,[],1);
            FileData = double(FileData);
            bgndL2C2 = [bgndL2C2 FileData];
        elseif (ImageNumber>length(F)/2 && mod(ImageNumber,8)==7)
            FileData = imread(Fname, 'tif');
            FileData = reshape(FileData,[],1);
            FileData = double(FileData);
            bgndL2L2 = [bgndL2L2 FileData];
        elseif (ImageNumber>length(F)/2 && mod(ImageNumber,8)==0)
            FileData = imread(Fname, 'tif');
            FileData = reshape(FileData,[],1);
            FileData = double(FileData);
            bgndL2L1 = [bgndL2L1 FileData];
        end
    end

    bgndMatC1C1=[bgndMatC1C1 bgndC1C1];bgndMatC2C1=[bgndMatC2C1 bgndC2C1];bgndMatL1C1=[bgndMatL1C1 bgndL1C1];bgndMatL2C1=[bgndMatL2C1 bgndL2C1];
    bgndMatC1C2=[bgndMatC1C2 bgndC1C2];bgndMatC2C2=[bgndMatC2C2 bgndC2C2];bgndMatL1C2=[bgndMatL1C2 bgndL1C2];bgndMatL2C2=[bgndMatL2C2 bgndL2C2];
    bgndMatC1L2=[bgndMatC1L2 bgndC1L2];bgndMatC2L2=[bgndMatC2L2 bgndC2L2];bgndMatL1L2=[bgndMatL1L2 bgndL1L2];bgndMatL2L2=[bgndMatL2L2 bgndL2L2];
    bgndMatC1L1=[bgndMatC1L1 bgndC1L1];bgndMatC2L1=[bgndMatC2L1 bgndC2L1];bgndMatL1L1=[bgndMatL1L1 bgndL1L1];bgndMatL2L1=[bgndMatL2L1 bgndL2L1];
end

%% %%%%%%%%%%%%%%%%%%%%%%%%% SHG Polarimetry %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize Data Matrices (You shouldn't have to change this)

ImageTag = [];
GroupTag = [];

CDcombined = [];
LDcombined = [];
Rcombined = [];
DCPcombined = [];
CIntcombined = []; 

for GroupNum = 1:length(GroupPaths)
    
    ImagePath1 = dir(GroupPaths{GroupNum});
    ImagePath2 = ImagePath1([ImagePath1(:).isdir]); 
    ImagePath2 = ImagePath2(~ismember({ImagePath2(:).name},{'.','..'}));
    ImagePaths = cell(1,length(ImagePath2));
    for i=1:length(ImagePath2)
        ImageName = ImagePath2(i).name;
        ImagePaths{i}=[GroupPaths{GroupNum} '\' ImageName];
    end

    for ImageNum = 1:length(ImagePaths)

        slashocc1 = strfind(ImagePaths{ImageNum},'\');
        slashocc2 = strfind(GroupPaths{GroupNum},'_');
        fprintf(['Performing polarimetric analysis on image ',ImagePaths{ImageNum}(slashocc1(end)+1:end),' in the ',GroupPaths{GroupNum}(slashocc2(end)+1:end),' group ...\n']);

        cd(char(ImagePaths(ImageNum)))

        % Variable Initialization (Notation: PSG-PSA-Wavelength)

        C1C1=[];C2C1=[];L1C1=[];L2C1=[];
        C1C2=[];C2C2=[];L1C2=[];L2C2=[];
        C1L2=[];C2L2=[];L1L2=[];L2L2=[];
        C1L1=[];C2L1=[];L1L1=[];L2L1=[];

        F = dir('*.tif');

        for iii = 1:length(F)
            Fname = F(iii).name; 
            ind1 = strfind(Fname,'_');
            ind2 = strfind(Fname,'.');
            ImageNumber = str2double(Fname(ind1(end)+1:ind2-1)); 
            if (ImageNumber<=length(F)/2 && mod(ImageNumber,8)==1)
                FileData = imread(Fname, 'tif');
                FileData = reshape(FileData,[],1);
                FileData = double(FileData);
                C1C1 = [C1C1 FileData];
            elseif (ImageNumber<=length(F)/2 && mod(ImageNumber,8)==2)
                FileData = imread(Fname, 'tif');
                FileData = reshape(FileData,[],1);
                FileData = double(FileData);
                C1C2 = [C1C2 FileData];
            elseif (ImageNumber<=length(F)/2 && mod(ImageNumber,8)==3)
                FileData = imread(Fname, 'tif');
                FileData = reshape(FileData,[],1);
                FileData = double(FileData);
                C1L2 = [C1L2 FileData];
            elseif (ImageNumber<=length(F)/2 && mod(ImageNumber,8)==4)
                FileData = imread(Fname, 'tif');
                FileData = reshape(FileData,[],1);
                FileData = double(FileData);
                C1L1 = [C1L1 FileData];
            elseif (ImageNumber<=length(F)/2 && mod(ImageNumber,8)==5)
                FileData = imread(Fname, 'tif');
                FileData = reshape(FileData,[],1);
                FileData = double(FileData);
                C2C1 = [C2C1 FileData];
            elseif (ImageNumber<=length(F)/2 && mod(ImageNumber,8)==6)
                FileData = imread(Fname, 'tif');
                FileData = reshape(FileData,[],1);
                FileData = double(FileData);
                C2C2 = [C2C2 FileData];
            elseif (ImageNumber<=length(F)/2 && mod(ImageNumber,8)==7)
                FileData = imread(Fname, 'tif');
                FileData = reshape(FileData,[],1);
                FileData = double(FileData);
                C2L2 = [C2L2 FileData];
            elseif (ImageNumber<=length(F)/2 && mod(ImageNumber,8)==0)
                FileData = imread(Fname, 'tif');
                FileData = reshape(FileData,[],1);
                FileData = double(FileData);
                C2L1 = [C2L1 FileData];   
            elseif (ImageNumber>length(F)/2 && mod(ImageNumber,8)==1)
                FileData = imread(Fname, 'tif');
                FileData = reshape(FileData,[],1);
                FileData = double(FileData);
                L1C1 = [L1C1 FileData];
            elseif (ImageNumber>length(F)/2 && mod(ImageNumber,8)==2)
                FileData = imread(Fname, 'tif');
                FileData = reshape(FileData,[],1);
                FileData = double(FileData);
                L1C2 = [L1C2 FileData];
            elseif (ImageNumber>length(F)/2 && mod(ImageNumber,8)==3)
                FileData = imread(Fname, 'tif');
                FileData = reshape(FileData,[],1);
                FileData = double(FileData);
                L1L2 = [L1L2 FileData];
            elseif (ImageNumber>length(F)/2 && mod(ImageNumber,8)==4)
                FileData = imread(Fname, 'tif');
                FileData = reshape(FileData,[],1);
                FileData = double(FileData);
                L1L1 = [L1L1 FileData];
            elseif (ImageNumber>length(F)/2 && mod(ImageNumber,8)==5)
                FileData = imread(Fname, 'tif');
                FileData = reshape(FileData,[],1);
                FileData = double(FileData);
                L2C1 = [L2C1 FileData];
            elseif (ImageNumber>length(F)/2 && mod(ImageNumber,8)==6)
                FileData = imread(Fname, 'tif');
                FileData = reshape(FileData,[],1);
                FileData = double(FileData);
                L2C2 = [L2C2 FileData];
            elseif (ImageNumber>length(F)/2 && mod(ImageNumber,8)==7)
                FileData = imread(Fname, 'tif');
                FileData = reshape(FileData,[],1);
                FileData = double(FileData);
                L2L2 = [L2L2 FileData];
            elseif (ImageNumber>length(F)/2 && mod(ImageNumber,8)==0)
                FileData = imread(Fname, 'tif');
                FileData = reshape(FileData,[],1);
                FileData = double(FileData);
                L2L1 = [L2L1 FileData];
            end
        end

        %% %%%%%%%%%%%%%%%%%%%%%%%% Averaging Image Data %%%%%%%%%%%%%%%%%%%%%%%%%%

        ImageName = str2num(ImagePaths{ImageNum}(slashocc1(end)+1:end));
        if any(ImageName == [1:15 17:26])==1
            bgndC1C1=bgndMatC1C1(:,1);bgndC2C1=bgndMatC2C1(:,1);
            bgndL1C1=bgndMatL1C1(:,1);bgndL2C1=bgndMatL2C1(:,1);
            bgndC1C2=bgndMatC1C2(:,1);bgndC2C2=bgndMatC2C2(:,1);
            bgndL1C2=bgndMatL1C2(:,1);bgndL2C2=bgndMatL2C2(:,1);
            bgndC1L2=bgndMatC1L2(:,1);bgndC2L2=bgndMatC2L2(:,1);
            bgndL1L2=bgndMatL1L2(:,1);bgndL2L2=bgndMatL2L2(:,1);
            bgndC1L1=bgndMatC1L1(:,1);bgndC2L1=bgndMatC2L1(:,1);
            bgndL1L1=bgndMatL1L1(:,1);bgndL2L1=bgndMatL2L1(:,1);
        else
            bgndC1C1=bgndMatC1C1(:,2);bgndC2C1=bgndMatC2C1(:,2);
            bgndL1C1=bgndMatL1C1(:,2);bgndL2C1=bgndMatL2C1(:,2);
            bgndC1C2=bgndMatC1C2(:,2);bgndC2C2=bgndMatC2C2(:,2);
            bgndL1C2=bgndMatL1C2(:,2);bgndL2C2=bgndMatL2C2(:,2);
            bgndC1L2=bgndMatC1L2(:,2);bgndC2L2=bgndMatC2L2(:,2);
            bgndL1L2=bgndMatL1L2(:,2);bgndL2L2=bgndMatL2L2(:,2);
            bgndC1L1=bgndMatC1L1(:,2);bgndC2L1=bgndMatC2L1(:,2);
            bgndL1L1=bgndMatL1L1(:,2);bgndL2L1=bgndMatL2L1(:,2);
        end

        % Average Background

        bgndC1C1=mean(bgndC1C1,2);bgndC2C1=mean(bgndC2C1,2);
        bgndL1C1=mean(bgndL1C1,2);bgndL2C1=mean(bgndL2C1,2);
        bgndC1C2=mean(bgndC1C2,2);bgndC2C2=mean(bgndC2C2,2);
        bgndL1C2=mean(bgndL1C2,2);bgndL2C2=mean(bgndL2C2,2);
        bgndC1L2=mean(bgndC1L2,2);bgndC2L2=mean(bgndC2L2,2);
        bgndL1L2=mean(bgndL1L2,2);bgndL2L2=mean(bgndL2L2,2);
        bgndC1L1=mean(bgndC1L1,2);bgndC2L1=mean(bgndC2L1,2);
        bgndL1L1=mean(bgndL1L1,2);bgndL2L1=mean(bgndL2L1,2);

        % Average Images 

        C1C1=mean(C1C1,2);C2C1=mean(C2C1,2);
        L1C1=mean(L1C1,2);L2C1=mean(L2C1,2);
        C1C2=mean(C1C2,2);C2C2=mean(C2C2,2);
        L1C2=mean(L1C2,2);L2C2=mean(L2C2,2);
        C1L2=mean(C1L2,2);C2L2=mean(C2L2,2);
        L1L2=mean(L1L2,2);L2L2=mean(L2L2,2);
        C1L1=mean(C1L1,2);C2L1=mean(C2L1,2);
        L1L1=mean(L1L1,2);L2L1=mean(L2L1,2);

        % Subtract Background

        C1C1=C1C1-bgndC1C1;C2C1=C2C1-bgndC2C1;
        L1C1=L1C1-bgndL1C1;L2C1=L2C1-bgndL2C1;
        C1C2=C1C2-bgndC1C2;C2C2=C2C2-bgndC2C2;
        L1C2=L1C2-bgndL1C2;L2C2=L2C2-bgndL2C2;
        C1L2=C1L2-bgndC1L2;C2L2=C2L2-bgndC2L2;
        L1L2=L1L2-bgndL1L2;L2L2=L2L2-bgndL2L2;
        C1L1=C1L1-bgndC1L1;C2L1=C2L1-bgndC2L1;
        L1L1=L1L1-bgndL1L1;L2L1=L2L1-bgndL2L1;

        % Remove Negative Intensity

        C1C1(C1C1<0)=0;C2C1(C2C1<0)=0;L1C1(L1C1<0)=0;L2C1(L2C1<0)=0;
        C1C2(C1C2<0)=0;C2C2(C2C2<0)=0;L1C2(L1C2<0)=0;L2C2(L2C2<0)=0;
        C1L2(C1L2<0)=0;C2L2(C2L2<0)=0;L1L2(L1L2<0)=0;L2L2(L2L2<0)=0;
        C1L1(C1L1<0)=0;C2L1(C2L1<0)=0;L1L1(L1L1<0)=0;L2L1(L2L1<0)=0; 

        %% %%%%%%%%%%%%%%%%%%%%%% Pixel Intensity Cutoff %%%%%%%%%%%%%%%%%%%%%%%%%%

        PixelValueCutoff = 0.107786*(SNRcutoff+1.8003)^2.98597;


        C1C1Final = C1C1; C1C1Final(C1C1<PixelValueCutoff)=nan;
        C2C1Final = C2C1; C2C1Final(C2C1<PixelValueCutoff)=nan;
        L1C1Final = L1C1; L1C1Final(L1C1<PixelValueCutoff)=nan;
        L2C1Final = L2C1; L2C1Final(L2C1<PixelValueCutoff)=nan;

        C1C2Final = C1C2; C1C2Final(C1C2<PixelValueCutoff)=nan;
        C2C2Final = C2C2; C2C2Final(C2C2<PixelValueCutoff)=nan;
        L1C2Final = L1C2; L1C2Final(L1C2<PixelValueCutoff)=nan;
        L2C2Final = L2C2; L2C2Final(L2C2<PixelValueCutoff)=nan;

        C1L2Final = C1L2; C1L2Final(C1L2<PixelValueCutoff)=nan;
        C2L2Final = C2L2; C2L2Final(C2L2<PixelValueCutoff)=nan;
        L1L2Final = L1L2; L1L2Final(L1L2<PixelValueCutoff)=nan;
        L2L2Final = L2L2; L2L2Final(L2L2<PixelValueCutoff)=nan;

        C1L1Final = C1L1; C1L1Final(C1L1<PixelValueCutoff)=nan;
        C2L1Final = C2L1; C2L1Final(C2L1<PixelValueCutoff)=nan;
        L1L1Final = L1L1; L1L1Final(L1L1<PixelValueCutoff)=nan;
        L2L1Final = L2L1; L2L1Final(L2L1<PixelValueCutoff)=nan;

        %% %%%%%%%%%%%%%%%%%%%% Stokes Vector Calculations %%%%%%%%%%%%%%%%%%%%%%%%

        % C1 = LCP       (Equation: s0/2 + s3/2)
        % C2 = RCP       (Equation: s0/2 - s3/2)
        % L1 = FullWave  (Equation: s0/2 + s1/2)
        % L2 = HalfWave  (Equation: s0/2 - s1/2)

        % For each PSG polarization input, the SHG Stokes vector is calculated.
        % We have 4 equations and 3 unknowns, we then have 4 different combinations that are averaged.

        %%%%% Combo 1: C2, C1, L1
        EqSystemC2C1L1 = [1,0,-1;1,0,1;1,1,0]; % only contains s0, s1, and s3

        % C1 Input
        C1MatrixC2C1L1 = vertcat(C1C2Final',C1C1Final',C1L1Final');
        C1sVecC2C1L1 = EqSystemC2C1L1\C1MatrixC2C1L1;
        C1s0C2C1L1 = C1sVecC2C1L1(1,:);C1s1C2C1L1 = C1sVecC2C1L1(2,:);C1s3C2C1L1 = C1sVecC2C1L1(3,:);

        % C2 Input
        C2MatrixC2C1L1 = vertcat(C2C2Final',C2C1Final',C2L1Final');
        C2sVecC2C1L1 = EqSystemC2C1L1\C2MatrixC2C1L1;
        C2s0C2C1L1 = C2sVecC2C1L1(1,:);C2s1C2C1L1 = C2sVecC2C1L1(2,:);C2s3C2C1L1 = C2sVecC2C1L1(3,:);

        % L1 Input
        L1MatrixC2C1L1 = vertcat(L1C2Final',L1C1Final',L1L1Final');
        L1sVecC2C1L1 = EqSystemC2C1L1\L1MatrixC2C1L1;
        L1s0C2C1L1 = L1sVecC2C1L1(1,:);L1s1C2C1L1 = L1sVecC2C1L1(2,:);L1s3C2C1L1 = L1sVecC2C1L1(3,:);

        % L2 Input
        L2MatrixC2C1L1 = vertcat(L2C2Final',L2C1Final',L2L1Final');
        L2sVecC2C1L1 = EqSystemC2C1L1\L2MatrixC2C1L1;
        L2s0C2C1L1 = L2sVecC2C1L1(1,:);L2s1C2C1L1 = L2sVecC2C1L1(2,:);L2s3C2C1L1 = L2sVecC2C1L1(3,:);

        %%%%% Combo 2: C2, C1, L2
        EqSystemC2C1L2 = [1,0,-1;1,0,1;1,-1,0]; % only contains s0, s1, and s3

        % C1 Input
        C1MatrixC2C1L2 = vertcat(C1C2Final',C1C1Final',C1L2Final');
        C1sVecC2C1L2 = EqSystemC2C1L2\C1MatrixC2C1L2;
        C1s0C2C1L2 = C1sVecC2C1L2(1,:);C1s1C2C1L2 = C1sVecC2C1L2(2,:);C1s3C2C1L2 = C1sVecC2C1L2(3,:);

        % C2 Input
        C2MatrixC2C1L2 = vertcat(C2C2Final',C2C1Final',C2L2Final');
        C2sVecC2C1L2 = EqSystemC2C1L2\C2MatrixC2C1L2;
        C2s0C2C1L2 = C2sVecC2C1L2(1,:);C2s1C2C1L2 = C2sVecC2C1L2(2,:);C2s3C2C1L2 = C2sVecC2C1L2(3,:);

        % L1 Input
        L1MatrixC2C1L2 = vertcat(L1C2Final',L1C1Final',L1L2Final');
        L1sVecC2C1L2 = EqSystemC2C1L2\L1MatrixC2C1L2;
        L1s0C2C1L2 = L1sVecC2C1L2(1,:);L1s1C2C1L2 = L1sVecC2C1L2(2,:);L1s3C2C1L2 = L1sVecC2C1L2(3,:);

        % L2 Input
        L2MatrixC2C1L2 = vertcat(L2C2Final',L2C1Final',L2L2Final');
        L2sVecC2C1L2 = EqSystemC2C1L2\L2MatrixC2C1L2;
        L2s0C2C1L2 = L2sVecC2C1L2(1,:);L2s1C2C1L2 = L2sVecC2C1L2(2,:);L2s3C2C1L2 = L2sVecC2C1L2(3,:);

        %%%%% Combo 3: C2, L1, L2
        EqSystemC2L1L2 = [1,0,-1;1,1,0;1,-1,0]; % only contains s0, s1, and s3

        % C1 Input
        C1MatrixC2L1L2 = vertcat(C1C2Final',C1L1Final',C1L2Final');
        C1sVecC2L1L2 = EqSystemC2L1L2\C1MatrixC2L1L2;
        C1s0C2L1L2 = C1sVecC2L1L2(1,:);C1s1C2L1L2 = C1sVecC2L1L2(2,:);C1s3C2L1L2 = C1sVecC2L1L2(3,:);

        % C2 Input
        C2MatrixC2L1L2 = vertcat(C2C2Final',C2L1Final',C2L2Final');
        C2sVecC2L1L2 = EqSystemC2L1L2\C2MatrixC2L1L2;
        C2s0C2L1L2 = C2sVecC2L1L2(1,:);C2s1C2L1L2 = C2sVecC2L1L2(2,:);C2s3C2L1L2 = C2sVecC2L1L2(3,:);

        % L1 Input
        L1MatrixC2L1L2 = vertcat(L1C2Final',L1L1Final',L1L2Final');
        L1sVecC2L1L2 = EqSystemC2L1L2\L1MatrixC2L1L2;
        L1s0C2L1L2 = L1sVecC2L1L2(1,:);L1s1C2L1L2 = L1sVecC2L1L2(2,:);L1s3C2L1L2 = L1sVecC2L1L2(3,:);

        % L2 Input
        L2MatrixC2L1L2 = vertcat(L2C2Final',L2L1Final',L2L2Final');
        L2sVecC2L1L2 = EqSystemC2L1L2\L2MatrixC2L1L2;
        L2s0C2L1L2 = L2sVecC2L1L2(1,:);L2s1C2L1L2 = L2sVecC2L1L2(2,:);L2s3C2L1L2 = L2sVecC2L1L2(3,:);

        %%%%% Combo 4: C1, L1, L2
        EqSystemC1L1L2 = [1,0,1;1,1,0;1,-1,0]; % only contains s0, s1, and s3

        % C1 Input
        C1MatrixC1L1L2 = vertcat(C1C1Final',C1L1Final',C1L2Final');
        C1sVecC1L1L2 = EqSystemC1L1L2\C1MatrixC1L1L2;
        C1s0C1L1L2 = C1sVecC1L1L2(1,:);C1s1C1L1L2 = C1sVecC1L1L2(2,:);C1s3C1L1L2 = C1sVecC1L1L2(3,:);

        % C2 Input
        C2MatrixC1L1L2 = vertcat(C2C1Final',C2L1Final',C2L2Final');
        C2sVecC1L1L2 = EqSystemC1L1L2\C2MatrixC1L1L2;
        C2s0C1L1L2 = C2sVecC1L1L2(1,:);C2s1C1L1L2 = C2sVecC1L1L2(2,:);C2s3C1L1L2 = C2sVecC1L1L2(3,:);

        % L1 Input
        L1MatrixC1L1L2 = vertcat(L1C1Final',L1L1Final',L1L2Final');
        L1sVecC1L1L2 = EqSystemC1L1L2\L1MatrixC1L1L2;
        L1s0C1L1L2 = L1sVecC1L1L2(1,:);L1s1C1L1L2 = L1sVecC1L1L2(2,:);L1s3C1L1L2 = L1sVecC1L1L2(3,:);

        % L2 Input
        L2MatrixC1L1L2 = vertcat(L2C1Final',L2L1Final',L2L2Final');
        L2sVecC1L1L2 = EqSystemC1L1L2\L2MatrixC1L1L2;
        L2s0C1L1L2 = L2sVecC1L1L2(1,:);L2s1C1L1L2 = L2sVecC1L1L2(2,:);L2s3C1L1L2 = L2sVecC1L1L2(3,:);

        % Averaged Stokes Vector Elements

        C1s0 = mean([C1s0C1L1L2;C1s0C2L1L2;C1s0C2C1L2;C1s0C2C1L1]);
        C1s1 = mean([C1s1C1L1L2;C1s1C2L1L2;C1s1C2C1L2;C1s1C2C1L1]);
        C1s3 = mean([C1s3C1L1L2;C1s3C2L1L2;C1s3C2C1L2;C1s3C2C1L1]);

        C2s0 = mean([C2s0C1L1L2;C2s0C2L1L2;C2s0C2C1L2;C2s0C2C1L1]);
        C2s1 = mean([C2s1C1L1L2;C2s1C2L1L2;C2s1C2C1L2;C2s1C2C1L1]);
        C2s3 = mean([C2s3C1L1L2;C2s3C2L1L2;C2s3C2C1L2;C2s3C2C1L1]);

        L1s0 = mean([L1s0C1L1L2;L1s0C2L1L2;L1s0C2C1L2;L1s0C2C1L1]);
        L1s1 = mean([L1s1C1L1L2;L1s1C2L1L2;L1s1C2C1L2;L1s1C2C1L1]);
        L1s3 = mean([L1s3C1L1L2;L1s3C2L1L2;L1s3C2C1L2;L1s3C2C1L1]);

        L2s0 = mean([L2s0C1L1L2;L2s0C2L1L2;L2s0C2C1L2;L2s0C2C1L1]);
        L2s1 = mean([L2s1C1L1L2;L2s1C2L1L2;L2s1C2C1L2;L2s1C2C1L1]);
        L2s3 = mean([L2s3C1L1L2;L2s3C2L1L2;L2s3C2C1L2;L2s3C2C1L1]);

        %% %%%%%%%%%%%%%%%% Polarimteric Parameter Calculations %%%%%%%%%%%%%%%%%%%

        CD = reshape(2*(C1s0-C2s0)./(C1s0+C2s0),sqrt(length(C1s0)),sqrt(length(C1s0)));
        CDabs = abs(CD);
        LD = reshape(2*(L1s0-L2s0)./(L1s0+L2s0),sqrt(length(L1s0)),sqrt(length(L1s0)));
        LDabs = abs(LD);
        ACC = reshape((C1s3-C2s3)./(C1s0+C2s0),sqrt(length(C1s0)),sqrt(length(C1s0)));
        Rratio = ACC;

        % Uses the positive R equation for ACC less than 0, and negative R equation for ACC larger than 0
        for i = 1:size(ACC,1)
            for j = 1:size(ACC,2)
                if ACC(i,j)<0
                    Rratio(i,j) = real(1+2./ACC(i,j)+2*sqrt(1./ACC(i,j).^2-1));
                elseif ACC(i,j)>0
                    Rratio(i,j) = real(1+2./ACC(i,j)-2*sqrt(1./ACC(i,j).^2-1));
                elseif ACC(i,j)==0
                    Rratio(i,j)=1;
                end
            end
        end

        DCP1 = reshape(abs(C1s3)./C1s0,sqrt(length(C1s0)),sqrt(length(C1s0)));
        DCP2 = reshape(abs(C2s3)./C2s0,sqrt(length(C2s0)),sqrt(length(C2s0)));
        DCP = (DCP1+DCP2)/2;

        AllCircInt = [C1C1Final,C1C2Final,C1L2Final,C1L1Final,C2C1Final,C2C2Final,C2L2Final,C2L1Final];
        CInt = mean(AllCircInt,2);
        CInt = reshape(CInt,sqrt(length(CInt)),[]);
        CInt1 = (1/2)*(C1s0+C2s0);
        CInt1 = reshape(CInt1,sqrt(length(CInt1)),[]);

        % Pixel Matching of Different Parameters
        CD(isnan(LD)==1)=nan;CD(isnan(Rratio)==1)=nan;CD(isnan(CInt)==1)=nan;CD(isnan(DCP)==1)=nan;
        LD(isnan(CD)==1)=nan;LD(isnan(Rratio)==1)=nan;LD(isnan(CInt)==1)=nan;LD(isnan(DCP)==1)=nan;
        Rratio(isnan(LD)==1)=nan;Rratio(isnan(CD)==1)=nan;Rratio(isnan(CInt)==1)=nan;Rratio(isnan(DCP)==1)=nan;
        CInt(isnan(LD)==1)=nan;CInt(isnan(CD)==1)=nan;CInt(isnan(Rratio)==1)=nan;CInt(isnan(DCP)==1)=nan;
        DCP(isnan(LD)==1)=nan;DCP(isnan(Rratio)==1)=nan;DCP(isnan(CInt)==1)=nan;DCP(isnan(CD)==1)=nan;
        CInt1(isnan(LD)==1)=nan;CInt1(isnan(CD)==1)=nan;CInt1(isnan(Rratio)==1)=nan;CInt1(isnan(DCP)==1)=nan;

        %% %%%%%%%%%%%%%%%%%%%%% Scan Name Initialization %%%%%%%%%%%%%%%%%%%%%%%%%

        % Folder Name
        CurrentFpath = pwd;
        SlashOccurances = strfind(CurrentFpath,'\');
        CurrentFolderName = CurrentFpath(SlashOccurances(end)+1:end);

        % Image Names
        CDImageName=sprintf('CD_%s.tiff',CurrentFolderName);
        LDImageName=sprintf('LD_%s.tiff',CurrentFolderName);
        RImageName=sprintf('R_%s.tiff',CurrentFolderName);
        DCPImageName=sprintf('DCP_%s.tiff',CurrentFolderName);
        PixelDensityImageName=sprintf('PixelDensity_%s.tiff',CurrentFolderName);
        LinearIntImageSameScaleName=sprintf('LinearIntensity_%s.tiff',CurrentFolderName);
        LogIntBWImageSameScaleName=sprintf('BWLogIntensity_%s.tiff',CurrentFolderName);
        LogIntColoredImageSameScaleName=sprintf('ColoredLogIntensity_%s.tiff',CurrentFolderName);

        % Histogram Names
        CDHistName=sprintf('CDhist_%s.tiff',CurrentFolderName);
        LDHistName=sprintf('LDhist_%s.tiff',CurrentFolderName);
        RHistName=sprintf('Rhist_%s.tiff',CurrentFolderName);
        DCPHistName=sprintf('DCPhist_%s.tiff',CurrentFolderName);
        Int2HistName=sprintf('INThistNORMLOG_%s.tiff',CurrentFolderName);
        Int1HistName=sprintf('INThistLOGLOG_%s.tiff',CurrentFolderName);

        %% %%%%%%%%%%%%%%%%%%% Parameter Matrix Concatenation %%%%%%%%%%%%%%%%%%%%%

        ImageTag = [ImageTag str2num(CurrentFolderName)];
        GroupTag = [GroupTag GroupNumber(GroupNum)];
        CDcombined = [CDcombined reshape(CD,[],1)];
        LDcombined = [LDcombined reshape(LD,[],1)];
        Rcombined = [Rcombined reshape(Rratio,[],1)];
        DCPcombined = [DCPcombined reshape(DCP,[],1)];
        CIntcombined = [CIntcombined reshape(CInt,[],1)];

        %% %%%%%%%%%%%%%%%%%%%%%%% Plot Parameter Images %%%%%%%%%%%%%%%%%%%%%%%%%%

        % Color Bar Initialization

        LDcmap = [0,0.2,1;0,0.4,1;0,0.6,1;0,0.8,1;0,0.95,1;1,1,1;1,0.95,0;1,0.8,0;1,0.6,0;1,0.4,0;1,0.2,0];
        CDcmap = LDcmap.^3;
        Rcmap =  [0.171,0.307,1;0.245,0.592,1;0.002,1,1;0.002,1,0.539;0.186,1,0.186;0.829,1,0.002;0.971,0.971,0.002;1,0.7,0.002;1,0.406,0.037;0.951,0.002,0.002;0.759,0.002,0.327];
        DCPcmap = spring;

        if DisplayImages == 1
            set(groot,'DefaultFigureVisible','off')
            figure('Name',CDImageName);
            HeatmapMaker(CD,2048,CDcmap,minCD,maxCD)

            figure('Name',LDImageName);
            HeatmapMaker(LD,2048,LDcmap,minLD,maxLD)

            figure('Name',RImageName);
            HeatmapMaker(Rratio,2048,Rcmap,minR,maxR)

            figure('Name',DCPImageName);
            HeatmapMaker(DCP,2048,DCPcmap,minDCP,maxDCP)

            figure('Name',PixelDensityImageName);
            HeatmapMaker(CInt,2048,gray,0,1)

            figure('Name',LinearIntImageSameScaleName);
            HeatmapMaker(CInt,2048,gray,MinSHGInt,MaxSHGInt)

            figure('Name',LogIntColoredImageSameScaleName);
            HeatmapMaker(log10(CInt),2048,jet,log10(MinSHGInt),log10(MaxSHGInt))

            figure('Name',LogIntBWImageSameScaleName);
            HeatmapMaker(log10(CInt),2048,gray,log10(MinSHGInt),log10(MaxSHGInt))
            set(groot,'DefaultFigureVisible','on');
            Allfigs = findall(groot,'Type','figure');
            for i=1:length(Allfigs)
                Allfigs(i).Visible='on';
            end
        end  

        %% %%%%%%%%%%%%%%%%%%%%% Plot Parameter Histograms %%%%%%%%%%%%%%%%%%%%%%%%

        if DisplayHistograms == 1
            set(groot,'DefaultFigureVisible','off')
            % CD
            [CDbinfreq,CDbinval] = histcounts(CD,'BinMethod','fd');
            CDbinval = CDbinval(1:end-1) + diff(CDbinval)/2;

            figure('Name',CDHistName);
            plot(CDbinval,CDbinfreq./max(max(CDbinfreq)),'k','LineWidth',2)
            xlabel('SHG-CD') 
            ylabel('Normalized Occurrence')
            set(gcf, 'Color', 'w');
            set(gca,'FontSize',18);
            pbaspect([1 1 1])
            axis([minCD maxCD 0 1])
            ax = gca;
            ax.XTick = [minCD minCD/2 0 maxCD/2 maxCD];

            % LD
            [LDbinfreq,LDbinval] = histcounts(LD,'BinMethod','fd');
            LDbinval = LDbinval(1:end-1) + diff(LDbinval)/2;

            figure('Name',LDHistName);
            plot(LDbinval,LDbinfreq./max(max(LDbinfreq)),'k','LineWidth',2)
            xlabel('SHG-LD') 
            ylabel('Normalized Occurrence')
            set(gcf, 'Color', 'w');
            set(gca,'FontSize',18);
            pbaspect([1 1 1])
            axis([minLD maxLD 0 1])
            ax = gca;
            ax.XTick = [minLD minLD/2 0 maxLD/2 maxLD];

            % R
            [Rbinfreq,Rbinval] = histcounts(Rratio,'BinMethod','fd');
            Rbinval = Rbinval(1:end-1) + diff(Rbinval)/2;

            figure('Name',RHistName);
            plot(Rbinval,Rbinfreq./max(max(Rbinfreq)),'k','LineWidth',2)
            xlabel('R-Ratio') 
            ylabel('Normalized Occurrence')
            set(gcf, 'Color', 'w');
            set(gca,'FontSize',18);
            pbaspect([1 1 1])
            axis([minR maxR 0 1])
            ax = gca;
            ax.XTick = minR:0.2:maxR;

            % DCP
            [DCPbinfreq,DCPbinval] = histcounts(DCP,'BinMethod','fd');
            DCPbinval = DCPbinval(1:end-1) + diff(DCPbinval)/2;

            figure('Name',DCPHistName);
            plot(DCPbinval,DCPbinfreq./max(max(DCPbinfreq)),'k','LineWidth',2)
            xlabel('Degree of Circular Polarization') 
            ylabel('Normalized Occurrence')
            set(gcf, 'Color', 'w');
            set(gca,'FontSize',18);
            pbaspect([1 1 1])
            axis([minDCP maxDCP 0 1])
            ax = gca;
            ax.XTick = minDCP:0.1:maxDCP;

            % Circular Intensity
            [CIntbinfreq,CIntbinval] = histcounts(log10(CInt),'BinMethod','fd');
            CIntbinval = CIntbinval(1:end-1) + diff(CIntbinval)/2;

            figure('Name',Int2HistName);
            loglog(10.^CIntbinval,CIntbinfreq,'k','LineWidth',2)
            xlabel('Intensity') 
            ylabel('Logarithmic Occurrence')
            set(gcf, 'Color', 'w');
            set(gca,'FontSize',18);
            pbaspect([1 1 1])
            axis([min(10.^CIntbinval) max(10.^CIntbinval) 0 Inf])
            ax = gca;
            ax.XTick = [1 10 100 1000 10000 100000];
            ax.YTick = [1 10 100 1000 10000 100000];

            figure('Name',Int1HistName);
            semilogx(10.^CIntbinval,CIntbinfreq./max(max(CIntbinfreq)),'k','LineWidth',2)
            xlabel('Intensity') 
            ylabel('Normalized Occurrence')
            set(gcf, 'Color', 'w');
            set(gca,'FontSize',18);
            pbaspect([1 1 1])
            axis([MinSHGInt MaxSHGInt 0 Inf])
            ax = gca;
            ax.XTick = [2 20 100 1000 10000];
        end
        set(groot,'DefaultFigureVisible','on');
        Allfigs = findall(groot,'Type','figure');
        for i=1:length(Allfigs)
            Allfigs(i).Visible='on';
        end
    end
end

fprintf('Polarimetric Analysis Complete!\n ----------------------------- \n')
clearvars -except ImageTag GroupTag CDcombined LDcombined Rcombined DCPcombined CIntcombined numOfSubImagesVec d MultipleComparisons TextureAnalysis Classification Boxplots

%% %%%%%%%%%%%%%%%%%%%%%%%%%% Texture Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%

if TextureAnalysis==1

    BoundsWithNgMat = [28.1118250000000,6990.98857500000,551.080000000000,1.15446760100000,2.93192463400000,526.360000000000,0.204406967000000,1.00314596900000,256.200000000000,-1.30466150700000,1.31916167200000,368.660000000000,-1.48694373800000,1.45630734000000,127;...
                    25.4881000000000,7320.27007500000,330.420000000000,1.13839967000000,2.98020291400000,379,0.190913765000000,1.01588850400000,232.220000000000,-1.33393594800000,1.35486755100000,184.610000000000,-1.50365614700000,1.50868326200000,105.060000000000;...
                    24.0510375000000,8715.38700000000,260.710000000000,1.06773819400000,2.98723743100000,177.260000000000,0.143326847000000,1.02289010600000,116,-1.38490485000000,1.40965903300000,122,-1.55380778900000,1.54214911400000,75;...
                    22.9785250000000,7458.31837500000,172.740000000000,0.996002434000000,2.99142220600000,124.740000000000,0.109720336000000,1.03471726000000,109,-1.43107587800000,1.45706883800000,108,-1.59706447600000,1.58374931900000,68;...
                    22.1250000000000,6477.88550000000,115,0.934585072000000,2.99346692500000,82,0.0919467720000000,1.04265046600000,54,-1.45773309100000,1.48485594900000,69,-1.62662288500000,1.59145249800000,53;...
                    22.3396000000000,5574.30220000000,68,0.921623898000000,2.99389782900000,55,0.0840917320000000,1.03842372000000,48,-1.47527037900000,1.51296577400000,51,-1.63651982100000,1.61587142600000,35;...
                    22.5764500000000,5049.81707500000,38,0.871701537000000,2.99419564100000,39,0.0739499470000000,1.03902800800000,33,-1.51920778800000,1.53299664000000,31,-1.64828660100000,1.62612413200000,31;...
                    23,4266.12267500000,20,0.844337113000000,2.99426994800000,24,0.0650000000000000,1.03873299000000,24,-1.54219242700000,1.55366748500000,22,-1.64412008500000,1.62980573100000,23;...
                    24.8750000000000,3633.09915000000,12,0.893695322000000,2.99187267700000,15,0.0772988850000000,1.02751436600000,15,-1.51361804700000,1.52312228800000,12,-1.60624814100000,1.58437285800000,14];


    SubINg = [];
    meanPD = [];
    meanFOS = [];
    meanFeat = [];
    madPD = [];
    madFOS = [];
    madFeat = [];
    SubimageRuntime = [];
    AllSubimageResults = cell(1,length(numOfSubImagesVec));

    for i = 1:length(numOfSubImagesVec)
        tic
        numOfSubImages = numOfSubImagesVec(i);

        fprintf(['Performing statistical and texture analysis with ',num2str(numOfSubImages),' subimages ...\n']);

        if numOfSubImages == 1
            BoundsWithNg = BoundsWithNgMat(1,:);
        elseif numOfSubImages == 4
            BoundsWithNg = BoundsWithNgMat(2,:);
        elseif numOfSubImages == 16
            BoundsWithNg = BoundsWithNgMat(3,:);
        elseif numOfSubImages == 64
            BoundsWithNg = BoundsWithNgMat(4,:);
        elseif numOfSubImages == 256
            BoundsWithNg = BoundsWithNgMat(5,:);
        elseif numOfSubImages == 1024
            BoundsWithNg = BoundsWithNgMat(6,:);
        elseif numOfSubImages == 4096
            BoundsWithNg = BoundsWithNgMat(7,:);
        elseif numOfSubImages == 16384
            BoundsWithNg = BoundsWithNgMat(8,:);
        elseif numOfSubImages == 65536
            BoundsWithNg = BoundsWithNgMat(9,:);
        end

        LbINT = BoundsWithNg(1);
        UbINT = BoundsWithNg(2);
        NgINT = BoundsWithNg(3);

        LbR = BoundsWithNg(4);
        UbR = BoundsWithNg(5);
        NgR = BoundsWithNg(6);

        LbDCP = BoundsWithNg(7);
        UbDCP = BoundsWithNg(8);
        NgDCP = BoundsWithNg(9);

        LbCD = BoundsWithNg(10);
        UbCD = BoundsWithNg(11);
        NgCD = BoundsWithNg(12);

        LbLD = BoundsWithNg(13);
        UbLD = BoundsWithNg(14);
        NgLD = BoundsWithNg(15);



        AllParamMatrices  = cell(1,5);
        AllParamMatrices{1} = CDcombined;
        AllParamMatrices{2} = DCPcombined;
        AllParamMatrices{3} = CIntcombined;
        AllParamMatrices{4} = LDcombined;
        AllParamMatrices{5} = Rcombined;

        % Initialize data arrays
        groupTag = [];      % Group tag, eg. normal is group 1
        coreNum = [];       % Breast core number
        PD      = [];       % Pixel density
        Feat = [];          % Textural features
        FOS = [];           % First order statistics

        ParamFeat = [];   % Textural features of one parameter
        ParamFOS = [];    % FOS                   "       "   

        for j = 1:length(AllParamMatrices)
            Iinter = AllParamMatrices{j};

            % Initialize results array
            ImageSubFeat = [];  
            ImageFOS = [];     

            for k = 1:size(Rcombined,2)

                % Reshape Data and Set Percentile Bounds
                I = Iinter(:,k);
                % Reshape data into image (for image division purposes)
                len = sqrt(size(I,1));
                I = reshape(I,  [len, len])';

                if j==1
                    imMin = LbCD; imMax = UbCD; Ng = round(NgCD);
                elseif j==2
                    imMin = LbDCP; imMax = UbDCP; Ng = round(NgDCP);
                elseif j==3
                    imMin = LbINT; imMax = UbINT; Ng = round(NgINT);
                elseif j==4
                    imMin = LbLD; imMax = UbLD; Ng = round(NgLD);
                elseif j==5
                    imMin = LbR; imMax = UbR; Ng = round(NgR);
                end

                % Divide into subimages
                blocks = DivideImage(I, numOfSubImages);
                blocks = reshape(blocks, [numOfSubImages,1]);

                if j==1
                % Get group tag and image number
                    if GroupTag(k) == 1
                        groupTag = [groupTag; 1*ones(numOfSubImages,1) GroupTag(k)*ones(numOfSubImages,1) [1:numOfSubImages]'];
                    elseif GroupTag(k) == 2
                        groupTag = [groupTag; 2*ones(numOfSubImages,1) GroupTag(k)*ones(numOfSubImages,1) [1:numOfSubImages]'];
                    elseif GroupTag(k) == 3
                        groupTag = [groupTag; 2*ones(numOfSubImages,1) GroupTag(k)*ones(numOfSubImages,1) [1:numOfSubImages]'];
                    elseif GroupTag(k) == 4
                        groupTag = [groupTag; 2*ones(numOfSubImages,1) GroupTag(k)*ones(numOfSubImages,1) [1:numOfSubImages]'];
                    end


                    coreNum = [coreNum; ImageTag(k)*ones(numOfSubImages,1)]; 

                    for im = 1:numOfSubImages
                        PD = [PD; nnz(~isnan(blocks{im}))];
                    end
                end

                % Perform texture analysis on subimages
                for m=1:length(blocks)
                    if nnz(~isnan(blocks{m}))>0
                        ImageFOS = [ImageFOS; nanmean(blocks{m}(:)) nanmedian(blocks{m}(:)) ...
                            nanstd(blocks{m}(:)) mad(blocks{m}(:)) mad(blocks{m}(:),1)];
                        ImageSubFeat = [ImageSubFeat; GLCMFeat(imMin,imMax,Ng,blocks{m},d)];
                    else
                        ImageFOS = [ImageFOS; nan nan nan nan nan];
                        ImageSubFeat = [ImageSubFeat; nan nan nan nan nan];
                    end
                end 

            end
            % Update array of parameters
            ParamFeat = [ParamFeat ImageSubFeat];
            ParamFOS = [ParamFOS ImageFOS];

        end
        Feat = [Feat; ParamFeat];
        FOS = [FOS; ParamFOS];

        Matrix = [groupTag coreNum PD Feat(:,1:5:25) Feat(:,2:5:25) Feat(:,3:5:25) Feat(:,4:5:25) Feat(:,5:5:25) FOS];
        AllSubimageResults{i} = Matrix;

        timelapse = toc;
        SubimageRuntime = [SubimageRuntime timelapse];

    end

    fprintf('Texture Analysis Complete!\n ----------------------------- \n')
    clearvars -except ImageTag GroupTag CDcombined LDcombined Rcombined DCPcombined CIntcombined numOfSubImagesVec AllParamMatrices AllSubimageResults SubimageRuntime MultipleComparisons Classification Boxplots
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%% Classification %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Classification==1
    % Defining the name of each column
    Xnametag = {'2-Group ID','4-Group ID','Subimage Number','Core Number','Pixel Density',...
        'SHG-CD Contrast','DCP Contrast','SHG Intensity Contrast','SHG-LD Contrast','R-Ratio Contrast',...
        'SHG-CD Correlation','DCP Correlation','SHG Intensity Correlation','SHG-LD Correlation','R-Ratio Correlation',...
        'SHG-CD Entropy','DCP Entropy','SHG Intensity Entropy','SHG-LD Entropy','R-Ratio Entropy',...
        'SHG-CD ASM','DCP ASM','SHG Intensity ASM','SHG-LD ASM','R-Ratio ASM',...
        'SHG-CD IDM','DCP IDM','SHG Intensity IDM','SHG-LD IDM','R-Ratio IDM',...
        'SHG-CD Mean','SHG-CD Median','SHG-CD Standard Deviation','SHG-CD Mean Absolute Deviation','SHG-CD Median Absolute Deviation',...
        'DCP Mean','DCP Median','DCP Standard Deviation','DCP Mean Absolute Deviation','DCP Median Absolute Deviation',...
        'SHG Intensity Mean','SHG Intensity Median','SHG Intensity Standard Deviation','SHG Intensity Mean Absolute Deviation','SHG Intensity Median Absolute Deviation',...
        'SHG-LD Mean','SHG-LD Median','SHG-LD Standard Deviation','SHG-LD Mean Absolute Deviation','SHG-LD Median Absolute Deviation',...
        'R-Ratio Mean','R-Ratio Median','R-Ratio Standard Deviation','R-Ratio Mean Absolute Deviation','R-Ratio Median Absolute Deviation'};

    % Array initialization
    ComputationTime = [];
    ClassificationSD = [];
    ClassificationMN = [];
    PredictionResults = [];

    % Perform statistical testing and classification for all selected subdivision levels
    for p = 1:length(numOfSubImagesVec)
        tic
        X = AllSubimageResults{p};

        % Removing corner subimages 
        % (not necessary for 4 subimages and more than 256 subimages, since most pixels would be background and removed in the next step of data preparation)
        if  numOfSubImagesVec(p)==64
            borderimages = [1 2 7 8 9 16 49 56 57 58 63 64];
            for k=1:length(borderimages)
                X(any(X(:,3)==borderimages(k), 2),:)=[];
            end
        elseif numOfSubImagesVec(p)==256
            borderimages = [[1 2 3 4 13 14 15 16 17 18 19 30 31 32 33 34 47 48 49 64] 257-[1 2 3 4 13 14 15 16 17 18 19 30 31 32 33 34 47 48 49 64]];
            for k=1:length(borderimages)
                X(any(X(:,3)==borderimages(k), 2),:)=[];
            end
        end


        % Remove subimages without bright pixels (NaNs)
        X(any(isnan(X)==1, 2),:) = [];

        % Separate a normal and a tumor core for future prediction
        NormalCoreForPred = 35;
        TumorCoreForPred = 13;
        coretumor = X(X(:,4) == TumorCoreForPred,:);
        corenormal = X(X(:,4) == NormalCoreForPred,:);
        X(any(X(:,4) == TumorCoreForPred, 2),:)=[];
        X(any(X(:,4) == NormalCoreForPred, 2),:)=[];

        % Select data between 1st and 99th percentiles to avoid extreme outliers
        PercentileValue = 99;
        Xnew = [];
        for j=1:4
            Xmat = X(X(:,2)==j,:);
            ubnormal = prctile(Xmat,PercentileValue);
            lbnormal = prctile(Xmat,100-PercentileValue);
            filtermat = Xmat;
            for i=1:size(Xmat,2)
                filtermat(:,i) = Xmat(:,i)>ubnormal(i) | Xmat(:,i)<lbnormal(i);
            end
            Xmat(any(filtermat==1, 2),:)=[];
            Xnew = [Xnew;Xmat];
        end

        % Reorder the columns and prepare for classification
        Yreordered = [Xnew(:,[1 2 3 4]) Xnew(:,[41 44 8 13 18 23 28]) Xnew(:,[51 54 10 15 20 25 30]) Xnew(:,[36 39 7 12 17 22 27]) Xnew(:,[31 34 6 11 16 21 26]) Xnew(:,[46 49 9 14 19 24 29]) Xnew(:,5)];
        Yclassification = Yreordered;

        coretumorReordered = [coretumor(:,[1 2 3 4]) coretumor(:,[41 44 8 13 18 23 28]) coretumor(:,[51 54 10 15 20 25 30]) coretumor(:,[36 39 7 12 17 22 27]) coretumor(:,[31 34 6 11 16 21 26]) coretumor(:,[46 49 9 14 19 24 29]) coretumor(:,5)];
        corenormalReordered = [corenormal(:,[1 2 3 4]) corenormal(:,[41 44 8 13 18 23 28]) corenormal(:,[51 54 10 15 20 25 30]) corenormal(:,[36 39 7 12 17 22 27]) corenormal(:,[31 34 6 11 16 21 26]) corenormal(:,[46 49 9 14 19 24 29]) corenormal(:,5)];
        coretumorClassification = coretumorReordered; 
        corenormalClassification = corenormalReordered; 

        ReorderedNameTags = [Xnametag([1 2 3 4]) Xnametag([41 44 8 13 18 23 28]) Xnametag([51 54 10 15 20 25 30]) Xnametag([36 39 7 12 17 22 27]) Xnametag([31 34 6 11 16 21 26]) Xnametag([46 49 9 14 19 24 29]) Xnametag(5)];
        ReorderedNameTagsClassification = ReorderedNameTags; 


        % % % % 4-group statistics
        % % % normal = Yreordered(Yreordered(:,1)==1,5:end); mad(mmm)
        % % % mmm = Yreordered(Yreordered(:,2)==2,5:end);
        % % % ppm = Yreordered(Yreordered(:,2)==3,5:end);
        % % % ppp = Yreordered(Yreordered(:,2)==4,5:end);
        % % % groupsstats = [mean(normal)' mad(normal)' mean(mmm)' mad(mmm)' mean(ppm)' mad(ppm)' mean(ppp)' mad(ppp)'];


        %%%%% Classification using data subsets (uncomment one at a time) %%%%%

        % Texture only
        % DiscardedCores = [2:6 12:13 18 19:20 26:27 33:34 40]; 

        % Pixel density
        % DiscardedCores = [2:size(Y,2)-1]; 

        % Intensity (Mean + MAD)
        % DiscardedCores = [2 3 4 7:size(Y,2)]; 

        % Polarimetric Paramters (Mean + MAD)
        % DiscardedCores = [2 3 4 7:11 14:18 21:25 26 28:32 33 35:40]; 

        % Intensity (Mean + MAD) + Pixel density
        % DiscardedCores = [2 3 4 7:size(Y,2)-1]; 

        % Polarimetric Paramters (Mean + MAD) + pixel density
        % DiscardedCores = [2 3 4 7:11 14:18 21:25 26 28:32 33 35:39]; 

        % Intensity (Mean + MAD + Texture)
        % DiscardedCores = [2 3 4 12:size(Y,2)]; 

        % Polarimetric Paramters (Mean + MAD + texture)
        % DiscardedCores = [2 3 4 18 26 33 40]; 

        % Intensity (Mean + MAD + Texture) + pixel density
        % DiscardedCores = [2 3 4 12:size(Y,2)-1];  

        % All parameters
        DiscardedCores = [2 3 4 18 26 33]; 

        % Discard unselected columns
        Yclassification(:,DiscardedCores)=[];
        coretumorClassification(:,DiscardedCores)=[];
        corenormalClassification(:,DiscardedCores)=[];
        ReorderedNameTagsClassification(DiscardedCores)=[];

        % Repeated cross-validation of logistic regression
        ClassMetrics = [];
        for ii = 1:10
            [TrainedClassifier, ClassifierPerformance] = LogisticRegCode(Yclassification);
            ClassMetrics = [ClassMetrics; ClassifierPerformance]; %performance for all repetitions
        end

        % mean and standard deviation of repeated classification
        ClassificationSD = [ClassificationSD; numOfSubImagesVec(p) std(ClassMetrics)];
        ClassificationMN = [ClassificationMN; numOfSubImagesVec(p) mean(ClassMetrics)];
        ComputationTime = [ComputationTime; numOfSubImagesVec(p) toc];

        %%%%%%%%%%%%%%% Prediction of Normal and Tumor Cores %%%%%%%%%%%%%%

        FitModel=TrainedClassifier.GeneralizedLinearModel;

        [probtumor,tumorci] = predict(FitModel,coretumorClassification(:,2:end));
        [probnormal,normalci] = predict(FitModel,corenormalClassification(:,2:end));

        PredictedCoresProbs = [probnormal ; probtumor];
        PredictedCoresBinary = PredictedCoresProbs;
        PredictedCoresBinary(PredictedCoresBinary>0.5)=2;
        PredictedCoresBinary(PredictedCoresBinary<=0.5)=1;
        PredictedCoresResp = [1*ones(length(probnormal),1) ; 2*ones(length(probtumor),1)];


        PostProbMatrix = PredictedCoresProbs;

        PostProbMatrix(PredictedCoresResp==2)=PostProbMatrix(PredictedCoresResp==2)-1;

        % Brier Score
        BrierScore = sum((PostProbMatrix).^2)/length(PostProbMatrix);

        % AUROC
        [X4,Y4,T4,AUC] = perfcurve(PredictedCoresResp,PredictedCoresProbs,2);

        % Confusion Matrix Metrics
        TP = sum(PredictedCoresBinary==2 & PredictedCoresResp==2);
        TN = sum(PredictedCoresBinary==1 & PredictedCoresResp==1);
        FP = sum(PredictedCoresBinary==2 & PredictedCoresResp==1);
        FN = sum(PredictedCoresBinary==1 & PredictedCoresResp==2);

        TPR = TP/(TP+FN);
        TNR = TN/(TN+FP);
        FPR = FP/(FP+TN);
        FNR = FN/(FN+TP);

        PPV = TP/(TP+FP);
        NPV = TN/(TN+FN);
        FDR = FP/(FP+TP);
        FOR = FN/(FN+TN);

        F1score = 2*PPV*TPR/(PPV+TPR);
        Accuracy = (TP+TN)/(TP+TN+FP+FN);

        ConfuseMat = [TN TP FN FP TNR TPR FNR FPR NPV PPV FOR FDR Accuracy F1score BrierScore AUC];
        PredictionResults = [PredictionResults; numOfSubImagesVec(p) ConfuseMat];
    end
    RepeatedCrossValidationMean = array2table(ClassificationMN,'VariableNames',{'Subimages','TN', 'TP', 'FN', 'FP', 'TNR', 'TPR', 'FNR', 'FPR', 'NPV', 'PPV', 'FOR', 'FDR', 'Accuracy', 'F1-score', 'BrierScore', 'AUROC'});
    RepeatedCrossValidationStdev = array2table(ClassificationSD,'VariableNames',{'Subimages','TN', 'TP', 'FN', 'FP', 'TNR', 'TPR', 'FNR', 'FPR', 'NPV', 'PPV', 'FOR', 'FDR', 'Accuracy', 'F1-score', 'BrierScore', 'AUROC'});
    PredictionTable = array2table(PredictionResults,'VariableNames',{'Subimages','TN', 'TP', 'FN', 'FP', 'TNR', 'TPR', 'FNR', 'FPR', 'NPV', 'PPV', 'FOR', 'FDR', 'Accuracy', 'F1-score', 'BrierScore', 'AUROC'});
    clearvars -except ImageTag GroupTag CDcombined LDcombined Rcombined DCPcombined CIntcombined numOfSubImagesVec AllParamMatrices AllSubimageResults SubimageRuntime RepeatedCrossValidationMean RepeatedCrossValidationStdev PredictionTable MultipleComparisons Yreordered Yclassification ReorderedNameTags ReorderedNameTagsClassification Boxplots Xnew
    fprintf('Classification Complete!\n ----------------------------- \n')
end

%% %%%%%%%%%%%%%%%%%% Kruskal-Wallis + Dunn-Bonferroni %%%%%%%%%%%%%%%%%%%%

if MultipleComparisons==1

    ChisqVec = [];
    Pvals = [];

    for i = 5:size(Yreordered,2)
        [p,tbl,stats] = kruskalwallis(Yreordered(:,i),Yreordered(:,2),'off');
        ChisqVec = [ChisqVec cell2mat(tbl(2,5))];
        multres = multcompare(stats, 'CType', 'bonferroni','Display','off');
        Pvals = vertcat(Pvals,multres(:,6)');
    end

    MultipleComps = [table(char(ReorderedNameTags(5:end)'),'VariableName',{'Parameter Name'}) array2table(Pvals,'VariableNames',{'Normal & -/-/-','Normal & +/+/-','Normal & +/+/+','-/-/- & +/+/-','-/-/- & +/+/+','+/+/- & +/+/+'})];
    clearvars -except ImageTag GroupTag CDcombined LDcombined Rcombined DCPcombined CIntcombined numOfSubImagesVec AllParamMatrices AllSubimageResults SubimageRuntime RepeatedCrossValidationMean RepeatedCrossValidationStdev PredictionTable Yreordered Yclassification ReorderedNameTags ReorderedNameTagsClassification Boxplots Xnew MultipleComps
    fprintf('Statistical Analysis Complete!\n ----------------------------- \n')
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%% Batch Boxplots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Boxplots==1
    
    XnametagWithLog = {'2 Groups','4 Groups','SubimageTag','Core Number','Log_{10}(Pixel Density)',...
        'SHG-CD Contrast','DCP Contrast','SHG Intensity Contrast','SHG-LD Contrast','R-Ratio Contrast',...
        'SHG-CD Correlation','DCP Correlation','SHG Intensity Correlation','SHG-LD Correlation','R-Ratio Correlation',...
        'SHG-CD Entropy','DCP Entropy','SHG Intensity Entropy','SHG-LD Entropy','R-Ratio Entropy',...
        'Log_{10}(SHG-CD ASM)','Log_{10}(DCP ASM)','Log_{10}(SHG Intensity ASM)','Log_{10}(SHG-LD ASM)','Log_{10}(R-Ratio ASM)',...
        'SHG-CD IDM','DCP IDM','SHG Intensity IDM','SHG-LD IDM','R-Ratio IDM',...
        'SHG-CD Mean','SHG-CD Median','SHG-CD Standard Deviation','SHG-CD MAD','SHG-CD Median Absolute Deviation',...
        'DCP Mean','DCP Median','DCP Standard Deviation','DCP MAD','DCP Median Absolute Deviation',...
        'Log_{10}(SHG Intensity Mean)','Log_{10}(SHG Intensity Median)','Log_{10}(SHG Intensity Standard Deviation)','Log_{10}(SHG Intensity MAD)','Log_{10}(SHG Intensity Median Absolute Deviation)',...
        'SHG-LD Mean','SHG-LD Median','SHG-LD Standard Deviation','SHG-LD MAD','SHG-LD Median Absolute Deviation',...
        'R-Ratio Mean','R-Ratio Median','R-Ratio Standard Deviation','R-Ratio MAD','R-Ratio Median Absolute Deviation'};

    GroupLabels = {'Normal','-/-/-','+/+/-','+/+/+'};

    colors1 = [237,248,251;179,205,227;65, 137, 196;71, 65, 135]/255;
    colors2 = [255, 211, 181;245, 136, 73;224, 20, 20;156, 0, 41]/255;
    colors3 = [253, 255, 143;169, 214, 90;81, 171, 36;3, 102, 0]/255;

    for i = 5:size(Xnew,2)
        if i==5
            colors = colors2; 
        elseif i>5 && i<31
            colors = colors3;
        else 
            colors = colors1;
        end
        plotvar = Xnew(:,i);
        grouping = Xnew(:,2);

        if i==5 || i==21 || i==22 || i==23 || i==24 || i==25 || i==41 || i==42 || i==43 || i==44 || i==45
            plotvar = log10(plotvar); 
        end

        figure;
        h2 = boxplot(plotvar,grouping,'Widths',0.8,'Notch','on','Symbol','o','OutlierSize',4,'Colors','k','Labels',GroupLabels);
        h = findobj(gca,'Tag','Box');
        for j=1:length(h)
            patch(get(h(j),'XData'),get(h(j),'YData'),colors(end+1-j,:),'FaceAlpha',1);
        end
        hold on 
        h2 = boxplot(plotvar,grouping,'Widths',0.8,'Notch','on','Symbol','o','OutlierSize',4,'Colors','k','Labels',GroupLabels);
        haspect = 1;
        vaspect = 1;
        pbaspect([haspect vaspect 1])
        ylabel(XnametagWithLog{i})
        set(gcf, 'Color', 'w');
        if haspect==2 && vaspect==1
            set(gca,'FontSize',10);
        else
            set(gca,'FontSize',20);
        end
        ax = gca;
        tickvalues = ax.YTick; tick = diff(tickvalues);
        ylim([min(plotvar)-0.15*tick(end) max(plotvar)+2*tick(end)])
        box off

    end
end

%%

clearvars -except ImageTag GroupTag CDcombined LDcombined Rcombined DCPcombined CIntcombined AllSubimageResults RepeatedCrossValidationMean RepeatedCrossValidationStdev PredictionTable Yreordered Yclassification ReorderedNameTags ReorderedNameTagsClassification MultipleComps


% Turn warnings back on
warning('on','all');
warning('query','all');

% Print Results
fprintf('\n\n\n ****************************** Multiple Comparisons with Kruskal Wallis Test + Bonferroni Correction ****************************** \n\n\n')
disp(MultipleComps)
fprintf('\n\n ***************************************** Repeated Cross-validated Classification Results  ****************************************** \n\n')
fprintf('\n\n\n ******************************************** Classification Performance Metrics Mean  ********************************************* \n\n\n')
disp(RepeatedCrossValidationMean)
fprintf('\n\n\n ************************************* Classification Performance Metrics Standard Deviation  ************************************** \n\n\n')
disp(RepeatedCrossValidationStdev)
fprintf('\n\n\n ************************************************ Prediction Performance Metrics  ************************************************** \n\n\n')
disp(PredictionTable)