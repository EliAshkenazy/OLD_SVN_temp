%TO DO:
%Check results with optimized CCM versus the old results (which were
%calculated with simple inverse and also were done on 24 patches.)
%
%Add CCM which is optimized by dE and CCM size (first find the best dE then
%define an accptable range of dE and look for minimum of CCM size)


addpath(genpath('C:\Users\user\Desktop\Eli_SVN\Unispectral Program'))
imageDir = 'C:\Users\user\Google Drive\SIS Shared Files\Experiment results and measurements\LM Piezo Filter\LM_Piezo2_80nm results\Images\20-1-16\';
imageDir = [imageDir,'25- LED illum']
% imageDir = [imageDir,'20 - Color filters_closeUp'];

% imageDir = [imageDir,'22- Piezo filter_closeUp'];

% files = dir(imageDir);
files = dir([imageDir,'\*.bmp']);
for ii=1:3
    index = [6,8,9];
%     index = [1,2,3];

    fileName = files(index(ii)).name
    temp = imread([imageDir,'\',fileName]);
    RGBW_temp(:,:,ii) = temp(:,:,1);
end
imshow(RGBW_temp);
RGBW = RGBW_temp(:,:,[2,1,3]);
% RGBW = RGBW_temp;
RGBW = imcrop(RGBW);
close all
%%
% filename = uigetfile('*.tif');
% I = im2double(imread(filename));
% load('image for CC.mat');
% Iwhitebalanced = DynamicRangeGPU(RGBW,1);

% box = 2; clear RGBW2
% RGBW2(:,:,1)= medfilt2(RGBW(:,:,1), [box,box]);
%     RGBW2(:,:,2)= medfilt2(RGBW(:,:,2), [box,box]);
%     RGBW2(:,:,3)= medfilt2(RGBW(:,:,3), [box,box]);
%     ShowTwoImages(RGBW, RGBW2, 'CC corrected', 'Denoiser'); 
RGBW2=RGBW;

Iwhitebalanced = im2double(  RGBW2);
% Iwhitebalanced = DynamicRange(AvrResize(RGBW2,20),1);
% Iwhitebalanced = DynamicRange(RGBW2,1);

I = Iwhitebalanced;     %   Normalized image

% hist(reshape(I,[],3),100);
% colormap([1 0 0; 0 1 0; 0 0 1]);
%%
s = size(I);    
figure;imshow(I)
title('Linear image','fontsize',20)
title('Click the center of the top left patch and press [Enter]','fontsize',20)
[x1, y1]    = ginput;
title('Click the center of the bottom right patch and press [Enter]','fontsize',20)
[x2, y2]    = ginput;
close all;

% LOAD COLOR CHART DATA
% This example is given for a Macbeth ColorChecker - you can load other charts data saved in the
% same format. To see the format, open MacbethColorChecker.m.

load MacbethColorCheckerData.mat

%   Loading the measured
%%

% load (  'ColorCheckerLabOptimal.mat');
load (  'CCLabLEDlight.mat');

% MAKE MASKS FOR THE PATCHES OF THE COLOR CHART
% This script works for a Macbeth chart but it can be modified to work for other charts.
% It places the masks for each patch on the image, and waits for the user to drag each mask over the correct patch.
% Once each mask is aligned its corresponding patch, the user should double click the first patch in the chart
% -- that will accept the masks for all patches. In a Macbeth ColorChecker, the first patch is the Dark Skin.
%%
masks = makeChartMask(I,chart,colors, 1, x1, y1, x2, y2);

% APPLY WHITE BALANCE
% The locatio of the WS and DS are for a Macbeth chart, modify as needed.
% % Iwhitebalanced = I;
% Ilinearized = I;
% whitePatch = [4 1];
% darkPatch = [4 6];
% Iwhitebalanced = whiteBalance(Ilinearized,masks.(colors{whitePatch(1),whitePatch(2)}).mask,masks.(colors{darkPatch(1),darkPatch(2)}).mask,0.8);
% figure;imshow(Iwhitebalanced)
% title('white balanced image','fontsize',20)

%% DERIVE THE TRANSFORMATION MATRIX T


RGB = getChartRGBvalues(Iwhitebalanced,masks,colors);

% C   = makecform('lab2xyz');
% XYZ = applycform (ColorCheckerLabOptimal, C); 
% 
% % XYZ = getChartXYZvalues(chart,colors);   %    Original CC values in XYZ
% T   = XYZ' * pinv(RGB)'  %   Color correction Matrix (3x3)
% disp(['CCM size: ',num2str(sqrt(sum(T(:).^2)))])
% % save ('CC Matrix.mat', 'T');
% 
% CorrectedCCXYZ = T*RGB';
% LabSpace = rgb2lab (CorrectedCCXYZ');
% dE = mean(sqrt( (LabSpace(:,2)-ColorCheckerLabOptimal(:,2)).^2 ...
%                         +   (LabSpace(:,3)-ColorCheckerLabOptimal(:,3)).^2 ...
%                         +   (LabSpace(:,1)+28-ColorCheckerLabOptimal(:,1)).^2));
% disp(['dE: ',num2str(dE)])
% 
% if false
%     T_temp = T;
%     T(1,:) = T(1,:)/sum(T(1,:));
%     T(2,:) = T(2,:)/sum(T(2,:));
%     T(3,:) = T(3,:)/sum(T(3,:));
%     
% end
% APPLY T TO A NOVEL IMAGE

[CCM, XYZ] = OptimizedColorCorrection( ColorCheckerLabOptimal,RGB );
T = CCM;

Ixyz = reshape((T*reshape(Iwhitebalanced,[s(1)*s(2) 3])')',[s(1) s(2) 3]);

Irgb = XYZ2ProPhoto(Ixyz); % ProPhoto is a wide gamut RGB space that won't clip most colors.
% Irgb = xyz2rgb(Ixyz); % ProPhoto is a wide gamut RGB space that won't clip most colors.

% clc
% WP      = whitepoint('D65');
% C       = makecform('xyz2srgb');   % Color space transformation for the Lab plot lables RGB colors
% Irgb1   = applycform (uint16(Ixyz.*2^8), C);
% Irgb    = DynamicRangeGPU(Irgb,1);
% Irgb1   = DynamicRange(Irgb1,1);
% Irgb = wbalanceGPU(Irgb);
ShowTwoImages(Irgb, Iwhitebalanced, 'CC corrected', 'Original'); 
% imshow(Irgb);title('Color transformed image','fontsize',20)
%%
[FileName,PathName] = uiputfile([imageDir,'\*.bmp'],'Save Image As');
if FileName~=0
    imwrite(Irgb,[PathName,FileName]);
    imwrite(Iwhitebalanced,[PathName,'Original_Cropped.bmp']);
end
return
% box = 2;
% Irgb_deonised(:,:,1)= medfilt2(Irgb(:,:,1), [box,box]);
%     Irgb_deonised(:,:,2)= medfilt2(Irgb(:,:,2), [box,box]);
%     Irgb_deonised(:,:,3)= medfilt2(Irgb(:,:,3), [box,box]);
%     ShowTwoImages(Irgb, Irgb_deonised, 'CC corrected', 'Denoiser'); 
