function [CCM_sizeVec,dEVec,T_Vec,combs] = runOnGaps()
addpath(genpath('C:\Users\user\Desktop\Eli_SVN\Unispectral Program'))
imageDir = 'C:\Users\user\Google Drive\SIS Shared Files\Experiment results and measurements\LM Piezo Filter\LM_Piezo2_80nm results\Images\20-1-16\';
imageDir = [imageDir,'25- LED illum'];
% imageDir = [imageDir,'22- Piezo filter_closeUp'];
chart = load('MacbethColorCheckerData.mat');
colors = chart.colors;
chart = chart.chart;

% ColorCheckerLabOptimal = struct2array(   load (  'ColorCheckerLabOptimal.mat'));
ColorCheckerLabOptimal = struct2array(   load (  'CCLabLEDlight.mat'));

% colors = {'darkskin','lightskin','bluesky','foliage','blueflower','bluishgreen';'orange','purplishblue','moderatered','purple','yellowgreen','orangeyellow';'blue','green','red','yellow','magenta','cyan';'white','neutral8','neutral65','neutral5','neutral35','black'};

files = dir([imageDir,'\*.bmp']);

for ii=1:length(files)
    fileName = files(ii).name;
    temp = imread([imageDir,'\',fileName]);
    RGBW_temp(:,:,ii) = temp(:,:,1);
end
% imshow(RGBW_temp(:,:,1:3));
RGBW = RGBW_temp(:,:,[1,2,3]);
% RGBW = RGBW_temp;
% RGBW = imcrop(RGBW);
% close all


Iwhitebalanced = im2double(  RGBW_temp);
% Iwhitebalanced = DynamicRange(AvrResize(RGBW2,20),1);
% Iwhitebalanced = DynamicRange(RGBW2,1);
I = Iwhitebalanced;     %   Normalized image

% hist(reshape(I,[],3),100);
% colormap([1 0 0; 0 1 0; 0 0 1]);
%%
s = size(RGBW);
figure;imshow(RGBW)
title('Linear image','fontsize',20)
title('Click the center of the top left patch and press [Enter]','fontsize',20)
[x1, y1]    = ginput;
title('Click the center of the bottom right patch and press [Enter]','fontsize',20)
[x2, y2]    = ginput;
close all;

masks = makeChartMask(RGBW,chart,colors, 1, x1, y1, x2, y2);
%%
combs = combinator(size(RGBW_temp,3),3,'c');
combs = mat2cell(combs,ones(length(combs),1),3);
[CCM_sizeVec,dEVec,T_Vec] = cellfun(@calculateCCM,combs,'UniformOutput', false);

    function [CCM_size,dE,T] = calculateCCM(comb)
        IMAGE = Iwhitebalanced(:,:,comb);
        RGB = getChartRGBvalues(IMAGE,masks,colors);
        
        C   = makecform('lab2xyz');
        XYZ = applycform (ColorCheckerLabOptimal, C);   
        T   = XYZ' * pinv(RGB)';  %  
        CCM_size = (sqrt(sum(T(:).^2)));
%         disp(['CCM size: ',num2str(sqrt(sum(T(:).^2)))])
        % save ('CC Matrix.mat', 'T');
        
        CorrectedCCXYZ = T*RGB';
        LabSpace = rgb2lab (CorrectedCCXYZ');
        dE = mean(sqrt( (LabSpace(:,2)-ColorCheckerLabOptimal(:,2)).^2 ...
            +   (LabSpace(:,3)-ColorCheckerLabOptimal(:,3)).^2 ...
            +   (LabSpace(:,1)-ColorCheckerLabOptimal(:,1)).^2));
%         disp(['dE: ',num2str(dE)])
        
        if false
            T_temp = T;
            T(1,:) = T(1,:)/sum(T(1,:));
            T(2,:) = T(2,:)/sum(T(2,:));
            T(3,:) = T(3,:)/sum(T(3,:));
            
        end
        % APPLY T TO A NOVEL IMAGE
        
%         Ixyz = reshape((T*reshape(Iwhitebalanced,[s(1)*s(2) 3])')',[s(1) s(2) 3]);
%         
%         Irgb = XYZ2ProPhoto(Ixyz); % ProPhoto is a wide gamut RGB space that won't clip most colors.
%         % Irgb = xyz2rgb(Ixyz); % ProPhoto is a wide gamut RGB space that won't clip most colors.
%         
%         % clc
%         % WP      = whitepoint('D65');
%         % C       = makecform('xyz2srgb');   % Color space transformation for the Lab plot lables RGB colors
%         % Irgb1   = applycform (uint16(Ixyz.*2^8), C);
%         % Irgb    = DynamicRangeGPU(Irgb,1);
%         % Irgb1   = DynamicRange(Irgb1,1);
%         % Irgb = wbalanceGPU(Irgb);
%         ShowTwoImages(Irgb, Iwhitebalanced, 'CC corrected', 'Original');
        % imshow(Irgb);title('Color transformed image','fontsize',20)
%         %%
%         [FileName,PathName] = uiputfile([imageDir,'\*.bmp'],'Save Image As');
%         imwrite(Irgb,[PathName,FileName]);
%         imwrite(Iwhitebalanced,[PathName,'Original_Cropped.bmp']);
%         
%         return
        % box = 2;
        % Irgb_deonised(:,:,1)= medfilt2(Irgb(:,:,1), [box,box]);
        %     Irgb_deonised(:,:,2)= medfilt2(Irgb(:,:,2), [box,box]);
        %     Irgb_deonised(:,:,3)= medfilt2(Irgb(:,:,3), [box,box]);
        %     ShowTwoImages(Irgb, Irgb_deonised, 'CC corrected', 'Denoiser');
    end
end