function [CCM, XYZ] = OptimizedColorCorrection( ccLabOptimal,SensorOutput )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% 
% 1. OPtimal dE
% 2. Optimal CCMsize
% 3. Sum of elemsnt=1
% 4. only 18 patches
SubOptimized = 0;
ccmColorSpace = 'xyz';

% OptimalCCLab = ccLabOptimal(1:18,:);
OptimalCCLab = ccLabOptimal(1:24,:);

% Get Color Checker optimal XYZ values
ColorCheckerXYZOptimal = lab2xyz (ccLabOptimal);
switch ccmColorSpace
    case 'xyz'
%         ColorCheckerXYZOptimal = lab2xyz (ColorCheckerLabOptimal);
    case 'srgb'
        ColorCheckerXYZOptimal = lab2rgb (ccLabOptimal);
    case 'lrgb'
        ColorCheckerXYZOptimal = lab2xyz (ccLabOptimal);
        ColorCheckerXYZOptimal = my_xyz2rgb(ColorCheckerXYZOptimal);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% OptimalCCLab = ColorCheckerLabOptimal(1:18,:);

%XYZ
CCMinit =  SensorOutput\ColorCheckerXYZOptimal;

if false
    CCMinit
    CorrectedCCLab1 = SensorOutput*CCMinit;
    CorrectedCCLab1 = CorrectedCCLab1(:,1:18)';
    
    %XYZ
    switch ccmColorSpace
        case 'xyz'
            CorrectedCCLab1 = xyz2lab (CorrectedCCLab1);
        case 'srgb'
            CorrectedCCLab1 = rgb2lab (CorrectedCCLab1);
        case 'lrgb'
            CorrectedCCLab1 = my_rgb2xyz(CorrectedCCLab1);
            CorrectedCCLab1 = xyz2lab (CorrectedCCLab1);
    end
    
    CorrectedCCLab1 = xyz2lab (CorrectedCCLab1);
    
    d1 = OptimalCCLab - CorrectedCCLab1;
    d1 = d1.^2;
    s1 = sum(d1');
    s1 = sqrt(s1);
    delta1 = mean(s1)
    size1 = sqrt(sum(CCMinit(:).*CCMinit(:)))
end

%==========================================================================
%   Optimization by minimal dE

%   sets optimization parameters
opt = optimset('Display','off','TolFun',0.05 );

if (SubOptimized)
    CCMVec = CCMinit;
    delta = LabError(CCMinit);
else
    [CCMVec, delta,ext,msg] = fminsearch(@LabError, CCMinit,opt);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function delta = LabError(CCMatrix)
        
        % XYZ
        CorrectedCCXYZ = SensorOutput*CCMatrix;
        
        switch ccmColorSpace
            case 'xyz'
                CorrectedCCLab = xyz2lab (CorrectedCCXYZ);
            case 'srgb'
                CorrectedCCLab = rgb2lab (CorrectedCCXYZ);
            case 'lrgb'
                CorrectedCCXYZ = my_rgb2xyz(CorrectedCCXYZ);
                CorrectedCCLab = xyz2lab (CorrectedCCXYZ);
        end                
        
%         CorrectedCCLab = CorrectedCCLab(1:18,:);
                CorrectedCCLab = CorrectedCCLab(1:24,:);

        d = OptimalCCLab - CorrectedCCLab;
        d = d.^2;
        s = sum(d');
        s = sqrt(s);
        delta = mean(s);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%==========================================================================
CCM = reshape(CCMVec, 3,3);
XYZ = SensorOutput*CCM;

if true
    delta
    ext
    msg
    sizeOpt = sqrt(sum(CCM(:).*CCM(:)))
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % XYZ = getChartXYZvalues(chart,colors);   %    Original CC values in XYZ
% T   = XYZ' * pinv(ccMeasured)'  %   Color correction Matrix (3x3)
% disp(['CCM size: ',num2str(sqrt(sum(T(:).^2)))])
% % save ('CC Matrix.mat', 'T');
% 
% CorrectedCCXYZ = T*ccMeasured';
% LabSpace = rgb2lab (CorrectedCCXYZ');
% dE = mean(sqrt( (LabSpace(:,2)-ccLabOptimal(:,2)).^2 ...
%                         +   (LabSpace(:,3)-ccLabOptimal(:,3)).^2 ...
%                         +   (LabSpace(:,1)+28-ccLabOptimal(:,1)).^2));
%                     
% end
end
