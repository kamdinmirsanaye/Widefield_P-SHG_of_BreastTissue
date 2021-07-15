function feat = GLCMFeat(imin,imax,Ng,image,D)
%SUMMARY: Calculates features of an image
%   This function calculates the glcms of the image passed in, then
%   computes and returns the features contrast, correlation, and entropy.
%
%   Params:
%   imin (double)           minimum value of the image
%   imax (double)           max value of the image
%   image (double)          background + noise corrected image
%   Ng (double)             Number of grey levels to use (should use 128)
%   D (double)              distance between neighbours
%   Returns:
%   feat (double)           array of contrast, correlation, entropy, ASM, IDM

% directions of analysis
if D == 1
    offsets = [0 1; -1 1; -1 0; -1 -1];
elseif D == 2
    offsets = [0 2; -1 2; -2 2; -2 1; -2 0; -2 -1; -2 -2; -1 -2];
elseif D == 3
    offsets = [0 3; -1 3; -2 3; -3 3; -3 2; -3 1; -3 0; -3 -1; -3 -2; ...
        -3 -3; -2 -3; -1 -3];
end

% Calculate glcms
totalGLCMs = graycomatrix(image,'Offset',offsets,'GrayLimits',...
    [imin imax],'NumLevels',Ng,'Symmetric',true);  % calculate GLCM matrices

% Normalize glcms for entropy and IDM calculation
for p = 1:size(totalGLCMs,3)
    normGLCMs(:,:,p) = Normalize(totalGLCMs(:,:,p));
end

% Calculate all features
feat = [];
stats = graycoprops(totalGLCMs);
feat(1) = mean(stats.Contrast);
feat(2) = mean(stats.Correlation);

entropy = -sum(sum( normGLCMs .* log2(normGLCMs + eps))); % Calculate entropy
feat(3) = mean(entropy);

feat(4) = mean(stats.Energy);   % ASM is also known as energy
feat(5) = mean(IDM(normGLCMs));

end

function glcm = Normalize(glcm)
%SUMMARY: This function normalizes any GLCMs passed into it

if any(glcm(:))
  glcm = glcm ./ sum(glcm(:));
end
% matSum = sum(sum(GLCMs));
% if matSum == 0
%     normGLCMs = 0;
% elseif matSum ~= 0
%     normGLCMs = GLCMs ./ matSum;
% end
end


function IDMarray = IDM(glcm)
%SUMMARY: This function calculates the inverse difference moment
%
%   Params:
%   glcm        NORMALIZED glcms

IDMarray = zeros(1,size(glcm,3));

for k = 1:size(glcm,3)
    for i = 1:size(glcm,1)
        for j = 1:size(glcm,2)
            IDMarray(k) = IDMarray(k) + (glcm(i,j,k)/( 1 + (i - j)^2));
        end
    end
end

end
