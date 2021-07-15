function blocks = DivideImage(image,numOfSubimages)
% SUMMARY: Divides an image into smaller sections (blocks)
%   This function divides an image into subimages, where number of
%   images is determined by the user. The function assumes the image is
%   evenly divisible into the number of subimages, and all subimages are
%   square (this can be modified)
% 
%   Params:
%   image                   background + noise corrected image to be divided
%   numOfSubimages          number of subimages to divide image into
%   Returns:
%   blocks (cell)           cell array containing all image sections

% Size of image
imageRows       = size(image,1);
imageColumns    = size(image,2);
imageSize       = imageRows * imageColumns;

% Determine size of *square* subimages
blockSize       = imageSize / numOfSubimages;
blockRows       = sqrt(blockSize);      % for a square block
blockCols       = blockRows;            % for a square block

numOfRows = floor(imageRows / blockRows);
numOfCols = floor(imageColumns / blockCols);
blockVectorR = blockRows * ones(1, numOfRows);
blockVectorC = blockCols * ones(1, numOfCols);


% Divide image into blocks
% cell array "blocks" contains a subimage in each cell
blocks = mat2cell(image, blockVectorR, blockVectorC);

end