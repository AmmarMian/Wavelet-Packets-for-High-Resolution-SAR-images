%% Description
% Read Sandia images
%% Specifiactions
% Other m-files required: none
% MAT-files required: none
%% Authors
% Authors: Ammar Mian, Ph.D., Signal processing at Supelec SONDRA
% Email address: ammar.mian@centralesupelec.fr
%% Copyright
% Copyright 2017 CentraleSupelec
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%     http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
%% ------------- BEGIN CODE --------------
function [Image, Header] = readSDMS(dataPath, Pass, resizeImage, cropIndexes)
% readSandia - Read and crop SANDIA Images
% Syntax:  [Image, Header] = readSandia(imagePath, resizeImage, cropIndexes)
%          WARNING: X = Range and Y = Azimuth
% Inputs:
%    - dataPath = path to folder containing SDMS Data. Ex:
%    'C:\Users\mian_amm\Dropbox\Thèse\Data\SDMS'.
%    - Pass = Number of the pass to read. Ex: 'FP120'.
%    - resizeImage = If we want to crop images
%    - ceopIndexes = Indexes to crop the image if resizeImage is true :
%     [XLeftUpCorner, YleftUpCorner, XLeftBottomCorner, YRightUpCorner]
%
% Outputs:
%   - Image = image array
%   - Header = image info

    % Init variable name for SDMS
    FolderCoPol = 'CCD-CP-CoPol'; % For HH and VV Polarisation
    FolderXPol = 'CCD-CP-XPol'; % For HV Polarisation
    FileName = 'c00007a283p50.mat';
    Image = zeros(4501,4501,3); % Init polar Image
    Header = struct;
    Header.RgCnt = 4501;
    Header.AzCnt = 4501;
    Header.RgPixelSz = 0.2;
    Header.AzPixelSz = 0.2;
    Header.numberOfPol = 3; 
    Header.NominalCenterFreq = 9.6e+9;
    imagePathHH = sprintf('%s\\%s\\HH\\%s\\%s', dataPath, FolderCoPol, Pass, FileName);
    imagePathVV = sprintf('%s\\%s\\VV\\%s\\%s', dataPath, FolderCoPol, Pass, FileName);
    imagePathHV = sprintf('%s\\%s\\HV\\%s\\%s', dataPath, FolderXPol, Pass, FileName);
    
    % Reading files
    data = load(imagePathHH);
    Image(:,:,1) = data.SARdataOut;
    data = load(imagePathVV);
    Image(:,:,2) = data.SARdataOut;
    data = load(imagePathHV);
    Image(:,:,3) = data.SARdataOut;
    Image = permute(Image, [2,1,3]);
    if resizeImage
        Image = Image(cropIndexes(1):cropIndexes(3), cropIndexes(2):cropIndexes(4),:); % Cropping the Image
        Header.AzCnt = cropIndexes(3) - cropIndexes(1) + 1; % Updating the Azimuth count
        Header.RgCnt = cropIndexes(4) - cropIndexes(2) + 1; % Updating the Range count
    end
    
    % Spatial vectors
    x = linspace(-Header.RgCnt * Header.RgPixelSz / 2, Header.RgCnt * Header.RgPixelSz / 2 - Header.RgPixelSz, Header.RgCnt); % Range
    y = linspace(-Header.AzCnt * Header.AzPixelSz / 2, Header.AzCnt * Header.AzPixelSz / 2 - Header.AzPixelSz, Header.AzCnt); % Azimuth
    c = 3e8; % Speed of light
    kudopcentral = 0;
    kcentral = 2 * Header.NominalCenterFreq / c;
    kx = linspace(kcentral - 1 / (2 * Header.RgPixelSz), kcentral + 1 / (2 * Header.RgPixelSz) - 1 / (Header.RgPixelSz * Header.RgCnt), Header.RgCnt);
    ky = linspace(kudopcentral - 1 / (2 * Header.AzPixelSz), kudopcentral + 1 / (2 * Header.AzPixelSz) - 1 / (Header.AzPixelSz * Header.AzCnt), Header.AzCnt);
    Header.x = x;
    Header.y = y;
    Header.kx = kx;
    Header.ky = ky;
    Header.Nx = length(x);
    Header.Ny = length(y);
end


%% ------------- END CODE ----------------
