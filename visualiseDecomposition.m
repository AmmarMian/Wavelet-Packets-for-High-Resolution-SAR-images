%% Description
% To visualise Decmposition (f, theta)
%% Specifiactions
% Other m-files required: target creation functions
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
%     clear variables
    clc
    close all
    set(0,'defaulttextinterpreter','latex')
    
%% Reading Image
    addpath('Hyperimages', 'Target Dictionary', 'Detectors')
    imagePath = '../../../../Data/SANDIA/MiniSAR20050519p0010image003.mat';
    [Image, Header] = readSandia(imagePath, true, [500, 1, 1000, 900]);
%     [Image, Header] = readSDMS('..\..\..\..\Data\SDMS', 'FP0121', true, [2200, 1000, 2701, 2001]);
%     Image = Image(:,:,2);
%          
%% Decomposition
    Decomposition = struct;
    Decomposition.method = 'Exponential';
    Decomposition.J = 1;
    Decomposition.numberOfBands = 1;
    Decomposition.numberOfLooks = 6;
    Decomposition.d1 = inf;
    Decomposition.d2 = inf;
    Decomposition.visualise = true;
    Decomposition.dyn_dB = 50;
    [ImageFTheta, k, theta] = decomposeFTheta(Image, Header, Decomposition);

    Decomposition.d1 = 3;
    Decomposition.d2 = 3;
    [ImageFTheta, k, theta] = decomposeFTheta(Image, Header, Decomposition);

%% RGB visualisation for numberofbands = 3
    Decomposition.visualise = false;
    Decomposition.numberOfBands = 3;
    Decomposition.numberOfLooks = 1;
    [ImageFTheta, k, theta] = decomposeFTheta(Image, Header, Decomposition);
    figure
    toShow = abs(ImageFTheta(:,1:3:end,:));
    toShow = toShow/max(toShow(:)) * 255;
    imagesc(Header.x, Header.y, toShow)
    xlabel('$x$', 'interpreter', 'latex')
    ylabel('$y$', 'interpreter', 'latex')
    title('Spectral behaviour', 'interpreter', 'latex')
    set(gca, 'ydir', 'normal')
    axis image
    
%% RGB visualisation for numberofLooks = 3
    Decomposition.visualise = false;
    Decomposition.numberOfBands = 1;
    Decomposition.numberOfLooks = 3;
    [ImageFTheta, k, theta] = decomposeFTheta(Image, Header, Decomposition);
    figure
    toShow = abs(ImageFTheta(1:3:end,:,:));
    toShow = toShow/max(toShow(:)) * 255;
    imagesc(Header.x, Header.y, toShow)
    xlabel('$x$', 'interpreter', 'latex')
    ylabel('$y$', 'interpreter', 'latex')
    title('Angular behaviour', 'interpreter', 'latex')
    set(gca, 'ydir', 'normal')
    axis image
%% ------------- END CODE ----------------
