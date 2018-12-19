%% Description
% Simulation with a target added through a sterring vector and test of
% detection
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
    clear variables
    clc
    close all
    set(0,'defaulttextinterpreter','latex')
    
%% Reading Image
    addpath('..', '../Hyperimages', '../Target Dictionary', '../Detectors')
    imagePath = '../../../../../Data/SANDIA/MiniSAR20050519p0010image002.mat';
    [Image, Header] = readSandia(imagePath, true, [1, 1, 301, 801]);
    %     [Image, Header] = readSDMS('C:\Users\mian_amm\Dropbox\Thèse\Data\SDMS', 'FP0120', true, [100, 500, 1100, 1500]);
%     Image = Image(:,:,1);
    
%% Target Creation

    N_k = 5;
    N_theta = 5;
    
    steeringVector1 = randn(N_theta,N_k) + 1i*randn(N_theta,N_k);
%     steeringVector = ones(N_theta,N_k);
%     steeringVector = reshape([1:13,12:-1:1], N_theta,N_k);
    Target1 = struct;
    Target1.position = [251, 251];
    Target1.Args = {steeringVector1, N_k, N_theta};
    Target1.SNR_dB = 50;
    Target1.function = @steeringVectorTarget;
    
    steeringVector2 = randn(N_theta,N_k) + 1i*randn(N_theta,N_k);
    Target2 = struct;
    Target2.position = [261, 261];
    Target2.Args = {steeringVector2, N_k, N_theta};
    Target2.SNR_dB = 70;
    Target2.function = @steeringVectorTarget;


    TargetStruct = struct;
    TargetStruct.list = [Target1];
    TargetStruct.number = length(TargetStruct.list);
    TargetStruct.Mx = 20;
    TargetStruct.My = 20;
    
    ImageTargets = createTargets(Image, Header, TargetStruct);
 
    % Steering vector
    p = steeringVector1.';
    p = p(:);
    p = p/sqrt(p'*p);

%% Decomposition
    Decomposition = struct;
    Decomposition.method = 'Exponential';
    Decomposition.J = 1;
    Decomposition.numberOfBands = N_k;
    Decomposition.numberOfLooks = N_theta;
    Decomposition.d1 = inf;
    Decomposition.d2 = inf;
    Decomposition.visualise = false;
    [ImageFTheta, k, theta] = decomposeFTheta(ImageTargets+Image, Header, Decomposition);

%     figure
%     plot(abs(p));
%     figure
%     plot(abs(reshape(ImageFTheta(gaussTarget.position(2), gaussTarget.position(1), :),N_k*N_theta,1)))
%     
    % Decimation
    ImageFTheta_decime = ImageFTheta(1:N_theta:end, 1:N_k:end,:);
    
%% Detectors
    DetectorsStruct = initTargetDetectors();
    Results = computeTestOnImage(ImageFTheta_decime, p, DetectorsStruct);
    thresholds = [2.202285750664898e+02;0.319801970816285;4.761946266981409e+06]; % For Pfa=10^-3
%% Plotting
    
% To visualize the image with decomposition
%     I = reshape(20*log10(abs(ImageFTheta)), Header.RgCnt, Header.AzCnt, N_k, N_theta);
%     essai_i4d(permute(I,[2,1,4,3]), Header.AzPixelSz, Header.AzCnt * Header.AzPixelSz, Header.RgPixelSz, Header.RgPixelSz*Header.RgCnt,...
%               k(1)*3e8/2/1.e9,k(N_k)*3e8/2/1.e9,theta(1)*180/pi,theta(N_theta)*180/pi,...
%               'transverse y (m)','radial x (m)','SAR Image (x,y)','f (GHz)',' \theta (deg)','Angular/Frequency behavior (f,\theta)'); 
%     
% 
    Amax = max([20*log10(abs(Image(:))); 20*log10(abs(ImageTargets(:)))]);
    Dyn_dB = 70;
    figure
    toShow = 20*log10(abs(Image));
    imagesc(Header.x, Header.y, toShow, [Amax - Dyn_dB,  Amax])
    axis image
    xlabel('x (m)', 'interpreter', 'latex')
    ylabel('y (m)', 'interpreter', 'latex')
    colorbar
    title('Image without targets', 'interpreter', 'latex')  
    colormap('bone')
    
    figure
    toShow = abs(fft2(ImageTargets));
    toShow = toShow*10^2/sqrt(sum(abs(toShow(:).^2)));
    imagesc(Header.kx, Header.ky, toShow)
    xlabel('$\mathbf{k}_x$ ($m^{-1}$)', 'interpreter', 'latex')
    ylabel('$\mathbf{k}_y$ ($m^{-1}$)', 'interpreter', 'latex')
    colorbar
    title('Spectro-Angular Behaviour of Targets', 'interpreter', 'latex')
    colormap('bone')
    
    figure
    toShow = abs(steeringVector1)/sqrt(sum(abs(steeringVector1(:).^2)));
    imagesc(k, theta, toShow)
    xlabel('$\|\mathbf{k}\|$ ($m^{-1}$)', 'interpreter', 'latex')
    ylabel('$\theta$ (rad)', 'interpreter', 'latex')
    colorbar
    title('Steering Vector', 'interpreter', 'latex')
    colormap('bone')
    
    figure
    toShow = 20*log10(abs(ImageTargets + Image));
    imagesc(Header.x, Header.y, toShow, [Amax - Dyn_dB,  Amax])
    axis image
    xlabel('x (m)', 'interpreter', 'latex')
    ylabel('y (m)', 'interpreter', 'latex')
    colorbar
    title('Image with targets', 'interpreter', 'latex')
    for tt=1:TargetStruct.number
        pos = TargetStruct.list(tt).position;
        circle(Header.x(pos(1)), Header.y(pos(2)), 3);
    end
    colormap('bone')
    
    figure
    toShow = 20*log10(abs(Image(1:N_theta:end, 1:N_k:end)));
    imagesc(Header.x, Header.y, toShow, [Amax - Dyn_dB,  Amax])
    axis image
    xlabel('x (m)', 'interpreter', 'latex')
    ylabel('y (m)', 'interpreter', 'latex')
    colorbar
    title('Image without targets decimated', 'interpreter', 'latex')
    colormap('bone')
    
    figure
    toShow = 20*log10(abs(sum(ImageFTheta_decime,3)));
    imagesc(Header.x, Header.y, toShow, [Amax - Dyn_dB,  Amax])
    axis image
    xlabel('x (m)', 'interpreter', 'latex')
    ylabel('y (m)', 'interpreter', 'latex')
    colorbar
    title('Image with targets decimated', 'interpreter', 'latex')
    for tt=1:TargetStruct.number
        pos = TargetStruct.list(tt).position;
        circle(Header.x(pos(1)), Header.y(pos(2)), 3);
    end   
   colormap('bone')
    


    
    for d=1:DetectorsStruct.number
        s = size(ImageFTheta_decime);
        toShow = zeros(s(1),s(2));
        [Mx, My] = size(DetectorsStruct.mask);
        toShow(floor(My/2)+1:end-floor(My/2), floor(Mx/2)+1:end-floor(Mx/2)) = Results{d};
        Amax = max(Results{d}(:));
        figure
        imagesc(Header.x, Header.y, toShow, [thresholds(d), Amax])
        colorbar
        title(DetectorsStruct.list(d).name, 'interpreter', 'latex')
        axis image
        for tt=1:TargetStruct.number
            pos = TargetStruct.list(tt).position;
            circle(Header.x(pos(1)), Header.y(pos(2)), 3);
        end
        xlabel('x (m)', 'interpreter', 'latex')
        ylabel('y (m)', 'interpreter', 'latex')
        colormap('bone')
    end
%     
%     for d=1:DetectorsStruct.number
%         s = size(ImageFTheta_decime);
%         [Mx, My] = size(DetectorsStruct.mask);
%         toShow = permute(Results{d}, [2,1]);
%         Amax = max(Results{d}(:));
%         Affiche_map(Header.x(floor(Mx/2)*N_theta:N_theta:end-floor(Mx/2)*N_theta),Header.y(floor(My/2)*N_k:end-floor(My/2)*N_k), 10*log10(toShow), DetectorsStruct.list(d).name)
%         axis image
%         for tt=1:TargetStruct.number
%             pos = TargetStruct.list(tt).position;
%             circle(Header.y(pos(2)), Header.x(pos(1)), 3);
%         end
%     end
%% ------------- END CODE ----------------
