%% Description
% Simulation with two close targets : one gaussian and one with steering vector
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
%     close all
    set(0,'defaulttextinterpreter','latex')
    
%% Reading Image
    addpath('..', '../Hyperimages', '../Target Dictionary', '../Detectors')
    imagePath = '../../../../../Data/SANDIA/MiniSAR20050519p0010image002.mat';
    [Image, Header] = readSandia(imagePath, true, [1, 1, 301, 501]);
%     [Image, Header] = readSDMS('C:\Users\mian_amm\Dropbox\Thèse\Data\SDMS', 'FP0120', true, [1, 1, 500, 500]);
    Image = Image(:,:,1);
        
%% Decomposition params    
    Decomposition = struct;
    Decomposition.method = 'Exponential';
    Decomposition.J = 1;
    Decomposition.numberOfBands = 5;
    Decomposition.numberOfLooks = 5;
    N_k = Decomposition.numberOfBands^Decomposition.J;
    N_theta = Decomposition.numberOfLooks^Decomposition.J;
    Decomposition.d1 = 10;
    Decomposition.d2 = 10;
    Decomposition.visualise = false;
    Decomposition.dyn_dB = 50;
    
%% Target Creation  
    
    Target1 = struct;
    Target1.position = [251, 151];
    mean = [110, -0.005];
    scale = [1.2, 0.012];
    Target1.Args = {mean, scale};
    Target1.SNR_dB = 70;
    Target1.function = @gaussianTarget;
    Target1.toCircle = false;
    
%     Target2 = struct;
%     Target2.position = [251, 201];
%     mean = [114, 0.005];
%     scale = [1.2, 0.012];
%     Target2.Args = {mean, scale};
%     Target2.SNR_dB = 22;
%     Target2.function = @gaussianTarget;
%     Target2.toCircle = true;
%     
%     % Steering vector
%     TargetStruct = struct;
%     TargetStruct.list = [Target2];
%     TargetStruct.number = length(TargetStruct.list);
%     TargetStruct.Mx = 20;
%     TargetStruct.My = 20;
%     ImageTargets = createTargets(Image, Header, TargetStruct);
%     steeringVector = spectreToSteeringVector(fft2(ImageTargets), Header.kx, Header.ky, Decomposition);
    
    N_k = 5;
    N_theta = 5;
    randn('seed',0)
    steeringVector = randn(N_theta,N_k) + 1i*randn(N_theta,N_k);
    Target2 = struct;
    Target2.position = [251, 181];
    Target2.Args = {steeringVector, N_k, N_theta};
    Target2.SNR_dB = 15;
    Target2.function = @steeringVectorTarget;
    Target2.toCircle = true;
    
    % Iamge of targets
    TargetStruct = struct;
    TargetStruct.list = [Target1, Target2];
    TargetStruct.number = length(TargetStruct.list);
    TargetStruct.Mx = 20;
    TargetStruct.My = 20;
    ImageTargets = createTargets(Image, Header, TargetStruct);

    
    
    
    
  

%% Decomposition
    [ImageFTheta, k, theta] = decomposeFTheta(ImageTargets + Image, Header, Decomposition);
  
    % Decimation
    ImageFTheta_decime = ImageFTheta(1:N_theta:end, 1:N_k:end,:);
    
%% Steering vector
    p = steeringVector.';
    p = p(:);
    p = p/sqrt(p'*p);    
    
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
    Dyn_dB = 70;
    figure
    toShow = 20*log10(abs(fft2(ImageTargets)));
    Amax = max(toShow(:));
    imagesc(Header.kx, Header.ky, toShow, [Amax - Dyn_dB,  Amax])
    xlabel('$\mathbf{k}_x$ ($m^{-1}$)', 'interpreter', 'latex')
    ylabel('$\mathbf{k}_y$ ($m^{-1}$)', 'interpreter', 'latex')
    colorbar
    title('Spectro-Angular Behaviour of Targets', 'interpreter', 'latex')
    colormap('bone')
    
    Amax = max([20*log10(abs(Image(:))); 20*log10(abs(ImageTargets(:)))]);
    
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
    toShow = abs(steeringVector);
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
        if TargetStruct.list(tt).toCircle
            pos = TargetStruct.list(tt).position;
            circle(Header.x(pos(1)), Header.y(pos(2)), 1.5);
        end
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
        if TargetStruct.list(tt).toCircle
            pos = TargetStruct.list(tt).position;
            circle(Header.x(pos(1)), Header.y(pos(2)), 1.5);
        end
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
%         colorbar
        title(DetectorsStruct.list(d).name, 'interpreter', 'latex')
        axis image
        for tt=1:TargetStruct.number
            if TargetStruct.list(tt).toCircle
                pos = TargetStruct.list(tt).position;
                circle(Header.x(pos(1)), Header.y(pos(2)), 1.5);
            end
        end   
        xlabel('x (m)', 'interpreter', 'latex')
        ylabel('y (m)', 'interpreter', 'latex')
        colormap('bone')
    end
    
%     for d=1:DetectorsStruct.number
%         s = size(ImageFTheta_decime);
%         [Mx, My] = size(DetectorsStruct.mask);
%         toShow = permute(Results{d}, [2,1]);
%         Amax = max(Results{d}(:));
%         Affiche_map(Header.x(floor(Mx/2)*N_theta:N_theta:end-floor(Mx/2)*N_theta),Header.y(floor(My/2)*N_k:end-floor(My/2)*N_k), 10*log10(toShow), DetectorsStruct.list(d).name)
%         for tt=1:TargetStruct.number
%             if TargetStruct.list(tt).toCircle
%                 pos = TargetStruct.list(tt).position;
%                 circle(Header.y(pos(2)), Header.x(pos(1)), 3);
%             end
%         end   
%     end
%% ------------- END CODE ----------------
