%% Description
% Simulation with a target with gaussian spectro-angular behaviour
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
    [Image, Header] = readSandia(imagePath, false, [1, 1, 501, 501]);
%     [Image, Header] = readSDMS('C:\Users\mian_amm\Dropbox\Thèse\Data\SDMS', 'FP0120', true, [100, 500, 1100, 1500]);
    Image = Image(:,:,1);
    
%% Target Creation  
    
    Target1 = struct;
    Target1.position = [51, 51];
    mean = [112, 0];
    scale = [1, 0.01];
    Target1.Args = {mean, scale};
    Target1.SNR_dB = 50;
    Target1.function = @gaussianTarget;
    


    TargetStruct = struct;
    TargetStruct.list = [Target1];
    TargetStruct.number = length(TargetStruct.list);
    TargetStruct.Mx = 20;
    TargetStruct.My = 20;
    
    ImageTargets = createTargets(Image, Header, TargetStruct);
 
  

%% Decomposition
    Decomposition = struct;
    Decomposition.method = 'Exponential';
    Decomposition.J = 1;
    Decomposition.numberOfBands = 5;
    Decomposition.numberOfLooks = 5;
    N_k = Decomposition.numberOfBands^Decomposition.J;
    N_theta = Decomposition.numberOfLooks^Decomposition.J;
    Decomposition.d1 = inf;
    Decomposition.d2 = inf;
    Decomposition.visualise = true;
    Decomposition.dyn_dB = 50;
    [ImageFTheta, k, theta] = decomposeFTheta(Image + ImageTargets, Header, Decomposition);
  
    % Decimation
    ImageFTheta_decime = ImageFTheta(1:N_theta:end, 1:N_k:end,:);
    
%% Steering vector
    steeringVector1 = spectreToSteeringVector(fft2(ImageTargets), Header.kx, Header.ky, Decomposition);
    p = steeringVector1.';
    p = p(:);
    p = p/sqrt(p'*p);    
    
%% Detectors
    DetectorsStruct = initTargetDetectors();
%     Results = computeTestOnImage(ImageFTheta_decime, p, DetectorsStruct);
    tic
   
    [Mx, My] = size(DetectorsStruct.mask);
    [Ny, Nx, ~] = size(ImageFTheta_decime);
    nbx = 6;
    nby = 4;
    numberOfWorkers = nbx*nby;
    ResParfor = cell(numberOfWorkers,1);
    SubImages = cell(numberOfWorkers,1);
    for index=1:numberOfWorkers
        %2d indexes for image
        subIntervalx = floor(linspace(0,Nx-Mx,nbx+1)) + floor(Mx/2)+1;
        subIntervaly = floor(linspace(0,Ny-My,nby+1)) + floor(My/2)+1;
        [i_y, i_x] = ind2sub([nby, nbx], index);
        indexesY = (subIntervaly(i_y) - floor(My/2):subIntervaly(i_y+1) + floor(My/2));
        indexesX = (subIntervalx(i_x) - floor(Mx/2):subIntervalx(i_x+1) + floor(Mx/2));

        % Transfering the needed data for the sub-image
        SubImages{index} = ImageFTheta_decime(indexesY, indexesX,:);
    end 
    
    
    % Compute results
%     hbar = parfor_progressbar(numberOfWorkers,'Computing...');  %create the progress bar
    parfor index=1:numberOfWorkers
        ResParfor{index} = computeTestOnImage(SubImages{index} , p, DetectorsStruct);
%         hbar.iterate(1);   % update progress by one iteration
    end
%     close(hbar);   %close progress bar

    Results = cell(DetectorsStruct.number,1);
    for d=1:DetectorsStruct.number
        ResTemp = zeros(Ny-My, Nx-Mx);
        for index=1:numberOfWorkers
            [i_y, i_x] = ind2sub([nby, nbx], index);
            subIntervalx = floor(linspace(0,Nx-Mx,nbx+1));
            subIntervaly = floor(linspace(0,Ny-My,nby+1));
            indexesY = (subIntervaly(i_y)+1:subIntervaly(i_y+1)+1);
            indexesX = (subIntervalx(i_x)+1:subIntervalx(i_x+1)+1);
            ResTemp(indexesY,indexesX) = ResParfor{index}{d};
        end
        Results{d} = ResTemp;
    end
 
    toc
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
    colormap('hot')
    
    figure
    toShow = abs(fft2(ImageTargets));    imagesc(Header.kx, Header.ky, toShow)
    xlabel('$\mathbf{k}_x$ ($m^{-1}$)', 'interpreter', 'latex')
    ylabel('$\mathbf{k}_y$ ($m^{-1}$)', 'interpreter', 'latex')
    colorbar
    title('Spectro-Angular Behaviour of Targets', 'interpreter', 'latex')
    colormap('hot')

    
    figure
    toShow = abs(steeringVector1);
    imagesc(k, theta, toShow)
    xlabel('$\|\mathbf{k}\|$ ($m^{-1}$)', 'interpreter', 'latex')
    ylabel('$\theta$ (rad)', 'interpreter', 'latex')
    colorbar
    title('Steering Vector', 'interpreter', 'latex')
    colormap('hot')
    
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
    colormap('hot')
    
    figure
    toShow = 20*log10(abs(Image(1:N_theta:end, 1:N_k:end)));
    imagesc(Header.x, Header.y, toShow, [Amax - Dyn_dB,  Amax])
    axis image
    xlabel('x (m)', 'interpreter', 'latex')
    ylabel('y (m)', 'interpreter', 'latex')
    colorbar
    title('Image without targets decimated', 'interpreter', 'latex')
    colormap('hot')
    
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
   colormap('hot')
    


    
    for d=1:DetectorsStruct.number
        s = size(ImageFTheta_decime);
        toShow = zeros(s(1),s(2));
        [Mx, My] = size(DetectorsStruct.mask);
        toShow(floor(My/2)+1:end-floor(My/2), floor(Mx/2)+1:end-floor(Mx/2)) = Results{d};
        Amax = max(Results{d}(:));
        figure
        imagesc(Header.x, Header.y, toShow)%, [DetectorsStruct.list(d).threshold, Amax])
        colorbar
        title(DetectorsStruct.list(d).name, 'interpreter', 'latex')
        axis image
        for tt=1:TargetStruct.number
            pos = TargetStruct.list(tt).position;
            circle(Header.x(pos(1)), Header.y(pos(2)), 3);
        end
        xlabel('x (m)', 'interpreter', 'latex')
        ylabel('y (m)', 'interpreter', 'latex')
    end
    
    for d=1:DetectorsStruct.number
        s = size(ImageFTheta_decime);
        [Mx, My] = size(DetectorsStruct.mask);
        toShow = permute(Results{d}, [2,1]);
        Amax = max(Results{d}(:));
        Affiche_map(Header.x(floor(Mx/2)*N_theta:N_theta:end-floor(Mx/2)*N_theta),Header.y(floor(My/2)*N_k:end-floor(My/2)*N_k), 10*log10(toShow), DetectorsStruct.list(d).name)
        axis image
        for tt=1:TargetStruct.number
            pos = TargetStruct.list(tt).position;
            circle(Header.y(pos(2)), Header.x(pos(1)), 3);
        end
    end
%% ------------- END CODE ----------------
