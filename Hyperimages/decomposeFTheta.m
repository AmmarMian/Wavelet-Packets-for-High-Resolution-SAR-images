%% Description
% Decomposing a SAR image through (f, theta) decomposition 
%% Specifiactions
% Other m-files required: gbellmfunction.m
% MAT-files required: none
%% Authors
% Authors: Ammar Mian, Ph.D., Signal processing at Supelec SONDRA
% Jean-Phillipe Ovarlez, Researcher at Supelec SONDRA
% Guillaume Ginolhac, Researcher at LISTIC Annecy
% Abdourrahmane M. ATTO, Researcher at LISTIC Annecy
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
function [ImageFTheta, k_vec, theta_vec, handleFTheta, handleXY, handleWavelet] = decomposeFTheta(Image, Header, Decomposition)
% decomposeFTheta - Do the decomposition in (f, theta) of the image and add
% a target
% Syntax:  [ImageFTheta, k, theta] = decomposeFTheta(Image, Header, Decomposition)
%
% Inputs:
%    - Image = The Image
%    - Header = A strcut with fields
%     * Nx, Ny = Size of image
%     * x, y = x and y vectors
%     * kx, ky = kx and ky vectors
%    - Decomposition params = A strcut with fields
%     * method = 'Shannon' or 'Exponential'
%     * d1, d2 = Exponential wavelet shape parameters
%     * J = decomposition level
%     * numberOfBands, numberOfLooks = number of divisions per level
%     * visualise = boolean for visualisation 
%
% Outputs:
%    ImageFTheta - The multivariate Image
%	 k_vec - The module vector in polar coordinates
%	 theta_vec - The angles vector in polar coordinates
%    handles - for visualisation
%

   
    % Some variables
    kx = Header.kx;
    ky = Header.ky;
    Nx = Header.Nx;
    Ny = Header.Ny;
    method = Decomposition.method;
    d1 = Decomposition.d1;
    d2 = Decomposition.d2;
    J = Decomposition.J;
    numberOfBands = Decomposition.numberOfBands;
    numberOfLooks = Decomposition.numberOfLooks;
    N_k = numberOfBands^J; % Number of bands k (module of scattering vector)
	N_theta = numberOfLooks^J; % Number of bands theta (angle of scattering vector)
    
    % Going to frequencies domains
    Kxy = fft2(Image);
    
    % Polar Coordinates
	[KX,KY] = meshgrid(kx, ky);
	k = sqrt(KY.^2+KX.^2);
	theta = atan2(KY, KX);

   % Decimate vectors in order to have the desired number of sub-bands
	kmin = min(k(:));
	kmax = max(k(:));
	theta_min = min(theta(:));
	theta_max = max(theta(:));
	sigma_t = (theta_max - theta_min) / N_theta;
	sigma_k = (kmax - kmin) / N_k;
	theta_vec = theta_min + (0:N_theta-1)*sigma_t + sigma_t/2;
	k_vec = kmin + (0:N_k-1)*sigma_k + sigma_k /2;

    % Initialise the multivariate images
    ImageFTheta = zeros(Ny, Nx, N_k * N_theta);

    % Doing the decomposition
    RadiusBell_k = sigma_k/(2*numberOfBands^(J-1));
    RadiusBell_theta_k = sigma_t/(2*numberOfLooks^(J-1));
    if ~strcmp(method, 'Shannon')
        if d1==inf
            method = 'Shannon';
        else
            method = 'Exponential';
        end
    end
    if Decomposition.visualise
        cmap = 'gray';
        Dyn_dB = Decomposition.dyn_dB;
        figure;
        Amax = max(max(20*log10(abs(Image))));
        imagesc(Header.x, Header.y, 20*log10(abs(Image)), [Amax-Dyn_dB, Amax])
        xlabel('$x$', 'interpreter', 'latex')
        ylabel('$y$', 'interpreter', 'latex')
        title('Image', 'interpreter', 'latex')
        set(gca, 'ydir', 'normal')
        axis image
        colormap(cmap)
        figure;
        Amax = max(max(20*log10(abs(Kxy))));
        imagesc(Header.kx, Header.ky, 20*log10(abs(Kxy)), [Amax-Dyn_dB, Amax])
        xlabel('$k_x$', 'interpreter', 'latex')
        ylabel('$k_y$', 'interpreter', 'latex')
        title('Spectre of image', 'interpreter', 'latex')
        set(gca, 'ydir', 'normal')
        axis image
        colormap(cmap)
        handleFTheta = figure('Position', [100, 100, 600, 400]);
        handleXY = figure('Position', [800, 100, 600, 400]);
        handleWavelet = figure('Position', [1500, 100, 600, 400]);
        handlesaxes = [];
    end
    for k_k = 1:N_k
        for k_theta = 1:N_theta
            % Projection / angular and spectral domains
            if ~strcmp(method,'Shannon')  % Methode Wavelet
                % Computation ---> FFT wavelet weights
                CenterBell_k = k_vec(k_k); %
                CenterBell_theta_k = theta_vec(k_theta);
                Hjmnk = gbellmfunction(k, [RadiusBell_k d1 CenterBell_k]);
                Hjmnt = gbellmfunction(theta, [RadiusBell_theta_k d2 CenterBell_theta_k]);
                Wavelet = Hjmnk.*Hjmnt;
                A = Wavelet.*Kxy;
                ImageFTheta(:,:,k_k+N_k*(k_theta-1)) = ifft2(A);
                

            else % Method Rectangular Window
                Wavelet = (abs(k-k_vec(k_k)) <= sigma_k/2) .* (abs(theta-theta_vec(k_theta)) <= sigma_t/2);
                A = Kxy .*  Wavelet;
                ImageFTheta(:,:,k_k+N_k*(k_theta-1)) = ifft2(A);

            end
            
            if Decomposition.visualise
                
                if d1==inf
                    chard1 = '\infty';
                else
                    chard1 = num2str(d1);
                end
                
                if d2==inf
                    chard2 = '\infty';
                else
                    chard2 = num2str(d2);
                end
                
                figure(handleWavelet);
                subplot(N_k, N_theta, k_theta+N_theta*(k_k-1))
                imagesc(Header.kx, Header.ky, abs(Wavelet))
                title(sprintf('$ \\Psi^{\\mathbf{S},[%s,%s]}_{%d,[%d,%d]} $', chard1, chard2, Decomposition.J, k_k, k_theta), 'interpreter', 'latex')
                xlabel('$k_x$', 'interpreter', 'latex')
                ylabel('$k_y$', 'interpreter', 'latex')
                set(gca, 'ydir', 'normal')
                axis image
                colormap(cmap)
                
                figure(handleFTheta);
                subplot(N_k, N_theta, k_theta+N_theta*(k_k-1))
                Amax = max(max(20*log10(abs(Kxy))));
                imagesc(Header.kx, Header.ky, 20*log10(abs(A)), [Amax-Dyn_dB, Amax])
                title(sprintf('$\\Psi^{\\mathbf{S},[%s,%s]}_{%d,[%d,%d]}\\times\\tilde{I}$', chard1, chard2, Decomposition.J, k_k, k_theta), 'interpreter', 'latex')
                xlabel('$k_x$', 'interpreter', 'latex')
                ylabel('$k_y$', 'interpreter', 'latex')
                set(gca, 'ydir', 'normal')
                axis image
                colormap(cmap)
                
                figure(handleXY);
                axtemp = subplot(N_k, N_theta, k_theta+N_theta*(k_k-1));
                handlesaxes = cat(2,handlesaxes,axtemp);
                Amax = max(max(20*log10(abs(Image))));
%                 Amax = max(max(20*log10(abs(ImageFTheta(1:N_k:end,1:N_theta:end,k_k+N_k*(k_theta-1))))));
                imagesc(Header.x, Header.y, 20*log10(abs(ImageFTheta(:,:,k_k+N_k*(k_theta-1)))), [Amax-Dyn_dB, Amax])
                title(sprintf('$\\mathbf{C}^{[%s,%s]}_{%d,[%d,%d]}$', chard1, chard2, Decomposition.J, k_k, k_theta), 'interpreter', 'latex')
                xlabel('$x$', 'interpreter', 'latex')
                ylabel('$y$', 'interpreter', 'latex')
                set(gca, 'ydir', 'normal')
                axis image
                colormap(cmap)
            end
            
        end
    end
    
    if Decomposition.visualise
        linkaxes(handlesaxes,'xy')
    end
    
    % k, theta vectors
    k = k(1,:);
    theta = theta(:,1);
end

%% ------------- END CODE ----------------
