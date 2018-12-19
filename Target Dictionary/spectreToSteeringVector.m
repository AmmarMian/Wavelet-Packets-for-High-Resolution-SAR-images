%% Description
% Spectre to steering vector function
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
function p = spectreToSteeringVector(Spectre, kx, ky, Decomposition)
% spectreToSteeringVector - discretise a 2d spectrum for steering vector
% purposes
% Syntax:  p = spectreToSteeringVector(Spectre, kx, ky, N_k, N_theta)
%
% Inputs:
%    - Spectre = the spectre image
%    - kx, ky = kx and ky vectors
%    - N_k, N_theta = the number of bands and looks
%    - Decomposition params = A strcut with fields
%     * method = 'Shannon' or 'Exponential'
%     * d1, d2 = Exponential wavelet shape parameters
%     * J = decomposition level
%     * numberOfBands, numberOfLooks = number of divisions per level
%     * visualise = boolean for visualisation 
% Outputs:
%    p - a N_k*N_theta matrix
%

 
    method = Decomposition.method;
    d1 = Decomposition.d1;
    d2 = Decomposition.d2;
    J = Decomposition.J;
    numberOfBands = Decomposition.numberOfBands;
    numberOfLooks = Decomposition.numberOfLooks;
    N_k = numberOfBands^J; % Number of bands k (module of scattering vector)
	N_theta = numberOfLooks^J; % Number of bands theta (angle of scattering vector)
    
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
    RadiusBell_k = sigma_k/(2*numberOfBands^(J-1));
    RadiusBell_theta_k = sigma_t/(2*numberOfLooks^(J-1));
    
    p = zeros(N_theta, N_k);
    for k_theta = 1:N_theta
        for k_k = 1:N_k
            % Projection / angular and spectral domains
            if ~strcmp(method,'Shannon')  % Methode Wavelet
                % Computation ---> FFT wavelet weights
                CenterBell_k = k_vec(k_k); %
                CenterBell_theta_k = theta_vec(k_theta);
                Hjmnk = gbellmfunction(k, [RadiusBell_k d1 CenterBell_k]);
                Hjmnt = gbellmfunction(theta, [RadiusBell_theta_k d2 CenterBell_theta_k]);
                Wavelet = Hjmnk.*Hjmnt;
                A = Wavelet.*Spectre;
            else % Method Rectangular Window
                Wavelet = (abs(k-k_vec(k_k)) <= sigma_k/2) .* (abs(theta-theta_vec(k_theta)) <= sigma_t/2);
                A = Spectre .*  Wavelet;
            end
            mask = isnan(A);
            p(k_theta, k_k) = sum(sum(abs(A(~mask)).^2));
        end
    end

end

%% ------------- END CODE ----------------
