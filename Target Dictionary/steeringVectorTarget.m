%% Description
% Target with given steering vectorGeneration
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
function ImageTarget = steeringVectorTarget(x, y, kx, ky, position, Args)
% ImageTarget - Create Gaussian Image Target
% Syntax:  ImageTarget = gaussianTarget(x, y, kx, ky, position, Args)
%
% Inputs:
%    - x, y = x and y vectors
%    - kx, ky = kx and ky vectors
%    - position = [xtarget, ytarget] position of the target on the image
%    - Args = {p, N_k, N_theta} = steering vector (size N_theta*N_k), number of Bands, number
%      of Looks
% Outputs:
%    ImageTarget - image of the target
%

    p = Args{1};
    N_k = Args{2};
    N_theta = Args{3};
    
    % Polar Coordinates
	[KX,KY] = meshgrid(kx, ky);
	k = sqrt(KY.^2+KX.^2);
	theta = atan2(KY, KX);
    
    % Target position
    xt = x(position(1));
    yt = y(position(2));
    phaseTerm = exp(-1i*2*pi*k.*((xt-x(1))*cos(theta)+(yt-y(1))*sin(theta)));

    % Decimate vectors in order to have the desired number of sub-bands
	kmin = min(k(:));
	kmax = max(k(:));
	theta_min = min(theta(:));
	theta_max = max(theta(:));
	sigma_t = (theta_max - theta_min) / N_theta;
	sigma_k = (kmax - kmin) / N_k;
	theta_0 = theta_min + (0:N_theta-1)*sigma_t + sigma_t/2;
	k_0 = kmin + (0:N_k-1)*sigma_k + sigma_k /2;
    
    
    % Steering vector
%     figure
    Spectre = zeros(size(k));
    for k_theta = 1:N_theta
        for k_k = 1:N_k
                R =  (abs(k-k_0(k_k)) <= sigma_k/2) .* (abs(theta-theta_0(k_theta)) <= sigma_t/2);
                Spectre = Spectre + R*p(k_theta,k_k);    
%                 imagesc(abs(R*p(k_theta,k_k)))
%                 pause
        end
    end

    ImageTarget = ifft2(Spectre.*phaseTerm);   

end

%% ------------- END CODE ----------------
