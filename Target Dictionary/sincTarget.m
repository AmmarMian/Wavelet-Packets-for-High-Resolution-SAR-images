%% Description
% sinus cardinal Target Generation
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
function ImageTarget = sincTarget(x, y, kx, ky, position, Args)
% ImageTarget - Create Gaussian Image Target
% Syntax:  ImageTarget = sincTarget(x, y, kx, ky, position, Args)
%
% Inputs:
%    - x, y = x and y vectors
%    - kx, ky = kx and ky vectors
%    - position = [xtarget, ytarget] position of the target on the image
%    - Args = {mean, scale} = position of sinc, scale of sinc
%
% Outputs:
%    ImageTarget - image of the target
%   
    
    mean = Args{1};
    scale = Args{2};
    
    % Polar Coordinates
	[KX,KY] = meshgrid(kx, ky);
	k = sqrt(KY.^2+KX.^2);
	theta = atan2(KY, KX);
    
    % Target position
    xt = x(position(1));
    yt = y(position(2));
    phaseTerm = exp(-1i*2*pi*k.*((xt-x(1))*cos(theta)+(yt-y(1))*sin(theta)));
    
    % Gaussian
    k_0 = mean(1);
    sigma_k = scale(1);
    theta_0 = mean(2);
    sigma_theta = scale(2);
    sincTerm = sinc( 2*(k-k_0)/(sigma_k) ) .* sinc( 2*(theta-theta_0)/(sigma_theta) );
    
    ImageTarget = ifft2(sincTerm.*phaseTerm);
end

%% ------------- END CODE ----------------
