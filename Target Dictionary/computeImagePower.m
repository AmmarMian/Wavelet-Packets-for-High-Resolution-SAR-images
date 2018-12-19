%% Description
% Estimate Image ppower on a window around a position
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

function imagePower = computeImagePower(Image, position, Mx, My)
% computeImagePower - Estimate Image ppower on a window around a position
% Syntax:  imagePower = computeImagePower(Image, pos_x, pos_y, N)
%
% Inputs:
%    - Image = the Image
%    - position = posiiton in pixels
%    - Mx, My = Window length for estimation
%       
% Outputs:
%    imagePower - the estimated power
%
    
    pos_y = position(1);
    pos_x = position(2);
    imagePower = sqrt(mean(mean(abs(Image(pos_y-floor(My/2):pos_y+floor(My/2), pos_x-floor(Mx/2):pos_x+floor(Mx/2)).^2))));
end

%% ------------- END CODE ----------------
