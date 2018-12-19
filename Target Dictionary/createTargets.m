%% Description
% Creation of an image of several targets
%% Specifiactions
% Other m-files required: computeImagePower.m
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

function ImageTargets = createTargets(Image, Header, TargetStruct)
% ImageTargets - Create several targets
% Syntax:  ImageTargets = createTargets(Image, Header, TargetStruct)
%
% Inputs:
%    - Image = the image on which we add the targets
%    - Header = the Header of the image
%    - TargetStruct = Specifications of targets
%
% Outputs:
%    ImageTargets - image of the targets

    x = Header.x;
    y = Header.y;
    kx = Header.kx;
    ky = Header.ky;
    ImageTargets = zeros(size(Image));
    for t=1:TargetStruct.number
        target = TargetStruct.list(t);
        imagePower = computeImagePower(Image, target.position, TargetStruct.Mx, TargetStruct.My);
        Temp = target.function(x, y, kx, ky, target.position, target.Args);
        Temp = Temp / sqrt(sum(abs(Temp(:).^2))) * imagePower * 10.^(target.SNR_dB/20);
        ImageTargets = ImageTargets + Temp;
    end
end

%% ------------- END CODE ----------------
