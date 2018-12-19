# Wavelet packet adapted to High-resolution SAR images

The work here correspond to the simulation done for the TGRS paper on Wavelets:

	"Design of New Wavelet Packets Adapted to High-Resolution SAR Images with an Application to Target Detection",
	Ammar Mian, Jean-Philippe Ovarlez , Abdourahmane M. Atto, Guillaume Ginolhac
	in Transactions on Geoscience and Remote Sensing
	Preprint available at: https://ammarmian.github.io/publication/tgrs-2019/

If you use any code of the repository, please consider citing the above mentionned reference.

## Files' organisation

This folder is organised as follows:
- **Target Detection/**: Main scripts generating figures presented in paper.
- **Target Dictionary/**: Functions allowing to create targets with given spectro-angular behaviour.
- **HyperImages/**: Wavelet decomposition functions.
- **Detectors/**: Functions for target detection.

## Requirements 

Matlab 2014 or above.

Datasets available at:
- https://www.sdms.afrl.af.mil/index.php?collection=ccd_challenge
- https://www.sandia.gov/radar/complex-data/

## Credits
**Author:** Ammar Mian, Ph.d student at SONDRA, CentraleSupélec
 - **E-mail:** ammar.mian@centralesupelec.fr
 - **Web:** https://ammarmian.github.io/
 
 Acknowledgements to:
 - [**Guillaume Ginolhac**](https://www.listic.univ-smb.fr/presentation/membres/enseignants-chercheurs/guillaume-ginolhac/), LISTIC, Université Savoie Mont-Blanc
 - [**Jean-Philippe Ovarlez**](http://www.jeanphilippeovarlez.com/), DEMR, ONERA , Université Paris-Saclay  & SONDRA, CentraleSupélec
 - [**Abdourahmane M. Atto**](http://am.atto.free.fr/), LISTIC, Université Savoie Mont-Blanc

 
## Copyright
 
 Copyright 2018 @CentraleSupelec

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.