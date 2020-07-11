# OptMet - Robust optimal design of 3D printed mechanical metamaterials

The Matlab-based software `OptMet` provides a rational approach for the robust and optimal design of 3D printed metamaterial beams. Uncertainties derived from the additive manufacturing process are rigorously quantified. Variability about the performance function is also considered, thus conferring the optimal design with robustness. The optimal design looks for maximizing vibration attenuation considering different number for added resonators with different masses.

## Example

As an example of the output of `OptMet`, some files are contained in this repository assuming: 

* The frequency range in which the attenuation is to be maximized is:  [280, 380] Hz
* Discretization of the mass variable of the vibration attenuation: 0.05 : 0.05 : 1 (ratio of the sum of resonator massesto the mass of the beam)
* Number of samples from the material properties PDFs (Young's modulus E ~ N(1621.7, 49.9^2) and Density ~ N(948.9, 7.4^2) ): 100
* Trade-off variable between expectation and variance of the performance function (A): 0.5
* Maximum number of resonators: 15
* Cost function defined with `pchip` Matlab function and the following interpolating points: (1,0), (8,0.3), (12,0.65) (15,0.9)

The software performs an exhaustive search in the discretrized mass variable for each number of resonators, up to the maximum number. Therefore, the optimal mass is obtained for each of these number of resonators. The output of each of these sub-steps are stored in the folder `res`. For example, for 4 resonators:

<img src="https://github.com/scanteroch/OptMet/blob/master/example_png/N_res4.png" alt="drawing" width="600"/>

## References

The `OptMet` Matlab-based software is based on the following references:

> H. Meng, D. Chronopoulos, A. T. Fabro, W. Elmadih, I. Maskery, Rainbow metamaterials for broadband multi-frequency vibration attenua-350tion: Numerical analysis and experimental validation, [Journal of Sound and Vibration 465 (2020) 115005](https://doi.org/10.1016/j.jsv.2019.115005).

>H. Meng, D. Chronopoulos, A. T. Fabro, I. Maskery, Y. Chen, Optimal design of rainbow elastic metamaterials, [International Journal of Mechanical Sciences (2019) 105185](https://doi.org/10.1016/j.ijmecsci.2019.105185).

>H. Meng,  D. Chronopoulos,  A. T. Fabro,  Numerical simulation data for the dynamic properties of rainbow metamaterials,  [Data in Brief (2019) 104772](https://doi.org/10.1016/j.dib.2019.104772).

## Acknowledgements

The authors would like to acknowledge the support acquired by the H2020 DiaMoND project (Grant Agreement ID:785859), FAPESP Thematic Grant ENVIBRO (Grant Agreement ID: 2018/15894-0), the Brazilian National Council of Research CNPq (Grant Agreement ID: 420304/2018-5), the Science and Technology Development Fund, Macau SAR (File no. SKL-IOTSC-2018-2020) and the Start-up Research Grant of University of Macau (File No. SRG2019-00194-IOTSC). 
