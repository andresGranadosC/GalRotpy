
# GalRotpy

A Python3-based tool for parametrizing the rotation curve and the galaxy potential of disk-like galaxies.


GalRotpy allows to model the dynamical mass of disk-like galaxies. It makes a parametric fit of the rotation curve by means the composed gravitational potential of the galaxy. This software is presented as a pedagogical tool for supporting the following research areas: **galactic dynamics, cosmology and gravitation**.

GalRotpy requires essentially the [galpy](https://github.com/jobovy/galpy) package and the three fundamental Python3 packages: [Numpy](http://www.numpy.org/), [Scipy](https://www.scipy.org/) and [Matplotlib](http://matplotlib.org/). GalRotpy doesn't support Python2.


## Available gravitational potentials:

GalRotpy can give a first approximation of the galaxy rotation curve using the following gravitational potential schemas:

 * bulge model uning a Miyamoto-Nagai potential,
 * stellar or gaseous disk:
   * thin or thick disks implementing Miyamoto-Nagai potentials, and/or
   * an exponential disk model
 * and finally the Dark Halo with a Navarro-Frenk-White (NFW) or with a Burkert potential.

## How to use GalRotpy

GalRotpy can be used in the following ways:

```sh
$ python3 GalRotpy rot_curve.txt [bulge halo] | [bulge disk halo] | [disk halo]
```
or
```sh
$ python3 GalRotpy rot_curve.txt --guess=init_guess_params.txt
```


 * **rot_curve.txt**  it is mandatory to specify the file that contains the rotation curve. There must be three columns separated by tabulation character:
    * r (galactocentric distance in ![](https://latex.codecogs.com/svg.latex?kpc) )
    * vel (the circular speed in ![](https://latex.codecogs.com/svg.latex?km/s) at each distance r )
    * e_vel (the uncertainty for the circular speed in ![](https://latex.codecogs.com/svg.latex?km/s) at each distance r )

|r| vel|    e_vel| 
| ------ | ------ | ------ |
|   0.24|  37.3|   6.2|
|   0.28|	37.9|	5.5|
|   0.46|	47.1|	2.8|
|   0.73|	55.1|	3.3|
|   ...|    ...|	...|

* **[bulge halo] | [bulge disk halo] | [disk halo]** are the available options to make a first guess of the potential composition to reproduce the given rotation curve: bulge and halo or bulge and disk and halo or disk and halo
* **--guess=init_guess_params.txt** if you don't enter any of the options above, GalRotpy needs to specify the guess table where earch component will be used:

|component| mass|   a (kpc)|	b (kpc)|	checked|
| ------ | ------ | ------ | ------ | ------ |
|'BULGE'|110000000.0|0.0|0.495| True |
|'THIN DISK'|3900000000.0|5.3|0.25| False |
|'THICK DISK'|39000000000.0|2.6|0.8| True |
|'EXP DISK'|500.0|5.3|------| False |
|'DARK HALO'|140000000000.0|13.0|------| True |
|'BURKERT HALO'|8000000.0|20.0|------| False |

* where the **mass** column is in units of ![](https://latex.codecogs.com/svg.latex?M_\odot) for 'BULGE', 'THIN DISK', 'THICK DISK', 'DARK HALO', surface mass density ![](https://latex.codecogs.com/svg.latex?M_\odot/pc^2) for 'EXP DISK', and ![](https://latex.codecogs.com/svg.latex?M_\odot/kpc^3) for 'BURKERT HALO'.
* **checked** column is the boolean value for including that gravitational potential component.


Finally there is shown the graphic composition of rotation curve. This is an interactive user aided tool for include or exclude the available components of the gravitational potential, then graphically recover the dynamic mass composition for the observed rotation curve and ccomputed ![](https://latex.codecogs.com/svg.latex?\chi%5E2) value:

[![N|Solid](https://github.com/andresGranadosC/GalRotpy/blob/master/docs/GalRotpy_example.png?raw=true)](https://github.com/andresGranadosC/GalRotpy/blob/master/docs/GalRotpy_example.png)

Once you have setted the initial parameters by the graphical reconstruction, press Start. GalRotpy will allow you to enter some parameters for MCMC fiting of all the dimensions setted previously and to find the dark halo mass. Then, there is required to set the number of times that GalRotpy must iterate (default=1).

[![N|Solid](https://github.com/andresGranadosC/GalRotpy/blob/master/docs/terminal3.png?raw=true)](https://github.com/andresGranadosC/GalRotpy/blob/master/docs/Terminal.png)

When GalRotpy finishes the MCMC processes, it will show up the results in independent plots interactively as is shown below.

[![N|Solid](https://github.com/andresGranadosC/GalRotpy/blob/master/docs/Parameter_fit.png?raw=true)](https://github.com/andresGranadosC/GalRotpy/blob/master/docs/Parameter_fit.png)

[![N|Solid](https://github.com/andresGranadosC/GalRotpy/blob/master/docs/Parameter_fit_2.png?raw=true)](https://github.com/andresGranadosC/GalRotpy/blob/master/docs/Parameter_fit_2.png)

[![N|Solid](https://github.com/andresGranadosC/GalRotpy/blob/master/docs/Parameter_fit_3.png?raw=true)](https://github.com/andresGranadosC/GalRotpy/blob/master/docs/Parameter_fit_3.png)

Finally, GalRotpy print out all the results of the MCMC parameter fitting including ![](https://latex.codecogs.com/svg.latex?\chi%5E2) value,
 the estimations of dark matter halo considering the redshift and cosmological overdensity setted previously.

[![N|Solid](https://github.com/andresGranadosC/GalRotpy/blob/master/docs/terminal4.png?raw=true)](https://github.com/andresGranadosC/GalRotpy/blob/master/docs/Final_fit.png)

**Conf_Regions.pdf** and **GalRotpy_fit.pdf** files

![Model View Controller](https://github.com/andresGranadosC/GalRotpy/blob/master/docs/Conf_Regions.png?raw=true)

![Model View Controller](https://github.com/andresGranadosC/GalRotpy/blob/master/docs/GalRotpy_fit.png?raw=true)

All the results are compiled in the following file called **final_params.txt** whose content is shown below:

[![N|Solid](https://github.com/andresGranadosC/GalRotpy/blob/master/docs/final_plot.png?raw=true)](https://github.com/andresGranadosC/GalRotpy/blob/master/docs/final_plot.png)

GalRotpy produces the following output files:

 * **Conf_Regions.pdf**
 * **final_params.txt**
 * **GalRotpy_fit.pdf**

## GalRotpy is a powerful method to:

 * checking the presence of an assumed mass type component in a observed rotation curve,
 * determine quantitatively the main mass contribution in a galaxy by means of the mass ratios of a given set of five potentials,
 * to bound the contribution of each mass component given its gravitational potential parameters.

The related pre-print reference: [GalRotpy: an educational tool to understand and parametrize the rotation curve and gravitational potential of disk-like galaxies](https://arxiv.org/abs/1705.01665)

