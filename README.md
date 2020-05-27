# GalRotpy

A Python3-based tool for parametrizing the rotation curve and the galaxy potential of disk-like galaxies.


GalRotpy allows to model the dynamical mass of disk-like galaxies. It makes a parametric fit of the rotation curve by means the composed gravitational potential of the galaxy. This software is presented as a pedagogical tool for supporting the following research areas: **galactic dynamics, cosmology and gravitation**.

GalRotpy requires essentially the [galpy](https://github.com/jobovy/galpy) package and the three fundamental Python3 packages: [Numpy](http://www.numpy.org/), [Scipy](https://www.scipy.org/) and [Matplotlib](http://matplotlib.org/). GalRotpy doesn't support Python2.


## About the rotation curve and mass component parametrization. Gravitational potentials available:

GalRotpy can give a first approximation of the galaxy rotation curve using the following gravitational potential schemas:

 * bulge model uning a Miyamoto-Nagai potential,
 * stellar or gaseous disk:
   * thin or thick disks implementing Miyamoto-Nagai potentials, and/or
   * an exponential disk model
 * and finally the Dark Halo with a Navarro-Frenk-White (NFW) or with a Burkert potential.

## How to use GalRotpy

Run GalRotpy with the following parameters:

```sh
$ python3 GalRotpy rot_curve.txt [bulge halo] | [bulge disk halo] | [disk halo]
```

 * **rot_curve.txt**  is the file that contains the rotation curve in three columns:
    * r (galactocentric distance in kpc)
    * vel (the speed in km/s at each distance r)
    * e_vel (the uncertainty for the circular velocity at each distance r)

|r| vel|    e_vel| 
| ------ | ------ | ------ |
|   0.24|  37.3|   6.2|
|   0.28|	37.9|	5.5|
|   0.46|	47.1|	2.8|
|   0.73|	55.1|	3.3|
|   ...|    ...|	...|

* **[bulge halo] | [bulge disk halo] | [disk halo]** optionally to set a first guess of the potential composition to reproduce the given rotation curve: bulge and halo or bulge and disk and halo or disk and halo
* Due to the last parameters are optional, if you don't enter any of the options above, GalRotpy will request you to introduce the following set of parameters:

|Component| Mass|   Parameter a|    Parameter b|
| ------ | ------ | ------ | ------ |
|'BULGE'|110000000.0|0.0|0.495|
|'THIN DISK'|3900000000.0|5.3|0.25|
|'THICK DISK'|39000000000.0|2.6|0.8|
|'EXP DISK'|500.0|5.3|------|
|'DARK HALO'|140000000000.0|13.0|------|
|'BURKERT HALO'|8000000.0|20.0|------|

Finally after the input parameters are checked, there is shown the graphic composition of rotation curves.



This is a simple and powerful method to:
 * checking the presence of an assumed mass type component in a observed rotation curve,
 * determine quantitatively the main mass contribution in a galaxy by means of the mass ratios of a given set of five potentials,
 * to bound the contribution of each mass component given its radial and height scales.

The related pre-print reference: [GalRotpy: an educational tool to understand and parametrize the rotation curve and gravitational potential of disk-like galaxies](https://arxiv.org/abs/1705.01665)
