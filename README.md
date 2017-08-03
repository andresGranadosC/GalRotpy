# GalRotpy

A Python-based tool for parametrizing the rotation curve and the galaxy potential of disk-like galaxies.


GalRotpy is intended to the understanding of the contributions of each mass component to the gravitational potential of disk-like galaxies and in the corresponding rotation curve. This software is presented as a pedagogical tool for supporting the following research areas: galactic dynamics, cosmology and gravitation.

GalRotpy requires essentially the [galpy](https://github.com/jobovy/galpy) package and the three fundamental Python packages: [Numpy](http://www.numpy.org/), [Scipy](https://www.scipy.org/) and [Matplotlib](http://matplotlib.org/). GalRotpy supports Python2.7.


## About the rotation curve and mass component parametrization

The GalRotpy tool can give a first approximation of the galaxy rotation curve using the following schemas:

 * bulge model uning a Miyamoto-Nagai potential,
 * stellar or gaseous disk:
   * thin or thick disks implementing Miyamoto-Nagai potentials, and/or
   * an exponential disk model
 * and finally the Dark Halo with a Navarro-Frenk-White (NFW) or with a Burkert potential.

The real time composition of rotation curves of disk-like galaxies is a simple and powerful method to:
 * checking the presence of an assumed mass type component in a observed rotation curve,
 * determine quantitatively the main mass contribution in a galaxy by means of the mass ratios of a given set of five potentials,
 * to bound the contribution of each mass component given its radial and height scales.

The related pre-print reference: [galpy](https://arxiv.org/abs/1705.01665)
