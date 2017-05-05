# -*- coding: utf-8 -*-
"""
    GalRotpy.py - a Python-based tool for parametrizing galaxy potential by rotation curve

    Copyright (c) 2016 Andr\'es Granados
    All rights reserved.

    Permission is hereby granted, free of charge, to any person obtaining
    a copy of this software and associated documentation files (the
    "Software"), to deal in the Software without restriction, including
    without limitation the rights to use, copy, modify, merge, publish,
    distribute, sublicense, and/or sell copies of the Software, and to
    permit persons to whom the Software is furnished to do so, subject to
    the following conditions:

    The above copyright notice and this permission notice shall be
    included in all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
    MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
    IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY 
    CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
    TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
    SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

    Created on Fri Jan 6 07:00:00 MST 2017
"""

from matplotlib.widgets import Slider, Button, RadioButtons, CheckButtons # Matplotlib widgets
from matplotlib import pyplot as plt # Plotting interface
import numpy as np # Array managing
from galpy.potential import MiyamotoNagaiPotential, NFWPotential, RazorThinExponentialDiskPotential # GALPY potentials
from galpy.potential import calcRotcurve # composed rotation curve calculation for plotting
from astropy import units # Physical/real units data managing
from astropy import table as Table # For fast and easy reading / writing with tables using numpy library

input_params=Table.Table.read('input_params.txt', format='ascii.tab')

# Definition of initial parameters for Bulge potential (for units see the potential definition below)
a1=input_params['a (kpc)'][0]
b1=input_params['b (kpc)'][0]
amp1=input_params['mass'][0]

# Initial parameters for Thin disk potential
a2=input_params['a (kpc)'][1]
b2=input_params['b (kpc)'][1]
amp2=input_params['mass'][1]

# Initial parameters for Thick disk potential
a3=input_params['a (kpc)'][2]
b3=input_params['b (kpc)'][2]
amp3=input_params['mass'][2]

# Initial parameters for Halo potential
a5=input_params['a (kpc)'][4]
amp5=input_params['mass'][4]

# Initial parameters for exponential disk potential
h_r=input_params['a (kpc)'][3]
amp4=input_params['mass'][3]

x_offset = 0 #0.35 # It defines a radial coordinate offset as user input
r_0=1*units.kpc
v_0=220*units.km/units.s

name='rot_curve.txt' # Set the name of the galaxy's txt file rotation curve
tt=Table.Table.read(name, format='ascii.tab')
r_data=tt['r']-x_offset # The txt file must contain the radial coordinate values in kpc
v_c_data=tt['vel'] # velocity in km/s
v_c_err_data = tt['e_vel'] # and velocity error in km/s


lista=np.linspace(0.00001, np.max(r_data), 10*len(r_data)) # radial coordinate for the rotation curve calculation
# Potentials definition using physical units (amplitude in Solar masses, scales in kpc and surface density in Solar masses / pc^2 )
MN_Bulge_p= MiyamotoNagaiPotential(amp=amp1*units.Msun,a=a1*units.kpc,b=b1*units.kpc,normalize=False,ro=r_0, vo=v_0)
MN_Thin_Disk_p= MiyamotoNagaiPotential(amp=amp2*units.Msun,a=a2*units.kpc,b=b2*units.kpc,normalize=False,ro=r_0, vo=v_0)
MN_Thick_Disk_p= MiyamotoNagaiPotential(amp=amp3*units.Msun,a=a3*units.kpc,b=b3*units.kpc,normalize=False,ro=r_0, vo=v_0)
EX_Disk_p = RazorThinExponentialDiskPotential(amp=amp4*(units.Msun/(units.pc**2)), hr=h_r*units.kpc, maxiter=20, tol=0.001, normalize=False, ro=r_0, vo=v_0, new=True, glorder=100)
NFW_p = NFWPotential(amp=amp5*units.Msun, a=a5*units.kpc, normalize=False, ro=r_0, vo=v_0)
# Circular velocities in km/s
MN_Bulge = calcRotcurve(MN_Bulge_p, lista, phi=None)*220
MN_Thin_Disk = calcRotcurve(MN_Thin_Disk_p, lista, phi=None)*220
MN_Thick_Disk = calcRotcurve(MN_Thick_Disk_p, lista, phi=None)*220
EX_Disk = calcRotcurve(EX_Disk_p, lista, phi=None)*220
NFW = calcRotcurve(NFW_p, lista, phi=None)*220
# Circular velocity for the composition of 5 potentials in km/s
v_circ_comp = calcRotcurve([MN_Bulge_p,MN_Thin_Disk_p,MN_Thick_Disk_p, EX_Disk_p, NFW_p], lista, phi=None)*220

# Global variables for the plot
fig = plt.figure()
ax = fig.add_axes((0.35, 0.1, 0.6, 0.85))
ax.set_xlim([0, r_data[-1]])
# A plot for each rotation curve with the colors indicated below
MN_b_plot, = ax.plot(lista, MN_Bulge, linestyle='--', c='gray')
MN_td_plot, = ax.plot(lista, MN_Thin_Disk, linestyle='--', c='red')
MN_tkd_plot, = ax.plot(lista, MN_Thick_Disk, linestyle='--', c='blue')
EX_d_plot, = ax.plot(lista, EX_Disk, linestyle='--', c='cyan')
NFW_plot, = ax.plot(lista, NFW, linestyle='--', c='green')
CV_galaxy = ax.errorbar(r_data, v_c_data, v_c_err_data,  c='k', fmt='', ls='none')
CV_galaxy_dot = ax.scatter(r_data, v_c_data, c='k')
# Composed rotation curve
v_circ_comp_plot, = ax.plot(lista, v_circ_comp, c='k')
ax.cla()

#fig = plt.figure()
#ax = fig.add_axes((0.35, 0.1, 0.6, 0.85))
ax.set_title(name, fontsize=30)
ax.set_xlabel(r'$R$ (kpc)', fontsize=20)
ax.set_ylabel(r'$v_c$ (km/s)', fontsize=20)
ax.set_xlim([0, r_data[-1]])
ax.set_ylim([0,np.max(v_c_data)*1.1])
MN_Bulge_p= MiyamotoNagaiPotential(amp=amp1*units.Msun,a=a1*units.kpc,b=b1*units.kpc,normalize=False,ro=r_0, vo=v_0)
MN_Thin_Disk_p= MiyamotoNagaiPotential(amp=amp2*units.Msun,a=a2*units.kpc,b=b2*units.kpc,normalize=False,ro=r_0, vo=v_0)
MN_Thick_Disk_p= MiyamotoNagaiPotential(amp=amp3*units.Msun,a=a3*units.kpc,b=b3*units.kpc,normalize=False,ro=r_0, vo=v_0)
EX_Disk_p = RazorThinExponentialDiskPotential(amp=amp4*(units.Msun/(units.pc**2)), hr=h_r*units.kpc, maxiter=20, tol=0.001, normalize=False, ro=r_0, vo=v_0, new=True, glorder=100)
NFW_p = NFWPotential(amp=amp5*units.Msun, a=a5*units.kpc, normalize=False)

MN_Bulge = calcRotcurve(MN_Bulge_p, lista, phi=None)*220
MN_Thin_Disk = calcRotcurve(MN_Thin_Disk_p, lista, phi=None)*220
MN_Thick_Disk = calcRotcurve(MN_Thick_Disk_p, lista, phi=None)*220
EX_Disk = calcRotcurve(EX_Disk_p, lista, phi=None)*220
NFW = calcRotcurve(NFW_p, lista, phi=None)*220
v_circ_comp = calcRotcurve([MN_Bulge_p,MN_Thin_Disk_p,MN_Thick_Disk_p, NFW_p], lista, phi=None)*220

MN_b_plot, = ax.plot(lista, MN_Bulge, linestyle='--', c='gray')
MN_td_plot, = ax.plot(lista, MN_Thin_Disk, linestyle='--', c='red')
MN_tkd_plot, = ax.plot(lista, MN_Thick_Disk, linestyle='--', c='blue')
EX_d_plot, = ax.plot(lista, EX_Disk, linestyle='--', c='cyan')
NFW_plot, = ax.plot(lista, NFW, linestyle='--', c='green')
CV_galaxy = ax.errorbar(r_data, v_c_data, v_c_err_data,  c='k', fmt='', ls='none')
CV_galaxy_dot = ax.scatter(r_data, v_c_data, c='k')
v_circ_comp_plot, = ax.plot(lista, v_circ_comp, c='k')

# Checkbox for selecting the potentials to compose the rotation
rax = plt.axes((0.05,0.8,0.2,0.15))
check = CheckButtons(rax, ('MN Bulge (GRAY)', 'MN Thin Disk (RED)', 'MN Thick Disk (BLUE)', 'Exp. Disk (CYAN)', 'NFW - Halo (GREEN)'), (True, True, True, True, True))
# Sliders for Bulge - gray
MN_b_amp_ax = fig.add_axes((0.05,0.75,0.17,0.03))
MN_b_amp_s = Slider(MN_b_amp_ax, r"M($M_\odot$)", input_params['mass'][0]/(10**input_params['threshold_mass'][0]), input_params['mass'][0]*(10**input_params['threshold_mass'][0]), valinit=input_params['mass'][0], color='gray', valfmt='%1.3E')
MN_b_a_ax = fig.add_axes((0.05,0.72,0.17,0.03))
MN_b_a_s = Slider(MN_b_a_ax, "a (kpc)", 0, 0.01*input_params['threshold_a'][0], valinit=0.0, color='gray')
MN_b_b_ax = fig.add_axes((0.05,0.69,0.17,0.03))
MN_b_b_s = Slider(MN_b_b_ax, "b (kpc)", input_params['b (kpc)'][0]*(1-0.01*input_params['threshold_b'][0]), input_params['b (kpc)'][0]*(1+0.01*input_params['threshold_b'][0]), valinit=input_params['b (kpc)'][0], color='gray')
# Sliders for Thin disk - red
MN_td_amp_ax = fig.add_axes((0.05,0.63,0.17,0.03))
MN_td_amp_s = Slider(MN_td_amp_ax, r"M($M_\odot$)", input_params['mass'][1]/(10**input_params['threshold_mass'][1]), input_params['mass'][1]*(10**input_params['threshold_mass'][1]), valinit=input_params['mass'][1], color='red', valfmt='%1.3E')
MN_td_a_ax = fig.add_axes((0.05,0.6,0.17,0.03))
MN_td_a_s = Slider(MN_td_a_ax, "a (kpc)", input_params['a (kpc)'][1]*(1-0.01*input_params['threshold_a'][1]), input_params['a (kpc)'][1]*(1+0.01*input_params['threshold_a'][1]), valinit=input_params['a (kpc)'][1], color='red')
MN_td_b_ax = fig.add_axes((0.05,0.57,0.17,0.03))
MN_td_b_s = Slider(MN_td_b_ax, "b (kpc)", input_params['b (kpc)'][1]*(1-0.01*input_params['threshold_b'][1]), input_params['b (kpc)'][1]*(1+0.01*input_params['threshold_b'][1]), valinit=input_params['b (kpc)'][1], color='red')
# Sliders for Thick disk - Blue
MN_tkd_amp_ax = fig.add_axes((0.05,0.51,0.17,0.03))
MN_tkd_amp_s = Slider(MN_tkd_amp_ax, r"M($M_\odot$)", input_params['mass'][2]/(10**input_params['threshold_mass'][2]), input_params['mass'][2]*(10**input_params['threshold_mass'][2]), valinit=input_params['mass'][2], color='blue', valfmt='%1.3E')
MN_tkd_a_ax = fig.add_axes((0.05,0.48,0.17,0.03))
MN_tkd_a_s = Slider(MN_tkd_a_ax, "a (kpc)", input_params['a (kpc)'][2]*(1-0.01*input_params['threshold_a'][2]), input_params['a (kpc)'][2]*(1+0.01*input_params['threshold_a'][2]), valinit=input_params['a (kpc)'][2], color='blue')
MN_tkd_b_ax = fig.add_axes((0.05,0.45,0.17,0.03))
MN_tkd_b_s = Slider(MN_tkd_b_ax, "b (kpc)", input_params['b (kpc)'][2]*(1-0.01*input_params['threshold_b'][2]), input_params['b (kpc)'][2]*(1+0.01*input_params['threshold_b'][2]), valinit=input_params['b (kpc)'][2], color='blue')

#MN_tkd_amp_ax.set_visible=FalseMN_tkd_b_s
# Sliders for exponential disk - Cyan
MN_ed_amp_ax = fig.add_axes((0.07,0.39,0.17,0.03))
MN_ed_amp_s = Slider(MN_ed_amp_ax, r"M($M_\odot/pc^2$)", input_params['mass'][3]/(10**input_params['threshold_mass'][3]), input_params['mass'][3]*(10**input_params['threshold_mass'][3]), valinit=input_params['mass'][3], color='cyan', valfmt='%1.3E')
MN_ed_a_ax = fig.add_axes((0.05,0.36,0.17,0.03))
MN_ed_a_s = Slider(MN_ed_a_ax, "h_r (kpc)", input_params['a (kpc)'][3]*(1-0.01*input_params['threshold_a'][3]), input_params['a (kpc)'][3]*(1+0.01*input_params['threshold_a'][3]), valinit=input_params['a (kpc)'][3], color='cyan')
# Sliders for Halo - green
NFW_amp_ax = fig.add_axes((0.05,0.27,0.17,0.03))
NFW_amp_s = Slider(NFW_amp_ax, r"M($M_\odot$)", input_params['mass'][4]/(10*input_params['threshold_mass'][4]), input_params['mass'][4]*(10**input_params['threshold_mass'][4]), valinit=input_params['mass'][4], color='green', valfmt='%1.3E')
NFW_a_ax = fig.add_axes((0.05,0.24,0.17,0.03))
NFW_a_s = Slider(NFW_a_ax, "a (kpc)", input_params['a (kpc)'][4]*(1-0.01*input_params['threshold_a'][4]), input_params['a (kpc)'][4]*(1+0.01*input_params['threshold_a'][4]), valinit=input_params['a (kpc)'][4], color='green')
# Function for setting new parameters for each potential
def MN_b_amp_s_func(val):
    if MN_b_plot.get_visible() == True:
        global MN_Bulge_p, amp1, a1, b1
        amp1=val*1
        #MN_b_amp_s.valtext.set_text(10**val)
        MN_Bulge_p = MiyamotoNagaiPotential(amp=val*units.Msun,a=a1*units.kpc,b=b1*units.kpc,normalize=False,ro=r_0, vo=v_0)   
        update_rot_curve()
def MN_b_a_s_func(val):
    if MN_b_plot.get_visible() == True:
        global MN_Bulge_p, amp1, a1, b1
        a1=val*1
        MN_Bulge_p = MiyamotoNagaiPotential(amp=amp1*units.Msun,a=val*units.kpc,b=b1*units.kpc,normalize=False,ro=r_0, vo=v_0)   
        update_rot_curve()
def MN_b_b_s_func(val):
    if MN_b_plot.get_visible() == True:
        global MN_Bulge_p, amp1, a1, b1
        b1=val*1
        MN_Bulge_p = MiyamotoNagaiPotential(amp=amp1*units.Msun,a=a1*units.kpc,b=val*units.kpc,normalize=False,ro=r_0, vo=v_0)   
        update_rot_curve()
def MN_td_amp_s_func(val):
    if MN_td_plot.get_visible() == True:
        global MN_Thin_Disk_p, amp2, a2, b2
        amp2=val*1
        MN_Thin_Disk_p= MiyamotoNagaiPotential(amp=val*units.Msun,a=a2*units.kpc,b=b2*units.kpc,normalize=False,ro=r_0, vo=v_0)  
        update_rot_curve()
def MN_td_a_s_func(val):
    if MN_td_plot.get_visible() == True:
        global MN_Thin_Disk_p, amp2, a2, b2
        a2=val*1
        MN_Thin_Disk_p= MiyamotoNagaiPotential(amp=amp2*units.Msun,a=val*units.kpc,b=b2*units.kpc,normalize=False,ro=r_0, vo=v_0)  
        update_rot_curve()
def MN_td_b_s_func(val):
    if MN_td_plot.get_visible() == True:
        global MN_Thin_Disk_p, amp2, a2, b2
        b2=val*1
        MN_Thin_Disk_p= MiyamotoNagaiPotential(amp=amp2*units.Msun,a=a2*units.kpc,b=val*units.kpc,normalize=False,ro=r_0, vo=v_0)  
        update_rot_curve()
def MN_tkd_amp_s_func(val):
    if MN_tkd_plot.get_visible() == True:
        global MN_Thick_Disk_p, amp3, a3, b3
        amp3=val*1
        MN_Thick_Disk_p= MiyamotoNagaiPotential(amp=val*units.Msun,a=a3*units.kpc,b=b3*units.kpc,normalize=False,ro=r_0, vo=v_0)  
        update_rot_curve()
def MN_tkd_a_s_func(val):
    if MN_tkd_plot.get_visible() == True:
        global MN_Thick_Disk_p, amp3, a3, b3
        a3=val*1
        MN_Thick_Disk_p= MiyamotoNagaiPotential(amp=amp3*units.Msun,a=val*units.kpc,b=b3*units.kpc,normalize=False,ro=r_0, vo=v_0)  
        update_rot_curve()
def MN_tkd_b_s_func(val):
    if MN_tkd_plot.get_visible() == True:
        global MN_Thick_Disk_p, amp3, a3, b3
        b3=val*1
        MN_Thick_Disk_p= MiyamotoNagaiPotential(amp=amp3*units.Msun,a=a3*units.kpc,b=val*units.kpc,normalize=False,ro=r_0, vo=v_0)  
        update_rot_curve()
def NFW_amp_s_func(val):
    if NFW_plot.get_visible() == True:
        global NFW_p, amp5,a5
        amp5=val*1
        NFW_p = NFWPotential(amp=val*units.Msun, a=a5*units.kpc, normalize=False, ro=r_0, vo=v_0)
        update_rot_curve()    
def NFW_a_s_func(val):
    if NFW_plot.get_visible() == True:
        global NFW_p, amp5,a5
        a5=val*1
        NFW_p = NFWPotential(amp=amp5*units.Msun, a=val*units.kpc, normalize=False, ro=r_0, vo=v_0)
        update_rot_curve()
def MN_ed_amp_s_func(val):
    if EX_d_plot.get_visible() == True:
        global EX_Disk_p, amp4,h_r
        amp4=val*1
        EX_Disk_p = RazorThinExponentialDiskPotential(amp=val*(units.Msun/(units.pc**2)), hr=h_r*units.kpc, maxiter=20, tol=0.001, normalize=False, ro=r_0, vo=v_0, new=True, glorder=100)
        update_rot_curve() 
def MN_ed_a_s_func(val):
    if EX_d_plot.get_visible() == True:
        global EX_Disk_p, amp4,h_r
        h_r=val*1
        EX_Disk_p = RazorThinExponentialDiskPotential(amp=amp4*(units.Msun/(units.pc**2)), hr=val*units.kpc, maxiter=20, tol=0.001, normalize=False, ro=r_0, vo=v_0, new=True, glorder=100)
        update_rot_curve()
    
# update rotation curve for the selected and the composed potential
def update_rot_curve():
    ax.cla()
    global MN_b_plot, MN_Bulge_p, MN_Thin_Disk_p,MN_Thick_Disk_p, MN_td_plot,MN_tkd_plot, NFW_p, NFW_plot, EX_d_plot, EX_Disk_p, CV_galaxy, CV_galaxy_dot
    composite_pot_array=[]
    ax.set_xlabel(r'$R$ (kpc)', fontsize=20)
    ax.set_ylabel(r'$v_c$ (km/s)', fontsize=20)
    ax.set_xlim([0, r_data[-1]])
    ax.set_ylim([0,np.max(v_c_data)*1.1])
    if MN_b_plot.get_visible() == True:
        MN_Bulge = calcRotcurve(MN_Bulge_p, lista, phi=None)*220
        MN_b_plot, = ax.plot(lista, MN_Bulge, linestyle='--', c='gray')
        composite_pot_array.append(MN_Bulge_p)
    if MN_td_plot.get_visible() == True:
        MN_Thin_Disk = calcRotcurve(MN_Thin_Disk_p, lista, phi=None)*220
        MN_td_plot, = ax.plot(lista, MN_Thin_Disk, linestyle='--', c='red')
        composite_pot_array.append(MN_Thin_Disk_p)
    if MN_tkd_plot.get_visible() == True:
        MN_Thick_Disk = calcRotcurve(MN_Thick_Disk_p, lista, phi=None)*220
        MN_tkd_plot, = ax.plot(lista, MN_Thick_Disk, linestyle='--', c='blue')
        composite_pot_array.append(MN_Thick_Disk_p)
    if NFW_plot.get_visible() == True:
        NFW = calcRotcurve(NFW_p, lista, phi=None)*220
        NFW_plot, = ax.plot(lista, NFW, linestyle='--', c='green')
        composite_pot_array.append(NFW_p)
    if EX_d_plot.get_visible() == True:
        EX_Disk = calcRotcurve(EX_Disk_p, lista, phi=None)*220
        EX_d_plot, = ax.plot(lista, EX_Disk, linestyle='--', c='cyan')
        composite_pot_array.append(EX_Disk_p)
    CV_galaxy = ax.errorbar(r_data, v_c_data, v_c_err_data,  c='k', fmt='', ls='none')
    CV_galaxy_dot = ax.scatter(r_data, v_c_data, c='k')
    v_circ_comp = calcRotcurve(composite_pot_array, lista, phi=None)*220
    v_circ_comp_plot, = ax.plot(lista, v_circ_comp, c='k')

MN_Bulge = calcRotcurve(MN_Bulge_p, lista, phi=None)*220
MN_Thin_Disk = calcRotcurve(MN_Thin_Disk_p, lista, phi=None)*220
MN_Thick_Disk = calcRotcurve(MN_Thick_Disk_p, lista, phi=None)*220
EX_Disk = calcRotcurve(EX_Disk_p, lista, phi=None)*220
NFW = calcRotcurve(NFW_p, lista, phi=None)*220
v_circ_comp = calcRotcurve([MN_Bulge_p,MN_Thin_Disk_p,MN_Thick_Disk_p, NFW_p], lista, phi=None)*220

ax.set_xlabel(r'$R$ (kpc)', fontsize=20)
ax.set_ylabel(r'$v_c$ (km/s)', fontsize=20)
ax.set_xlim([0, r_data[-1]])
ax.set_ylim([0,np.max(v_c_data)*1.1])
MN_b_plot, = ax.plot(lista, MN_Bulge, linestyle='--', c='gray')
MN_td_plot, = ax.plot(lista, MN_Thin_Disk, linestyle='--', c='red')
MN_tkd_plot, = ax.plot(lista, MN_Thick_Disk, linestyle='--', c='blue')
EX_d_plot, = ax.plot(lista, EX_Disk, linestyle='--', c='cyan')
NFW_plot, = ax.plot(lista, NFW, linestyle='--', c='green')
CV_galaxy = ax.errorbar(r_data, v_c_data, v_c_err_data,  c='k', fmt='', ls='none')
CV_galaxy_dot = ax.scatter(r_data, v_c_data, c='k')
v_circ_comp_plot, = ax.plot(lista, v_circ_comp, c='k')


# Define the sliders update functions
MN_b_amp_s.on_changed(MN_b_amp_s_func)
MN_b_a_s.on_changed(MN_b_a_s_func)
MN_b_b_s.on_changed(MN_b_b_s_func)
MN_td_amp_s.on_changed(MN_td_amp_s_func)
MN_td_a_s.on_changed(MN_td_a_s_func)
MN_td_b_s.on_changed(MN_td_b_s_func)
MN_tkd_amp_s.on_changed(MN_tkd_amp_s_func)
MN_tkd_a_s.on_changed(MN_tkd_a_s_func)
MN_tkd_b_s.on_changed(MN_tkd_b_s_func)
NFW_amp_s.on_changed(NFW_amp_s_func)
NFW_a_s.on_changed(NFW_a_s_func)
MN_ed_amp_s.on_changed(MN_ed_amp_s_func)
MN_ed_a_s.on_changed(MN_ed_a_s_func)
# Enable/disable the selected potential for the composed rotation curve
def check_on_clicked(label):
    if label == 'MN Bulge (GRAY)':
        MN_b_plot.set_visible(not MN_b_plot.get_visible())
        update_rot_curve()
    elif label == 'MN Thin Disk (RED)':
        MN_td_plot.set_visible(not MN_td_plot.get_visible())
        update_rot_curve()
    elif label == 'MN Thick Disk (BLUE)':
        MN_tkd_plot.set_visible(not MN_tkd_plot.get_visible())
        update_rot_curve()
    elif label == 'Exp. Disk (CYAN)':
        EX_d_plot.set_visible(not EX_d_plot.get_visible())
        update_rot_curve()
    elif label == 'NFW - Halo (GREEN)':
        NFW_plot.set_visible(not NFW_plot.get_visible())
        update_rot_curve()
    plt.draw()
    
# Plotting all the curves
ax.set_xlabel(r'$R$ (kpc)', fontsize=20)
ax.set_ylabel(r'$v_c$ (km/s)', fontsize=20)
ax.set_xlim([0, r_data[-1]])
ax.set_ylim([0,np.max(v_c_data)*1.1])
check.on_clicked(check_on_clicked)
plt.ylim([0,np.max(v_c_data)*1.1])
plt.show()

chk=[]
if MN_b_plot.get_visible() == True:
    chk.append(True)
else:
    chk.append(False)
if MN_td_plot.get_visible() == True:
    chk.append(True)
else:
    chk.append(False)
if MN_tkd_plot.get_visible() == True:
    chk.append(True)
else:
    chk.append(False)
if EX_d_plot.get_visible() == True:
    chk.append(True)
    print "La Masa total del modelo exponencial es:", 2*np.pi*(h_r**2.)*amp4*1.e6, " M_sun"
else:
    chk.append(False)
if NFW_plot.get_visible() == True:
    chk.append(True)
    print "La Masa virial del modelo NFW es: ", '%E' % NFW_p.mvir(), "M_sun"
    print "El parámetro de concentracion del modelo NFW es: ",'%E' % NFW_p.conc()
else:
    chk.append(False)

compnts = ['BULGE','THIN DISK','THICK DISK','EXP DISK', 'DARK HALO']
masses = [amp1, amp2, amp3, amp4, amp5]
aa = [a1, a2, a3, h_r, a5]
bb = [b1, b2, b3, 0, 0]
#chk = [True,True,False,False,False]


output_parameters = Table.Table([compnts,masses, aa,bb, chk], names=('component', 'mass', 'a (kpc)', 'b (kpc)', 'checked'), meta={'name': 'output_parameters'})

output_parameters.write('output_parameters.txt', format='ascii.tab', overwrite=True)
