"""
    RotPy.py - a Python-based tool for parametrizing galaxy potential by rotation curve

    Copyright (c) 2016 Andrés Granados
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

from matplotlib.widgets import Slider, Button, RadioButtons, CheckButtons
from matplotlib import pyplot as plt
import numpy as np
from galpy.potential import MiyamotoNagaiPotential, NFWPotential, RazorThinExponentialDiskPotential
from galpy.potential import calcRotcurve
import numpy as np
from astropy import units
from astropy import table as Table


a1=0.
b1=0.3
amp1=460*2.32e7

a2=5.3
b2=0.25
amp2=1700*2.32e7

a3=2.6
b3=0.8
amp3=1700*2.32e7

a5=14
amp5=6000*2.32e7

h_r=3.3
amp4=1e3

name='NGC6361'
tt=Table.Table.read('rot_curve_'+name+'.txt', format='ascii.tab')
r_califa=tt['r']
v_c_califa=tt['vel']
v_c_err_califa = tt['e_vel']


lista=np.linspace(0.00001, np.max(r_califa), 10*len(r_califa))

MN_Bulge_p= MiyamotoNagaiPotential(amp=amp1*units.Msun,a=a1*units.kpc,b=b1*units.kpc,normalize=False,ro=8*units.kpc, vo=220*units.km/units.s)
MN_Thin_Disk_p= MiyamotoNagaiPotential(amp=amp2*units.Msun,a=a2*units.kpc,b=b2*units.kpc,normalize=False,ro=8*units.kpc, vo=220*units.km/units.s)
MN_Thick_Disk_p= MiyamotoNagaiPotential(amp=amp3*units.Msun,a=a3*units.kpc,b=b3*units.kpc,normalize=False,ro=8*units.kpc, vo=220*units.km/units.s)
EX_Disk_p = RazorThinExponentialDiskPotential(amp=amp4*(units.Msun/(units.pc**2)), hr=h_r*units.kpc, maxiter=20, tol=0.001, normalize=False, ro=8*units.kpc, vo=220*units.km/units.s, new=True, glorder=100)
NFW_p = NFWPotential(amp=amp5*units.Msun, a=a5*units.kpc, normalize=False)

MN_Bulge = calcRotcurve(MN_Bulge_p, lista, phi=None)*220
MN_Thin_Disk = calcRotcurve(MN_Thin_Disk_p, lista, phi=None)*220
MN_Thick_Disk = calcRotcurve(MN_Thick_Disk_p, lista, phi=None)*220
EX_Disk = calcRotcurve(EX_Disk_p, lista, phi=None)*220
NFW = calcRotcurve(NFW_p, lista, phi=None)*220
v_circ_comp = calcRotcurve([MN_Bulge_p,MN_Thin_Disk_p,MN_Thick_Disk_p, EX_Disk_p, NFW_p], lista, phi=None)*220


fig = plt.figure()
ax = fig.add_axes((0.35, 0.05, 0.6, 0.85))

MN_b_plot, = ax.plot(lista, MN_Bulge, linestyle='--', c='gray')
MN_td_plot, = ax.plot(lista, MN_Thin_Disk, linestyle='--', c='red')
MN_tkd_plot, = ax.plot(lista, MN_Thick_Disk, linestyle='--', c='blue')
EX_d_plot, = ax.plot(lista, EX_Disk, linestyle='--', c='cyan')
NFW_plot, = ax.plot(lista, NFW, linestyle='--', c='green')
v_circ_comp_plot, = ax.plot(lista, v_circ_comp, c='k')




fig = plt.figure()
ax = fig.add_axes((0.35, 0.1, 0.6, 0.85))
ax.set_title(name, fontsize=30)
ax.set_xlabel(r'$R$ (kpc)', fontsize=20)
ax.set_ylabel(r'$v_c$ (km/s)', fontsize=20)
MN_Bulge_p= MiyamotoNagaiPotential(amp=amp1*units.Msun,a=a1*units.kpc,b=b1*units.kpc,normalize=False,ro=8*units.kpc, vo=220*units.km/units.s)
MN_Thin_Disk_p= MiyamotoNagaiPotential(amp=amp2*units.Msun,a=a2*units.kpc,b=b2*units.kpc,normalize=False,ro=8*units.kpc, vo=220*units.km/units.s)
MN_Thick_Disk_p= MiyamotoNagaiPotential(amp=amp3*units.Msun,a=a3*units.kpc,b=b3*units.kpc,normalize=False,ro=8*units.kpc, vo=220*units.km/units.s)
EX_Disk_p = RazorThinExponentialDiskPotential(amp=amp4*(units.Msun/(units.pc**2)), hr=h_r*units.kpc, maxiter=20, tol=0.001, normalize=False, ro=8*units.kpc, vo=220*units.km/units.s, new=True, glorder=100)
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
v_circ_comp_plot, = ax.plot(lista, v_circ_comp, c='k')

rax = plt.axes((0.05,0.8,0.2,0.15))
check = CheckButtons(rax, ('MN Bulge (GRAY)', 'MN Thin Disk (RED)', 'MN Thick Disk (BLUE)', 'Exp. Disk (CYAN)', 'NFW - Halo (GREEN)'), (True, True, True, True, True))

MN_b_amp_ax = fig.add_axes((0.05,0.75,0.17,0.03))
MN_b_amp_s = Slider(MN_b_amp_ax, r"M($M_\odot$)", 1e8, 1e11, valinit=460*2.32e7, color='gray', valfmt='%1.3E')
MN_b_a_ax = fig.add_axes((0.05,0.72,0.17,0.03))
MN_b_a_s = Slider(MN_b_a_ax, "a (kpc)", 0, 0.1, valinit=0.01, color='gray')
MN_b_b_ax = fig.add_axes((0.05,0.69,0.17,0.03))
MN_b_b_s = Slider(MN_b_b_ax, "b (kpc)", 0., 1., valinit=0.3, color='gray')

MN_td_amp_ax = fig.add_axes((0.05,0.63,0.17,0.03))
MN_td_amp_s = Slider(MN_td_amp_ax, r"M($M_\odot$)", 1e9, 1e11, valinit=1700*2.32e7, color='red', valfmt='%1.3E')
MN_td_a_ax = fig.add_axes((0.05,0.6,0.17,0.03))
MN_td_a_s = Slider(MN_td_a_ax, "a (kpc)", 1., 10., valinit=5.3, color='red')
MN_td_b_ax = fig.add_axes((0.05,0.57,0.17,0.03))
MN_td_b_s = Slider(MN_td_b_ax, "b (kpc)", 0.01, 1., valinit=0.25, color='red')

MN_tkd_amp_ax = fig.add_axes((0.05,0.51,0.17,0.03))
MN_tkd_amp_s = Slider(MN_tkd_amp_ax, r"M($M_\odot$)", 1e9, 1e11, valinit=1700*2.32e7, color='blue', valfmt='%1.3E')
MN_tkd_a_ax = fig.add_axes((0.05,0.48,0.17,0.03))
MN_tkd_a_s = Slider(MN_tkd_a_ax, "a (kpc)", 1., 10., valinit=2.6, color='blue')
MN_tkd_b_ax = fig.add_axes((0.05,0.45,0.17,0.03))
MN_tkd_b_s = Slider(MN_tkd_b_ax, "b (kpc)", 0.1, 15, valinit=8, color='blue')

#MN_tkd_amp_ax.set_visible=FalseMN_tkd_b_s

MN_ed_amp_ax = fig.add_axes((0.05,0.39,0.17,0.03))
MN_ed_amp_s = Slider(MN_ed_amp_ax, r"M($M_\odot/pc^2$)", 1e2, 1.3e3, valinit=1e3, color='cyan', valfmt='%1.3E')
MN_ed_a_ax = fig.add_axes((0.05,0.36,0.17,0.03))
MN_ed_a_s = Slider(MN_ed_a_ax, "h_r (kpc)", 2, 6, valinit=3.3, color='cyan')

NFW_amp_ax = fig.add_axes((0.05,0.27,0.17,0.03))
NFW_amp_s = Slider(NFW_amp_ax, r"M($M_\odot$)", 1e10, 1e12, valinit=6000*2.32e7, color='green', valfmt='%1.3E')
NFW_a_ax = fig.add_axes((0.05,0.24,0.17,0.03))
NFW_a_s = Slider(NFW_a_ax, "a (kpc)", 8, 30, valinit=14, color='green')

def MN_b_amp_s_func(val):
    if MN_b_plot.get_visible() == True:
        global MN_Bulge_p, amp1, a1, b1
        amp1=val*1
        MN_Bulge_p = MiyamotoNagaiPotential(amp=val*units.Msun,a=a1*units.kpc,b=b1*units.kpc,normalize=False,ro=8*units.kpc, vo=220*units.km/units.s)   
        update_rot_curve()
def MN_b_a_s_func(val):
    global MN_Bulge_p, amp1, a1, b1
    a1=val*1
    MN_Bulge_p = MiyamotoNagaiPotential(amp=amp1*units.Msun,a=val*units.kpc,b=b1*units.kpc,normalize=False,ro=8*units.kpc, vo=220*units.km/units.s)   
    update_rot_curve()
def MN_b_b_s_func(val):
    global MN_Bulge_p, amp1, a1, b1
    b1=val*1
    MN_Bulge_p = MiyamotoNagaiPotential(amp=amp1*units.Msun,a=a1*units.kpc,b=val*units.kpc,normalize=False,ro=8*units.kpc, vo=220*units.km/units.s)   
    update_rot_curve()
def MN_td_amp_s_func(val):
    global MN_Thin_Disk_p, amp2, a2, b2
    amp2=val*1
    MN_Thin_Disk_p= MiyamotoNagaiPotential(amp=val*units.Msun,a=a2*units.kpc,b=b2*units.kpc,normalize=False,ro=8*units.kpc, vo=220*units.km/units.s)  
    update_rot_curve()
def MN_td_a_s_func(val):
    global MN_Thin_Disk_p, amp2, a2, b2
    a2=val*1
    MN_Thin_Disk_p= MiyamotoNagaiPotential(amp=amp2*units.Msun,a=val*units.kpc,b=b2*units.kpc,normalize=False,ro=8*units.kpc, vo=220*units.km/units.s)  
    update_rot_curve()
def MN_td_b_s_func(val):
    global MN_Thin_Disk_p, amp2, a2, b2
    b2=val*1
    MN_Thin_Disk_p= MiyamotoNagaiPotential(amp=amp2*units.Msun,a=a2*units.kpc,b=val*units.kpc,normalize=False,ro=8*units.kpc, vo=220*units.km/units.s)  
    update_rot_curve()
def MN_tkd_amp_s_func(val):
    global MN_Thick_Disk_p, amp3, a3, b3
    amp3=val*1
    MN_Thick_Disk_p= MiyamotoNagaiPotential(amp=val*units.Msun,a=a3*units.kpc,b=b3*units.kpc,normalize=False,ro=8*units.kpc, vo=220*units.km/units.s)  
    update_rot_curve()
def MN_tkd_a_s_func(val):
    global MN_Thick_Disk_p, amp3, a3, b3
    a3=val*1
    MN_Thick_Disk_p= MiyamotoNagaiPotential(amp=amp3*units.Msun,a=val*units.kpc,b=b3*units.kpc,normalize=False,ro=8*units.kpc, vo=220*units.km/units.s)  
    update_rot_curve()
def MN_tkd_b_s_func(val):
    global MN_Thick_Disk_p, amp3, a3, b3
    b3=val*1
    MN_Thick_Disk_p= MiyamotoNagaiPotential(amp=amp3*units.Msun,a=a3*units.kpc,b=val*units.kpc,normalize=False,ro=8*units.kpc, vo=220*units.km/units.s)  
    update_rot_curve()
def NFW_amp_s_func(val):
    global NFW_p, amp5,a5
    amp5=val*1
    NFW_p = NFWPotential(amp=val*units.Msun, a=a5*units.kpc, normalize=False)
    update_rot_curve()    
def NFW_a_s_func(val):
    global NFW_p, amp5,a5
    a5=val*1
    NFW_p = NFWPotential(amp=amp5*units.Msun, a=val*units.kpc, normalize=False)
    update_rot_curve()
def MN_ed_amp_s_func(val):
    global EX_Disk_p, amp4,h_r
    amp4=val*1
    EX_Disk_p = RazorThinExponentialDiskPotential(amp=val*(units.Msun/(units.pc**2)), hr=h_r*units.kpc, maxiter=20, tol=0.001, normalize=False, ro=8*units.kpc, vo=220*units.km/units.s, new=True, glorder=100)
    update_rot_curve() 
def MN_ed_a_s_func(val):
    global EX_Disk_p, amp4,h_r
    h_r=val*1
    EX_Disk_p = RazorThinExponentialDiskPotential(amp=amp4*(units.Msun/(units.pc**2)), hr=val*units.kpc, maxiter=20, tol=0.001, normalize=False, ro=8*units.kpc, vo=220*units.km/units.s, new=True, glorder=100)
    update_rot_curve()
    

def update_rot_curve():
    ax.cla()
    global MN_b_plot, MN_Bulge_p, MN_Thin_Disk_p,MN_Thick_Disk_p, MN_td_plot,MN_tkd_plot, NFW_p, NFW_plot, EX_d_plot, EX_Disk_p
    composite_pot_array=[]
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
    v_circ_comp = calcRotcurve(composite_pot_array, lista, phi=None)*220
    v_circ_comp_plot, = ax.plot(lista, v_circ_comp, c='k')

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
v_circ_comp_plot, = ax.plot(lista, v_circ_comp, c='k')



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
    
    
check.on_clicked(check_on_clicked)
plt.ylim([0,300])
plt.show()