# -*- coding: utf-8 -*-
# Copyright (C) 2017 - Michael Baudin

import pylab as pl
import numpy as np
#

def alti(t,c):
    g=9.81
    z0=100.
    v0=55.
    m=80.
    tau=m/c
    vinf=-m*g/c
    z=z0+vinf*t+tau*(v0-vinf)*(1-np.exp(-t/tau))
    return np.maximum(z,0.)

def altiZeroDrag(t):
    g=9.81
    z0=100.
    v0=55.
    z=z0+v0*t-g*t**2
    return np.maximum(z,0.)

# 2. Graphique
pl.figure()
t=np.linspace(0,12,100)
z=alti(t,15.)
pl.plot(t,z,"-",label="c=15")
z=alti(t,1.)
pl.plot(t,z,"-",label="c=1")
z=altiZeroDrag(t)
pl.plot(t,z,"-",label="c=0")
if (True):
    pl.xlabel("Time (s)")
    pl.ylabel("Height (m)")
    pl.title("Freefall in a viscous fluid")
else:
    pl.xlabel("Temps (s)")
    pl.ylabel("Altitude (m)")
    pl.title("Chute dans un fluide visqueux")
pl.legend()
pl.show()
