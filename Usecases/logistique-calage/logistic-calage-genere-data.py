# -*- coding: utf-8 -*-
# Copyright (C) 2016 - Michael Baudin

import openturns as ot
from numpy import array, arange


'''
Logistic Growth

Reference
Differential equations, 4th ed., Braun, 1993, 
TAM  Chap.1, "First order differential equations"   
p. 28
'''
# Donnees reelles
ustime=arange(1790,2001,10);
uspop=array([3.9,5.3,7.2,9.6,13.,17.,23.,31.,39.,\
 50.,62.,76.,92.,106.,123.,132.,151.,179.,\
 203.,221.,250.,281.]);

sampleSize = len(ustime)

observedSample = ot.Sample(sampleSize,2)
observedSample.setDescription(["Date (Annees)","Population (Millions)"])
observedSample[:,0] = ustime.reshape((sampleSize,1))
observedSample[:,1] = uspop.reshape((sampleSize,1))

observedSample.exportToCSVFile("calage-logistique-observations.csv")
