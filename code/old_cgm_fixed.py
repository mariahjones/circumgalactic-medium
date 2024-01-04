#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 22 12:39:53 2018

@author: ebuieii
"""

import yt
from yt.units import *
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import h5py
import subprocess
import scipy.integrate as integrate



#Various constants
he_h = 10**10.925/(10**12)
c_h = 10**8.39/(10**12)
n_h = 10**7.86/(10**12)
o_h = 10**8.73/(10**12)
ne_h = 10**8.05/(10**12)
na_h = 10**6.29/(10**12)
mg_h = 10**7.54/(10**12)
si_h = 10**7.52/(10**12)
s_h = 10**7.16/(10**12)
ca_h = 10**6.31/(10**12)
fe_h = 10**7.46/(10**12)
meta = 1.0
ph54 = 4.171475061185105E-005
cj21 = 8.23e-3#2.4e-2#4.23e-3#515.0
mass = 1e12*1.98e33 #in grams
kboltz = 1.38065e-16
mh = 1.67262e-24

#yt fields that are contained in MAIHEM files, can be manipulated to find other quantities (e.g. numbder density, see ndens function)
st = "cell"
def _hspec(field,data):
    return data['h   ']
yt.add_field("hspec",function=_hspec,sampling_type=st)
def _hpspec(field,data):
    return data['hp  ']
yt.add_field("hpspec",function=_hpspec,sampling_type=st)

def _hespec(field,data):
    return data['he  ']
yt.add_field("hespec",function=_hespec,sampling_type=st)
def _hepspec(field,data):
    return data['hep ']
yt.add_field("hepspec",function=_hepspec,sampling_type=st)
def _he2pspec(field,data):
    return data['he2p']
yt.add_field("he2pspec",function=_he2pspec,sampling_type=st)

def _cspec(field,data):
    return data['c   ']
yt.add_field("cspec",function=_cspec,sampling_type=st)
def _cpspec(field,data):
    return data['cp  ']
yt.add_field("cpspec",function=_cpspec,sampling_type=st)
def _c2pspec(field,data):
    return data['c2p ']
yt.add_field("c2pspec",function=_c2pspec,sampling_type=st)
def _c3pspec(field,data):
    return data['c3p ']
yt.add_field("c3pspec",function=_c3pspec,sampling_type=st)
def _c4pspec(field,data):
    return data['c4p ']
yt.add_field("c4pspec",function=_c4pspec,sampling_type=st)
def _c5pspec(field,data):
    return data['c5p ']
yt.add_field("c5pspec",function=_c5pspec,sampling_type=st)

def _nspec(field,data):
    return data['n   ']
yt.add_field("nspec",function=_nspec,sampling_type=st)
def _npspec(field,data):
    return data['np  ']
yt.add_field("npspec",function=_npspec,sampling_type=st)
def _n2pspec(field,data):
    return data['n2p ']
yt.add_field("n2pspec",function=_n2pspec,sampling_type=st)
def _n3pspec(field,data):
    return data['n3p ']
yt.add_field("n3pspec",function=_n3pspec,sampling_type=st)
def _n4pspec(field,data):
    return data['n4p ']
yt.add_field("n4pspec",function=_n4pspec,sampling_type=st)
def _n5pspec(field,data):
    return data['n5p ']
yt.add_field("n5pspec",function=_n5pspec,sampling_type=st)
def _n6pspec(field,data):
    return data['n6p ']
yt.add_field("n6pspec",function=_n6pspec,sampling_type=st)

def _ospec(field,data):
    return data['o   ']
yt.add_field("ospec",function=_ospec,sampling_type=st)
def _opspec(field,data):
    return data['op  ']
yt.add_field("opspec",function=_opspec,sampling_type=st)
def _o2pspec(field,data):
    return data['o2p ']
yt.add_field("o2pspec",function=_o2pspec,sampling_type=st)
def _o3pspec(field,data):
    return data['o3p ']
yt.add_field("o3pspec",function=_o3pspec,sampling_type=st)
def _o4pspec(field,data):
    return data['o4p ']
yt.add_field("o4pspec",function=_o4pspec,sampling_type=st)
def _o5pspec(field,data):
    return data['o5p ']
yt.add_field("o5pspec",function=_o5pspec,sampling_type=st)
def _o6pspec(field,data):
    return data['o6p ']
yt.add_field("o6pspec",function=_o6pspec,sampling_type=st)
def _o7pspec(field,data):
    return data['o7p ']
yt.add_field("o7pspec",function=_o7pspec,sampling_type=st)

def _nespec(field,data):
    return data['ne  ']
yt.add_field("nespec",function=_nespec,sampling_type=st)
def _nepspec(field,data):
    return data['nep ']
yt.add_field("nepspec",function=_nepspec,sampling_type=st)
def _ne2pspec(field,data):
    return data['ne2p']
yt.add_field("ne2pspec",function=_ne2pspec,sampling_type=st)
def _ne3pspec(field,data):
    return data['ne3p']
yt.add_field("ne3pspec",function=_ne3pspec,sampling_type=st)
def _ne4pspec(field,data):
    return data['ne4p']
yt.add_field("ne4pspec",function=_ne4pspec,sampling_type=st)
def _ne5pspec(field,data):
    return data['ne5p']
yt.add_field("ne5pspec",function=_ne5pspec,sampling_type=st)
def _ne6pspec(field,data):
    return data['ne6p']
yt.add_field("ne6pspec",function=_ne6pspec,sampling_type=st)
def _ne7pspec(field,data):
    return data['ne7p']
yt.add_field("ne7pspec",function=_ne7pspec,sampling_type=st)
def _ne8pspec(field,data):
    return data['ne8p']
yt.add_field("ne8pspec",function=_ne8pspec,sampling_type=st)
def _ne9pspec(field,data):
    return data['ne9p']
yt.add_field("ne9pspec",function=_ne9pspec,sampling_type=st)

def _naspec(field,data):
    return data['na  ']
yt.add_field("naspec",function=_naspec,sampling_type=st)
def _napspec(field,data):
    return data['nap ']
yt.add_field("napspec",function=_napspec,sampling_type=st)
def _na2pspec(field,data):
    return data['na2p']
yt.add_field("na2pspec",function=_na2pspec,sampling_type=st)

def _mgspec(field,data):
    return data['mg  ']
yt.add_field("mgspec",function=_mgspec,sampling_type=st)
def _mgpspec(field,data):
    return data['mgp ']
yt.add_field("mgpspec",function=_mgpspec,sampling_type=st)
def _mg2pspec(field,data):
    return data['mg2p']
yt.add_field("mg2pspec",function=_mg2pspec,sampling_type=st)
def _mg3pspec(field,data):
    return data['mg3p']
yt.add_field("mg3pspec",function=_mg3pspec,sampling_type=st)

def _sispec(field,data):
    return data['si  ']
yt.add_field("sispec",function=_sispec,sampling_type=st)
def _sipspec(field,data):
    return data['sip ']
yt.add_field("sipspec",function=_sipspec,sampling_type=st)
def _si2pspec(field,data):
    return data['si2p']
yt.add_field("si2pspec",function=_si2pspec,sampling_type=st)
def _si3pspec(field,data):
    return data['si3p']
yt.add_field("si3pspec",function=_si3pspec,sampling_type=st)
def _si4pspec(field,data):
    return data['si4p']
yt.add_field("si4pspec",function=_si4pspec,sampling_type=st)
def _si5pspec(field,data):
    return data['si5p']
yt.add_field("si5pspec",function=_si5pspec,sampling_type=st)

def _sspec(field,data):
    return data['s   ']
yt.add_field("sspec",function=_sspec,sampling_type=st)
def _spspec(field,data):
    return data['sp  ']
yt.add_field("spspec",function=_spspec,sampling_type=st)
def _s2pspec(field,data):
    return data['s2p ']
yt.add_field("s2pspec",function=_s2pspec,sampling_type=st)
def _s3pspec(field,data):
    return data['s3p ']
yt.add_field("s3pspec",function=_s3pspec,sampling_type=st)
def _s4pspec(field,data):
    return data['s4p ']
yt.add_field("s4pspec",function=_s4pspec,sampling_type=st)

def _caspec(field,data):
    return data['ca  ']
yt.add_field("caspec",function=_caspec,sampling_type=st)
def _capspec(field,data):
    return data['cap ']
yt.add_field("capspec",function=_capspec,sampling_type=st)
def _ca2pspec(field,data):
    return data['ca2p']
yt.add_field("ca2pspec",function=_ca2pspec,sampling_type=st)
def _ca3pspec(field,data):
    return data['ca3p']
yt.add_field("ca3pspec",function=_ca3pspec,sampling_type=st)
def _ca4pspec(field,data):
    return data['ca4p']
yt.add_field("ca4pspec",function=_ca4pspec,sampling_type=st)

def _fespec(field,data):
    return data['fe  ']
yt.add_field("fespec",function=_fespec,sampling_type=st)
def _fepspec(field,data):
    return data['fep ']
yt.add_field("fepspec",function=_fepspec,sampling_type=st)
def _fe2pspec(field,data):
    return data['fe2p']
yt.add_field("fe2pspec",function=_fe2pspec,sampling_type=st)
def _fe3pspec(field,data):
    return data['fe3p']
yt.add_field("fe3pspec",function=_fe3pspec,sampling_type=st)
def _fe4pspec(field,data):
    return data['fe4p']
yt.add_field("fe4pspec",function=_fe4pspec,sampling_type=st)

def _elecspec(field,data):
    return data['elec']
yt.add_field("elecspec",function=_elecspec,sampling_type=st)
def _myabar(field,data):
    abar = (data['hspec'] + data['hpspec'])/1.0 + (data['hespec'] + data['hepspec'] + data['he2pspec'])/4.0 + (data['cspec'] + data['cpspec'] + data['c2pspec'] + data['c3pspec'] + data['c4pspec'] + data['c5pspec'])/12.0 + (data['nspec'] + data['npspec'] + data['n2pspec'] + data['n3pspec'] + data['n4pspec'] + data['n5pspec'] + data['n6pspec'])/14.0  + (data['ospec'] + data['opspec'] + data['o2pspec'] + data['o3pspec'] + data['o4pspec'] + data['o5pspec'] + data['o6pspec'] + data['o7pspec'])/16.0 + (data['nespec'] + data['nepspec'] + data['ne2pspec'] + data['ne3pspec'] + data['ne4pspec'] + data['ne5pspec'] + data['ne6pspec'] + data['ne7pspec'] + data['ne8pspec'] + data['ne9pspec'])/20.0 + (data['naspec'] + data['napspec'] + data['na2pspec'])/22.0 + (data['mgspec'] + data['mgpspec'] + data['mg2pspec'] + data['mg3pspec'])/24.0 + (data['sispec'] + data['sipspec'] + data['si2pspec'] + data['si3pspec'] + data['si4pspec'] + data['si5pspec'])/28.0  + (data['sspec'] + data['spspec'] + data['s2pspec'] + data['s3pspec'] + data['s4pspec'])/32.0  + (data['caspec'] + data['capspec'] + data['ca2pspec'] + data['ca3pspec'] + data['ca4pspec'])/40.0 + (data['fespec'] + data['fepspec'] + data['fe2pspec'] + data['fe3pspec'] + data['fe4pspec'])/56.0   + data['elecspec']/0.000549
    abar = 1.0/abar
    return abar
yt.add_field("myabar",function=_myabar,sampling_type=st)
def _mysoundspeed(field,data):
    cs = 1.66666667*kboltz*data['temp']/(mh*data['myabar'])
    cs = np.sqrt(cs)
    return cs
yt.add_field("mysoundspeed",function=_mysoundspeed,units="cm/s",sampling_type=st)
def _mymachnumber(field,data):
    v2 = data['velx']**2 + data['vely']**2 + data['velz']**2
    v = np.sqrt(v2)
    mn = v/data['mysoundspeed']
    return mn
yt.add_field("mymachnumber",function=_mymachnumber,sampling_type=st)
def _vinfall(field,data):
    v2 = data['velx']*data['x'] + data['vely']*data['y'] + data['velz']*data['z']
    v = v2/np.sqrt(data['x']**2+data['y']**2+data['z']**2)
    mn = v
    return mn
yt.add_field("vinfall",function=_vinfall,units="km/s",sampling_type=st)
def _sigma(field,data):
    v2 = data['velx']**2 + data['vely']**2 + data['velz']**2
    v = np.sqrt(v2-data['vinfall']**2)
    mn = v/np.sqrt(3)
    return mn
yt.add_field("sigma",function=_sigma,units="km/s",sampling_type=st)
def _velocity(field,data):
    v2 = data['velx']**2 + data['vely']**2 + data['velz']**2
    v = np.sqrt(v2)
    mn = v/np.sqrt(3)
    return mn
yt.add_field("velocity",function=_velocity,units="km/s",sampling_type=st)
def _ndens(field,data):
    return data['density'] / (data['myabar'] * mh * g )
yt.add_field("ndens",function=_ndens,units="cm**-3",sampling_type=st)
def _xvar(field,data):
    avgd = np.average(data['dens'])
    return np.log( data['dens'] / avgd )
yt.add_field("xvar",function=_xvar,sampling_type=st)
#def _radius(field,data):
#    rad = data['x']**2 + data['y']**2 + data['z']**2
#    ans = np.sqrt(rad)
#    return ans
#yt.add_field("radius",function=_radius, units="cm")
def _turb_kinetic(field,data):
    v2 = np.sqrt(3)*data['sigma']**2
    v = 0.5*v2
    return v
yt.add_field("turb_kinetic",function=_turb_kinetic,units="cm**2/s**2",sampling_type=st)
def _fall_kinetic(field,data):
    v2 = data['vinfall']**2
    v = 0.5*data['cell_mass']*v2
    return v
yt.add_field("fall_kinetic",function=_fall_kinetic,units="erg",sampling_type=st)
def _grav_accel(field,data):
    ans = G*1e12*1.98e33*gram/data['radius']**2
    return ans
yt.add_field("grav_accel",function=_grav_accel,units="cm/(s**2)",sampling_type=st)
def _turb_vel(field,data):
    vol = data['dx']*data['dy']*data['dz']
    #print len(vol), type(vol), type(np.sum(vol))
    v_vel = np.sum( vol * (data['velx']*data['velx']+data['vely']*data['vely']+data['velz']*data['velz']) ) / np.sum(vol)
    v_vel = np.sqrt(v_vel)
    #d_vel = np.sum( (dd['density']) * (dd['velx']*dd['velx']+dd['vely']*dd['vely']+dd['velz']*dd['velz']) ) / np.sum(dd['density'])
    #d_vel = np.sqrt(d_vel)
    #turb_accel = (data['accx']**2 + data['accy']**2 + data['accz']**2)*centimeter**2/second**4
    #ans = np.sqrt(turb_accel)#*1.59e16*second #this is now velocity and that number is the decay timescale for turbulence
    #ans = ans/np.sqrt(3)
    return v_vel
yt.add_field("turb_vel",function=_turb_vel,units="km/s",sampling_type=st)
def _v_rad(field,data):
    v2 = data['velx']*data['x'] + data['vely']*data['y'] + data['velz']*data['z']
    v = v2/np.sqrt(data['x']**2+data['y']**2+data['z']**2)
    mn = v
    return mn
yt.add_field("v_rad",function=_v_rad,units="km/s",sampling_type=st)
    
#----HERE ARE THE NORMAL MASS FRACTIONS-------!!!!!!!!!!!!!!!

def _hfrac(field,data):
	return data['hspec']/(data['hspec']+data['hpspec'])
yt.add_field('hfrac',function=_hfrac,sampling_type=st)
def _hpfrac(field,data):
	return data['hpspec']/(data['hspec']+data['hpspec'])
yt.add_field('hpfrac',function=_hpfrac,sampling_type=st)

def _hefrac(field,data):
	return he_h*data['hespec']/(data['hespec']+data['hepspec']+data['he2pspec'])
yt.add_field('hefrac',function=_hefrac,sampling_type=st)
def _hepfrac(field,data):
	return he_h*data['hepspec']/(data['hespec']+data['hepspec']+data['he2pspec'])
yt.add_field('hepfrac',function=_hepfrac,sampling_type=st)
def _he2pfrac(field,data):
	return he_h*data['he2pspec']/(data['hespec']+data['hepspec']+data['he2pspec'])
yt.add_field('he2pfrac',function=_he2pfrac,sampling_type=st)

def _cfrac(field,data):
	return meta*c_h*data['cspec']/(data['cspec']+data['cpspec']+data['c2pspec']+data['c3pspec']+data['c4pspec']+data['c5pspec'])
yt.add_field('cfrac',function=_cfrac,sampling_type=st)
def _cpfrac(field,data):
	return meta*c_h**data['cpspec']/(data['cspec']+data['cpspec']+data['c2pspec']+data['c3pspec']+data['c4pspec']+data['c5pspec'])
yt.add_field('cpfrac',function=_cpfrac,sampling_type=st)
def _c2pfrac(field,data):
	return meta*c_h*data['c2pspec']/(data['cspec']+data['cpspec']+data['c2pspec']+data['c3pspec']+data['c4pspec']+data['c5pspec'])
yt.add_field('c2frac',function=_c2pfrac,sampling_type=st)
def _c3pfrac(field,data):
	return meta*c_h*data['c3pspec']/(data['cspec']+data['cpspec']+data['c2pspec']+data['c3pspec']+data['c4pspec']+data['c5pspec'])
yt.add_field('c3pfrac',function=_c3pfrac,sampling_type=st)
def _c4pfrac(field,data):
	return meta*c_h*data['c4pspec']/(data['cspec']+data['cpspec']+data['c2pspec']+data['c3pspec']+data['c4pspec']+data['c5pspec'])
yt.add_field('c4pfrac',function=_c4pfrac,sampling_type=st)
def _c5pfrac(field,data):
	return meta*c_h*data['c5pspec']/(data['cspec']+data['cpspec']+data['c2pspec']+data['c3pspec']+data['c4pspec']+data['c5pspec'])
yt.add_field('c5pfrac',function=_c5pfrac,sampling_type=st)

def _nfrac(field,data):
	return meta*n_h*data['nspec']/(data['nspec']+data['npspec']+data['n2pspec']+data['n3pspec']+data['n4pspec']+data['n5pspec']+data['n6pspec'])
yt.add_field('nfrac',function=_nfrac,sampling_type=st)
def _npfrac(field,data):
	return meta*n_h*data['npspec']/(data['nspec']+data['npspec']+data['n2pspec']+data['n3pspec']+data['n4pspec']+data['n5pspec']+data['n6pspec'])
yt.add_field('npfrac',function=_npfrac,sampling_type=st)
def _n2pfrac(field,data):
	return meta*n_h*data['n2pspec']/(data['nspec']+data['npspec']+data['n2pspec']+data['n3pspec']+data['n4pspec']+data['n5pspec']+data['n6pspec'])
yt.add_field('n2pfrac',function=_n2pfrac,sampling_type=st)
def _n3pfrac(field,data):
	return meta*n_h*data['n3pspec']/(data['nspec']+data['npspec']+data['n2pspec']+data['n3pspec']+data['n4pspec']+data['n5pspec']+data['n6pspec'])
yt.add_field('n3pfrac',function=_n3pfrac,sampling_type=st)
def _n4pfrac(field,data):
	return meta*n_h*data['n4pspec']/(data['nspec']+data['npspec']+data['n2pspec']+data['n3pspec']+data['n4pspec']+data['n5pspec']+data['n6pspec'])
yt.add_field('n4pfrac',function=_n4pfrac,sampling_type=st)
def _n5pfrac(field,data):
	return meta*n_h*data['n5pspec']/(data['nspec']+data['npspec']+data['n2pspec']+data['n3pspec']+data['n4pspec']+data['n5pspec']+data['n6pspec'])
yt.add_field('n5pfrac',function=_n5pfrac,sampling_type=st)
def _n6pfrac(field,data):
	return meta*n_h*data['n6pspec']/(data['nspec']+data['npspec']+data['n2pspec']+data['n3pspec']+data['n4pspec']+data['n5pspec']+data['n6pspec'])
yt.add_field('n6pfrac',function=_n6pfrac,sampling_type=st)

def _ofrac(field,data):
	return meta*o_h*data['ospec']/(data['ospec']+data['opspec']+data['o2pspec']+data['o3pspec']+data['o4pspec']+data['o5pspec']+data['o6pspec']+data['o7pspec'])
yt.add_field('ofrac',function=_ofrac,sampling_type=st)
def _opfrac(field,data):
	return meta*o_h*data['opspec']/(data['ospec']+data['opspec']+data['o2pspec']+data['o3pspec']+data['o4pspec']+data['o5pspec']+data['o6pspec']+data['o7pspec'])
yt.add_field('opfrac',function=_opfrac,sampling_type=st)
def _o2pfrac(field,data):
	return meta*o_h*data['o2pspec']/(data['ospec']+data['opspec']+data['o2pspec']+data['o3pspec']+data['o4pspec']+data['o5pspec']+data['o6pspec']+data['o7pspec'])
yt.add_field('o2pfrac',function=_o2pfrac,sampling_type=st)
def _o3pfrac(field,data):
	return meta*o_h*data['o3pspec']/(data['ospec']+data['opspec']+data['o2pspec']+data['o3pspec']+data['o4pspec']+data['o5pspec']+data['o6pspec']+data['o7pspec'])
yt.add_field('o3pfrac',function=_o3pfrac,sampling_type=st)
def _o4pfrac(field,data):
	return meta*o_h*data['o4pspec']/(data['ospec']+data['opspec']+data['o2pspec']+data['o3pspec']+data['o4pspec']+data['o5pspec']+data['o6pspec']+data['o7pspec'])
yt.add_field('o4pfrac',function=_o4pfrac,sampling_type=st)
def _o5pfrac(field,data):
	return meta*o_h*data['o5pspec']/(data['ospec']+data['opspec']+data['o2pspec']+data['o3pspec']+data['o4pspec']+data['o5pspec']+data['o6pspec']+data['o7pspec'])
yt.add_field('o5pfrac',function=_o5pfrac,sampling_type=st)
def _o6pfrac(field,data):
	return meta*o_h*data['o6pspec']/(data['ospec']+data['opspec']+data['o2pspec']+data['o3pspec']+data['o4pspec']+data['o5pspec']+data['o6pspec']+data['o7pspec'])
yt.add_field('o6pfrac',function=_o6pfrac,sampling_type=st)
def _o7pfrac(field,data):
	return meta*o_h*data['o7pspec']/(data['ospec']+data['opspec']+data['o2pspec']+data['o3pspec']+data['o4pspec']+data['o5pspec']+data['o6pspec']+data['o7pspec'])
yt.add_field('o7pfrac',function=_o7pfrac,sampling_type=st)

def _nefrac(field,data):
	return meta*ne_h*data['nespec']/(data['nespec']+data['nepspec']+data['ne2pspec']+data['ne3pspec']+data['ne4pspec']+data['ne5pspec']+data['ne6pspec']+data['ne7pspec']+data['ne8pspec']+data['ne9pspec'])
yt.add_field('nefrac',function=_nefrac,sampling_type=st)
def _nepfrac(field,data):
	return meta*ne_h*data['nepspec']/(data['nespec']+data['nepspec']+data['ne2pspec']+data['ne3pspec']+data['ne4pspec']+data['ne5pspec']+data['ne6pspec']+data['ne7pspec']+data['ne8pspec']+data['ne9pspec'])
yt.add_field('nepfrac',function=_nepfrac,sampling_type=st)
def _ne2pfrac(field,data):
	return meta*ne_h*data['ne2pspec']/(data['nespec']+data['nepspec']+data['ne2pspec']+data['ne3pspec']+data['ne4pspec']+data['ne5pspec']+data['ne6pspec']+data['ne7pspec']+data['ne8pspec']+data['ne9pspec'])
yt.add_field('ne2pfrac',function=_ne2pfrac,sampling_type=st)
def _ne3pfrac(field,data):
	return meta*ne_h*data['ne3pspec']/(data['nespec']+data['nepspec']+data['ne2pspec']+data['ne3pspec']+data['ne4pspec']+data['ne5pspec']+data['ne6pspec']+data['ne7pspec']+data['ne8pspec']+data['ne9pspec'])
yt.add_field('ne3pfrac',function=_ne3pfrac,sampling_type=st)
def _ne4pfrac(field,data):
	return meta*ne_h*data['ne4pspec']/(data['nespec']+data['nepspec']+data['ne2pspec']+data['ne3pspec']+data['ne4pspec']+data['ne5pspec']+data['ne6pspec']+data['ne7pspec']+data['ne8pspec']+data['ne9pspec'])
yt.add_field('ne4pfrac',function=_ne4pfrac,sampling_type=st)
def _ne5pfrac(field,data):
	return meta*ne_h*data['ne5pspec']/(data['nespec']+data['nepspec']+data['ne2pspec']+data['ne3pspec']+data['ne4pspec']+data['ne5pspec']+data['ne6pspec']+data['ne7pspec']+data['ne8pspec']+data['ne9pspec'])
yt.add_field('ne5pfrac',function=_ne5pfrac,sampling_type=st)
def _ne6pfrac(field,data):
	return meta*ne_h*data['ne6pspec']/(data['nespec']+data['nepspec']+data['ne2pspec']+data['ne3pspec']+data['ne4pspec']+data['ne5pspec']+data['ne6pspec']+data['ne7pspec']+data['ne8pspec']+data['ne9pspec'])
yt.add_field('ne6pfrac',function=_ne6pfrac,sampling_type=st)
def _ne7pfrac(field,data):
	return meta*ne_h*data['ne7pspec']/(data['nespec']+data['nepspec']+data['ne2pspec']+data['ne3pspec']+data['ne4pspec']+data['ne5pspec']+data['ne6pspec']+data['ne7pspec']+data['ne8pspec']+data['ne9pspec'])
yt.add_field('ne7pfrac',function=_ne7pfrac,sampling_type=st)
def _ne8pfrac(field,data):
	return meta*ne_h*data['ne8pspec']/(data['nespec']+data['nepspec']+data['ne2pspec']+data['ne3pspec']+data['ne4pspec']+data['ne5pspec']+data['ne6pspec']+data['ne7pspec']+data['ne8pspec']+data['ne9pspec'])
yt.add_field('ne8pfrac',function=_ne8pfrac,sampling_type=st)
def _ne9pfrac(field,data):
	return meta*ne_h*data['ne9pspec']/(data['nespec']+data['nepspec']+data['ne2pspec']+data['ne3pspec']+data['ne4pspec']+data['ne5pspec']+data['ne6pspec']+data['ne7pspec']+data['ne8pspec']+data['ne9pspec'])
yt.add_field('ne9pfrac',function=_nepfrac,sampling_type=st)

def _nafrac(field,data):
	return meta*na_h*data['naspec']/(data['naspec']+data['napspec']+data['na2pspec'])
yt.add_field('nafrac',function=_nafrac,sampling_type=st)
def _napfrac(field,data):
	return meta*na_h*data['napspec']/(data['naspec']+data['napspec']+data['na2pspec'])
yt.add_field('napspec',function=_napfrac,sampling_type=st)
def _na2pfrac(field,data):
	return meta*na_h*data['na2pspec']/(data['naspec']+data['napspec']+data['na2pspec'])
yt.add_field('na2pspec',function=_na2pfrac,sampling_type=st)

def _mgfrac(field,data):
	return meta*mg_h*data['mgspec']/(data['mgspec']+data['mgpspec']+data['mg2pspec']+data['mg3pspec'])
yt.add_field('mgfrac',function=_mgfrac,sampling_type=st)
def _mgpfrac(field,data):
	return meta*mg_h*data['mgpspec']/(data['mgspec']+data['mgpspec']+data['mg2pspec']+data['mg3pspec'])
yt.add_field('mgpfrac',function=_mgpfrac,sampling_type=st)
def _mg2pfrac(field,data):
	return meta*mg_h*data['mg2pspec']/(data['mgspec']+data['mgpspec']+data['mg2pspec']+data['mg3pspec'])
yt.add_field('mg2pfrac',function=_mg2pfrac,sampling_type=st)
def _mg3pfrac(field,data):
	return meta*mg_h*data['mg3pspec']/(data['mgspec']+data['mgpspec']+data['mg2pspec']+data['mg3pspec'])
yt.add_field('mg3pfrac',function=_mg3pfrac,sampling_type=st)

def _sifrac(field,data):
	return meta*si_h*data['sispec']/(data['sispec']+data['sipspec']+data['si2pspec']+data['si3pspec']+data['si4pspec']+data['si5pspec'])
yt.add_field('sifrac',function=_sifrac,sampling_type=st)
def _sipfrac(field,data):
	return meta*si_h*data['sipspec']/(data['sispec']+data['sipspec']+data['si2pspec']+data['si3pspec']+data['si4pspec']+data['si5pspec'])
yt.add_field('sipfrac',function=_sipfrac,sampling_type=st)
def _si2pfrac(field,data):
	return meta*si_h*data['si2pspec']/(data['sispec']+data['sipspec']+data['si2pspec']+data['si3pspec']+data['si4pspec']+data['si5pspec'])
yt.add_field('si2pfrac',function=_si2pfrac,sampling_type=st)
def _si3pfrac(field,data):
	return meta*si_h*data['si3pspec']/(data['sispec']+data['sipspec']+data['si2pspec']+data['si3pspec']+data['si4pspec']+data['si5pspec'])
yt.add_field('si3pfrac',function=_si3pfrac,sampling_type=st)
def _si4pfrac(field,data):
	return meta*si_h*data['si4pspec']/(data['sispec']+data['sipspec']+data['si2pspec']+data['si3pspec']+data['si4pspec']+data['si5pspec'])
yt.add_field('si4pfrac',function=_si4pfrac,sampling_type=st)
def _si5pfrac(field,data):
	return meta*si_h*data['si5pspec']/(data['sispec']+data['sipspec']+data['si2pspec']+data['si3pspec']+data['si4pspec']+data['si5pspec'])
yt.add_field('si5pfrac',function=_si5pfrac,sampling_type=st)

def _sfrac(field,data):
	return meta*s_h*data['sspec']/(data['sspec']+data['spspec']+data['s2pspec']+data['s3pspec']+data['s4pspec'])
yt.add_field('sfrac',function=_sfrac,sampling_type=st)
def _spfrac(field,data):
	return meta*s_h*data['spspec']/(data['sspec']+data['spspec']+data['s2pspec']+data['s3pspec']+data['s4pspec'])
yt.add_field('spfrac',function=_spfrac,sampling_type=st)
def _s2pfrac(field,data):
	return meta*s_h*data['s2pspec']/(data['sspec']+data['spspec']+data['s2pspec']+data['s3pspec']+data['s4pspec'])
yt.add_field('s2pfrac',function=_s2pfrac,sampling_type=st)
def _s3pfrac(field,data):
	return meta*s_h*data['s3pspec']/(data['sspec']+data['spspec']+data['s2pspec']+data['s3pspec']+data['s4pspec'])
yt.add_field('s3pfrac',function=_s3pfrac,sampling_type=st)
def _s4pfrac(field,data):
	return meta*s_h*data['s4pspec']/(data['sspec']+data['spspec']+data['s2pspec']+data['s3pspec']+data['s4pspec'])
yt.add_field('s4pfrac',function=_s4pfrac,sampling_type=st)

def _cafrac(field,data):
	return meta*ca_h*data['caspec']/(data['caspec']+data['capspec']+data['ca2pspec']+data['ca3pspec']+data['ca4pspec'])
yt.add_field('cafrac',function=_cafrac,sampling_type=st)
def _capfrac(field,data):
	return meta*ca_h*data['capspec']/(data['caspec']+data['capspec']+data['ca2pspec']+data['ca3pspec']+data['ca4pspec'])
yt.add_field('capfrac',function=_capfrac,sampling_type=st)
def _ca2pfrac(field,data):
	return meta*ca_h*data['ca2pspec']/(data['caspec']+data['capspec']+data['ca2pspec']+data['ca3pspec']+data['ca4pspec'])
yt.add_field('ca2pfrac',function=_ca2pfrac,sampling_type=st)
def _ca3pfrac(field,data):
	return meta*ca_h*data['ca3pspec']/(data['caspec']+data['capspec']+data['ca2pspec']+data['ca3pspec']+data['ca4pspec'])
yt.add_field('ca3pfrac',function=_ca3pfrac,sampling_type=st)
def _ca4pfrac(field,data):
	return meta*ca_h*data['ca4pspec']/(data['caspec']+data['capspec']+data['ca2pspec']+data['ca3pspec']+data['ca4pspec'])
yt.add_field('ca4pfrac',function=_ca4pfrac,sampling_type=st)

def _fefrac(field,data):
	return meta*fe_h*data['fespec']/(data['fespec']+data['fepspec']+data['fe2pspec']+data['fe3pspec']+data['fe4pspec'])
yt.add_field('fefrac',function=_fefrac,sampling_type=st)
def _fepfrac(field,data):
	return meta*fe_h*data['fepspec']/(data['fespec']+data['fepspec']+data['fe2pspec']+data['fe3pspec']+data['fe4pspec'])
yt.add_field('fepfrac',function=_fepfrac,sampling_type=st)
def _fe2pfrac(field,data):
	return meta*fe_h*data['fe2pspec']/(data['fespec']+data['fepspec']+data['fe2pspec']+data['fe3pspec']+data['fe4pspec'])
yt.add_field('fe2pfrac',function=_fe2pfrac,sampling_type=st)
def _fe3pfrac(field,data):
	return meta*fe_h*data['fe3pspec']/(data['fespec']+data['fepspec']+data['fe2pspec']+data['fe3pspec']+data['fe4pspec'])
yt.add_field('fe3pfrac',function=_fe3pfrac,sampling_type=st)
def _fe4pfrac(field,data):
	return meta*fe_h*data['fe4pspec']/(data['fespec']+data['fepspec']+data['fe2pspec']+data['fe3pspec']+data['fe4pspec'])
yt.add_field('fe4pfrac',function=_fe4pfrac,sampling_type=st)
#--------HERE ARE THE NUMBER DENSITIES FOR MASS FRACTIONS, FOR COLUMN DENSITIES, MAKE PROJECTION PLOTS, 10^12 - 10^15 FOR THE IONS, UP TO 10^20 FOR H ----------------!!!!!!!!!!!!!!!!!!!!	

def _hdens(field,data):
	return data['ndens']*data['hspec']/(data['hspec']+data['hpspec'])
yt.add_field('hdens',function=_hdens,units="cm**-3",sampling_type=st)
def _hpdens(field,data):
	return data['ndens']*data['hpspec']/(data['hspec']+data['hpspec'])
yt.add_field('hpdens',function=_hpdens,units="cm**-3",sampling_type=st)

def _hedens(field,data):
	return data['ndens']*he_h*data['hespec']/(data['hespec']+data['hepspec']+data['he2pspec'])
yt.add_field('hedens',function=_hedens,units="cm**-3",sampling_type=st)
def _hepdens(field,data):
	return data['ndens']*he_h*data['hepspec']/(data['hespec']+data['hepspec']+data['he2pspec'])
yt.add_field('hepdens',function=_hepdens,units="cm**-3",sampling_type=st)
def _he2pdens(field,data):
	return data['ndens']*he_h*data['he2pspec']/(data['hespec']+data['hepspec']+data['he2pspec'])
yt.add_field('he2pdens',function=_he2pdens,units="cm**-3",sampling_type=st)

def _cdens(field,data):
	return data['ndens']*meta*c_h*data['cspec']/(data['cspec']+data['cpspec']+data['c2pspec']+data['c3pspec']+data['c4pspec']+data['c5pspec'])
yt.add_field('cdens',function=_cdens,units="cm**-3",sampling_type=st)
def _cpdens(field,data):
	return data['ndens']*meta*c_h*data['cpspec']/(data['cspec']+data['cpspec']+data['c2pspec']+data['c3pspec']+data['c4pspec']+data['c5pspec'])
yt.add_field('cpdens',function=_cpdens,units="cm**-3",sampling_type=st)
def _c2pdens(field,data):
	return data['ndens']*meta*c_h*data['c2pspec']/(data['cspec']+data['cpspec']+data['c2pspec']+data['c3pspec']+data['c4pspec']+data['c5pspec'])
yt.add_field('c2pdens',function=_c2pdens,units="cm**-3",sampling_type=st)
def _c3pdens(field,data):
	return data['ndens']*meta*c_h*data['c3pspec']/(data['cspec']+data['cpspec']+data['c2pspec']+data['c3pspec']+data['c4pspec']+data['c5pspec'])
yt.add_field('c3pdens',function=_c3pdens,units="cm**-3",sampling_type=st)
def _c4pdens(field,data):
	return data['ndens']*meta*c_h*data['c4pspec']/(data['cspec']+data['cpspec']+data['c2pspec']+data['c3pspec']+data['c4pspec']+data['c5pspec'])
yt.add_field('c4pdens',function=_c4pdens,units="cm**-3",sampling_type=st)
def _c5pdens(field,data):
	return data['ndens']*meta*c_h*data['c5pspec']/(data['cspec']+data['cpspec']+data['c2pspec']+data['c3pspec']+data['c4pspec']+data['c5pspec'])
yt.add_field('c5pdens',function=_c5pdens,units="cm**-3",sampling_type=st)

def _nidens(field,data):
	return data['ndens']*meta*n_h*data['nspec']/(data['nspec']+data['npspec']+data['n2pspec']+data['n3pspec']+data['n4pspec']+data['n5pspec']+data['n6pspec'])
yt.add_field('nidens',function=_nidens,units="cm**-3",sampling_type=st)
def _npdens(field,data):
	return data['ndens']*meta*n_h*data['npspec']/(data['nspec']+data['npspec']+data['n2pspec']+data['n3pspec']+data['n4pspec']+data['n5pspec']+data['n6pspec'])
yt.add_field('npdens',function=_npdens,units="cm**-3",sampling_type=st)
def _n2pdens(field,data):
	return data['ndens']*meta*n_h*data['n2pspec']/(data['nspec']+data['npspec']+data['n2pspec']+data['n3pspec']+data['n4pspec']+data['n5pspec']+data['n6pspec'])
yt.add_field('n2pdens',function=_n2pdens,units="cm**-3",sampling_type=st)
def _n3pdens(field,data):
	return data['ndens']*meta*n_h*data['n3pspec']/(data['nspec']+data['npspec']+data['n2pspec']+data['n3pspec']+data['n4pspec']+data['n5pspec']+data['n6pspec'])
yt.add_field('n3pdens',function=_n3pdens,units="cm**-3",sampling_type=st)
def _n4pdens(field,data):
	return data['ndens']*meta*n_h*data['n4pspec']/(data['nspec']+data['npspec']+data['n2pspec']+data['n3pspec']+data['n4pspec']+data['n5pspec']+data['n6pspec'])
yt.add_field('n4pdens',function=_n4pdens,units="cm**-3",sampling_type=st)
def _n5pdens(field,data):
	return data['ndens']*meta*n_h*data['n5pspec']/(data['nspec']+data['npspec']+data['n2pspec']+data['n3pspec']+data['n4pspec']+data['n5pspec']+data['n6pspec'])
yt.add_field('n5pfrac',function=_n5pdens,units="cm**-3",sampling_type=st)
def _n6pdens(field,data):
	return data['ndens']*meta*n_h*data['n6pspec']/(data['nspec']+data['npspec']+data['n2pspec']+data['n3pspec']+data['n4pspec']+data['n5pspec']+data['n6pspec'])
yt.add_field('n6pdens',function=_n6pdens,units="cm**-3",sampling_type=st)

def _odens(field,data):
	return data['ndens']*meta*o_h*data['ospec']/(data['ospec']+data['opspec']+data['o2pspec']+data['o3pspec']+data['o4pspec']+data['o5pspec']+data['o6pspec']+data['o7pspec'])
yt.add_field('odens',function=_odens,units="cm**-3",sampling_type=st)
def _opdens(field,data):
	return data['ndens']*meta*o_h*data['opspec']/(data['ospec']+data['opspec']+data['o2pspec']+data['o3pspec']+data['o4pspec']+data['o5pspec']+data['o6pspec']+data['o7pspec'])
yt.add_field('opdens',function=_opdens,units="cm**-3",sampling_type=st)
def _o2pdens(field,data):
	return data['ndens']*meta*o_h*data['o2pspec']/(data['ospec']+data['opspec']+data['o2pspec']+data['o3pspec']+data['o4pspec']+data['o5pspec']+data['o6pspec']+data['o7pspec'])
yt.add_field('o2pdens',function=_o2pdens,units="cm**-3",sampling_type=st)
def _o3pdens(field,data):
	return data['ndens']*meta*o_h*data['o3pspec']/(data['ospec']+data['opspec']+data['o2pspec']+data['o3pspec']+data['o4pspec']+data['o5pspec']+data['o6pspec']+data['o7pspec'])
yt.add_field('o3pdens',function=_o3pdens,units="cm**-3",sampling_type=st)
def _o4pdens(field,data):
	return data['ndens']*meta*o_h*data['o4pspec']/(data['ospec']+data['opspec']+data['o2pspec']+data['o3pspec']+data['o4pspec']+data['o5pspec']+data['o6pspec']+data['o7pspec'])
yt.add_field('o4pdens',function=_o4pdens,units="cm**-3",sampling_type=st)
def _o5pdens(field,data):
	return data['ndens']*meta*o_h*data['o5pspec']/(data['ospec']+data['opspec']+data['o2pspec']+data['o3pspec']+data['o4pspec']+data['o5pspec']+data['o6pspec']+data['o7pspec'])
yt.add_field('o5pdens',function=_o5pdens,units="cm**-3",sampling_type=st)
def _o6pdens(field,data):
	return data['ndens']*meta*o_h*data['o6pspec']/(data['ospec']+data['opspec']+data['o2pspec']+data['o3pspec']+data['o4pspec']+data['o5pspec']+data['o6pspec']+data['o7pspec'])
yt.add_field('o6pdens',function=_o6pdens,units="cm**-3",sampling_type=st)
def _o7pdens(field,data):
	return data['ndens']*meta*o_h*data['o7pspec']/(data['ospec']+data['opspec']+data['o2pspec']+data['o3pspec']+data['o4pspec']+data['o5pspec']+data['o6pspec']+data['o7pspec'])
yt.add_field('o7pdens',function=_o7pdens,units="cm**-3",sampling_type=st)

def _nedens(field,data):
	return data['ndens']*meta*ne_h*data['nespec']/(data['nespec']+data['nepspec']+data['ne2pspec']+data['ne3pspec']+data['ne4pspec']+data['ne5pspec']+data['ne6pspec']+data['ne7pspec']+data['ne8pspec']+data['ne9pspec'])
yt.add_field('nedens',function=_nedens,units="cm**-3",sampling_type=st)
def _nepdens(field,data):
	return data['ndens']*meta*ne_h*data['nepspec']/(data['nespec']+data['nepspec']+data['ne2pspec']+data['ne3pspec']+data['ne4pspec']+data['ne5pspec']+data['ne6pspec']+data['ne7pspec']+data['ne8pspec']+data['ne9pspec'])
yt.add_field('nepdens',function=_nepdens,units="cm**-3",sampling_type=st)
def _ne2pdens(field,data):
	return data['ndens']*meta*ne_h*data['ne2pspec']/(data['nespec']+data['nepspec']+data['ne2pspec']+data['ne3pspec']+data['ne4pspec']+data['ne5pspec']+data['ne6pspec']+data['ne7pspec']+data['ne8pspec']+data['ne9pspec'])
yt.add_field('ne2pdens',function=_ne2pdens,units="cm**-3",sampling_type=st)
def _ne3pdens(field,data):
	return data['ndens']*meta*ne_h*data['ne3pspec']/(data['nespec']+data['nepspec']+data['ne2pspec']+data['ne3pspec']+data['ne4pspec']+data['ne5pspec']+data['ne6pspec']+data['ne7pspec']+data['ne8pspec']+data['ne9pspec'])
yt.add_field('ne3pdens',function=_ne3pfrac,units="cm**-3",sampling_type=st)
def _ne4pdens(field,data):
	return data['ndens']*meta*ne_h*data['ne4pspec']/(data['nespec']+data['nepspec']+data['ne2pspec']+data['ne3pspec']+data['ne4pspec']+data['ne5pspec']+data['ne6pspec']+data['ne7pspec']+data['ne8pspec']+data['ne9pspec'])
yt.add_field('ne4pdens',function=_ne4pdens,units="cm**-3",sampling_type=st)
def _ne5pdens(field,data):
	return data['ndens']*meta*ne_h*data['ne5pspec']/(data['nespec']+data['nepspec']+data['ne2pspec']+data['ne3pspec']+data['ne4pspec']+data['ne5pspec']+data['ne6pspec']+data['ne7pspec']+data['ne8pspec']+data['ne9pspec'])
yt.add_field('ne5pdens',function=_ne5pdens,units="cm**-3",sampling_type=st)
def _ne6pdens(field,data):
	return data['ndens']*meta*ne_h*data['ne6pspec']/(data['nespec']+data['nepspec']+data['ne2pspec']+data['ne3pspec']+data['ne4pspec']+data['ne5pspec']+data['ne6pspec']+data['ne7pspec']+data['ne8pspec']+data['ne9pspec'])
yt.add_field('ne6pdens',function=_ne6pdens,units="cm**-3",sampling_type=st)
def _ne7pdens(field,data):
	return data['ndens']*meta*ne_h*data['ne7pspec']/(data['nespec']+data['nepspec']+data['ne2pspec']+data['ne3pspec']+data['ne4pspec']+data['ne5pspec']+data['ne6pspec']+data['ne7pspec']+data['ne8pspec']+data['ne9pspec'])
yt.add_field('ne7pdens',function=_ne7pdens,units="cm**-3",sampling_type=st)
def _ne8pdens(field,data):
	return data['ndens']*meta*ne_h*data['ne8pspec']/(data['nespec']+data['nepspec']+data['ne2pspec']+data['ne3pspec']+data['ne4pspec']+data['ne5pspec']+data['ne6pspec']+data['ne7pspec']+data['ne8pspec']+data['ne9pspec'])
yt.add_field('ne8pdens',function=_ne8pdens,units="cm**-3",sampling_type=st)
def _ne9pdens(field,data):
	return data['ndens']*meta*ne_h*data['ne9pspec']/(data['nespec']+data['nepspec']+data['ne2pspec']+data['ne3pspec']+data['ne4pspec']+data['ne5pspec']+data['ne6pspec']+data['ne7pspec']+data['ne8pspec']+data['ne9pspec'])
yt.add_field('ne9pdens',function=_ne9pdens,units="cm**-3",sampling_type=st)

def _nadens(field,data):
	return data['ndens']*meta*na_h*data['naspec']/(data['naspec']+data['napspec']+data['na2pspec'])
yt.add_field('nadens',function=_nadens,units="cm**-3",sampling_type=st)
def _napdens(field,data):
	return data['ndens']*meta*na_h*data['napspec']/(data['naspec']+data['napspec']+data['na2pspec'])
yt.add_field('napdens',function=_napdens,units="cm**-3",sampling_type=st)
def _na2pdens(field,data):
	return data['ndens']*meta*na_h*data['na2pspec']/(data['naspec']+data['napspec']+data['na2pspec'])
yt.add_field('na2pdens',function=_na2pdens,units="cm**-3",sampling_type=st)

def _mgdens(field,data):
	return data['ndens']*meta*mg_h*data['mgspec']/(data['mgspec']+data['mgpspec']+data['mg2pspec']+data['mg3pspec'])
yt.add_field('mgdens',function=_mgdens,units="cm**-3",sampling_type=st)
def _mgpdens(field,data):
	return data['ndens']*meta*mg_h*data['mgpspec']/(data['mgspec']+data['mgpspec']+data['mg2pspec']+data['mg3pspec'])
yt.add_field('mgpdens',function=_mgpdens,units="cm**-3",sampling_type=st)
def _mg2pdens(field,data):
	return data['ndens']*meta*mg_h*data['mg2pspec']/(data['mgspec']+data['mgpspec']+data['mg2pspec']+data['mg3pspec'])
yt.add_field('mg2pdens',function=_mg2pdens,units="cm**-3",sampling_type=st)
def _mg3pdens(field,data):
	return data['ndens']*meta*mg_h*data['mg3pspec']/(data['mgspec']+data['mgpspec']+data['mg2pspec']+data['mg3pspec'])
yt.add_field('mg3pdens',function=_mg3pdens,units="cm**-3",sampling_type=st)

def _sidens(field,data):
	return data['ndens']*meta*si_h*data['sispec']/(data['sispec']+data['sipspec']+data['si2pspec']+data['si3pspec']+data['si4pspec']+data['si5pspec'])
yt.add_field('sidens',function=_sidens,units="cm**-3",sampling_type=st)
def _sipdens(field,data):
	return data['ndens']*meta*si_h*data['sipspec']/(data['sispec']+data['sipspec']+data['si2pspec']+data['si3pspec']+data['si4pspec']+data['si5pspec'])
yt.add_field('sipdens',function=_sipdens,units="cm**-3",sampling_type=st)
def _si2pdens(field,data):
	return data['ndens']*meta*si_h*data['si2pspec']/(data['sispec']+data['sipspec']+data['si2pspec']+data['si3pspec']+data['si4pspec']+data['si5pspec'])
yt.add_field('si2pdens',function=_si2pdens,units="cm**-3",sampling_type=st)
def _si3pdens(field,data):
	return data['ndens']*meta*si_h*data['si3pspec']/(data['sispec']+data['sipspec']+data['si2pspec']+data['si3pspec']+data['si4pspec']+data['si5pspec'])
yt.add_field('si3pdens',function=_si3pdens,units="cm**-3",sampling_type=st)
def _si4pdens(field,data):
	return data['ndens']*meta*si_h*data['si4pspec']/(data['sispec']+data['sipspec']+data['si2pspec']+data['si3pspec']+data['si4pspec']+data['si5pspec'])
yt.add_field('si4pdens',function=_si4pdens,units="cm**-3",sampling_type=st)
def _si5pdens(field,data):
	return data['ndens']*meta*si_h*data['si5pspec']/(data['sispec']+data['sipspec']+data['si2pspec']+data['si3pspec']+data['si4pspec']+data['si5pspec'])
yt.add_field('si5pdens',function=_si5pdens,units="cm**-3",sampling_type=st)

def _sdens(field,data):
	return data['ndens']*meta*s_h*data['sspec']/(data['sspec']+data['spspec']+data['s2pspec']+data['s3pspec']+data['s4pspec'])
yt.add_field('sdens',function=_sdens,units="cm**-3",sampling_type=st)
def _spdens(field,data):
	return data['ndens']*meta*s_h*data['spspec']/(data['sspec']+data['spspec']+data['s2pspec']+data['s3pspec']+data['s4pspec'])
yt.add_field('spdens',function=_spdens,units="cm**-3",sampling_type=st)
def _s2pdens(field,data):
	return data['ndens']*meta*s_h*data['s2pspec']/(data['sspec']+data['spspec']+data['s2pspec']+data['s3pspec']+data['s4pspec'])
yt.add_field('s2pdens',function=_s2pdens,units="cm**-3",sampling_type=st)
def _s3pdens(field,data):
	return data['ndens']*meta*s_h*data['s3pspec']/(data['sspec']+data['spspec']+data['s2pspec']+data['s3pspec']+data['s4pspec'])
yt.add_field('s3pdens',function=_s3pdens,units="cm**-3",sampling_type=st)
def _s4pdens(field,data):
	return data['ndens']*meta*s_h*data['s4pspec']/(data['sspec']+data['spspec']+data['s2pspec']+data['s3pspec']+data['s4pspec'])
yt.add_field('s4pdens',function=_s4pdens,units="cm**-3",sampling_type=st)

def _cadens(field,data):
	return data['ndens']*meta*ca_h*data['caspec']/(data['caspec']+data['capspec']+data['ca2pspec']+data['ca3pspec']+data['ca4pspec'])
yt.add_field('cadens',function=_cadens,units="cm**-3",sampling_type=st)
def _capdens(field,data):
	return data['ndens']*meta*ca_h*data['capspec']/(data['caspec']+data['capspec']+data['ca2pspec']+data['ca3pspec']+data['ca4pspec'])
yt.add_field('capdens',function=_capdens,units="cm**-3",sampling_type=st)
def _ca2pdens(field,data):
	return data['ndens']*meta*ca_h*data['ca2pspec']/(data['caspec']+data['capspec']+data['ca2pspec']+data['ca3pspec']+data['ca4pspec'])
yt.add_field('ca2pdens',function=_ca2pdens,units="cm**-3",sampling_type=st)
def _ca3pdens(field,data):
	return data['ndens']*meta*ca_h*data['ca3pspec']/(data['caspec']+data['capspec']+data['ca2pspec']+data['ca3pspec']+data['ca4pspec'])
yt.add_field('ca3pdens',function=_ca3pdens,units="cm**-3",sampling_type=st)
def _ca4pdens(field,data):
	return data['ndens']*meta*ca_h*data['ca4pspec']/(data['caspec']+data['capspec']+data['ca2pspec']+data['ca3pspec']+data['ca4pspec'])
yt.add_field('ca4pdens',function=_ca4pdens,units="cm**-3",sampling_type=st)

def _fedens(field,data):
	return data['ndens']*meta*fe_h*data['fespec']/(data['fespec']+data['fepspec']+data['fe2pspec']+data['fe3pspec']+data['fe4pspec'])
yt.add_field('fedens',function=_fedens,units="cm**-3",sampling_type=st)
def _fepdens(field,data):
	return data['ndens']*meta*fe_h*data['fepspec']/(data['fespec']+data['fepspec']+data['fe2pspec']+data['fe3pspec']+data['fe4pspec'])
yt.add_field('fepdens',function=_fepdens,units="cm**-3",sampling_type=st)
def _fe2pdens(field,data):
	return data['ndens']*meta*fe_h*data['fe2pspec']/(data['fespec']+data['fepspec']+data['fe2pspec']+data['fe3pspec']+data['fe4pspec'])
yt.add_field('fe2pdens',function=_fe2pdens,units="cm**-3",sampling_type=st)
def _fe3pdens(field,data):
	return data['ndens']*meta*fe_h*data['fe3pspec']/(data['fespec']+data['fepspec']+data['fe2pspec']+data['fe3pspec']+data['fe4pspec'])
yt.add_field('fe3pdens',function=_fe3pdens,units="cm**-3",sampling_type=st)
def _fe4pdens(field,data):
	return data['ndens']*meta*fe_h*data['fe4pspec']/(data['fespec']+data['fepspec']+data['fe2pspec']+data['fe3pspec']+data['fe4pspec'])
yt.add_field('fe4pdens',function=_fe4pdens,units="cm**-3",sampling_type=st)

def _Cool(field,data):
    return 9.4e16*cm**2*gram/second**2*(data['lh  ']+data['lhe ']+data['lc  ']+data['ln  ']+data['lo  ']+data['lne ']+data['lna ']+data['lmg ']+data['lsi ']+data['ls  ']+data['lca ']+data['lfe ']+data['lfl ']+data['lph '])#*erg/second
yt.add_field('Cool',function=_Cool,units="erg",sampling_type=st)

def _Upar(field,data):
    return cj21*ph54/data['ndens']
yt.add_field('Upar',function=_Upar,units="cm**3",sampling_type=st)

def _Cool(field,data):
    return (data['lh  ']+data['lhe ']+data['lc  ']+data['ln  ']+data['lo  ']+data['lne ']+data['lna ']+data['lmg ']+data['lsi ']+data['ls  ']+data['lca ']+data['lfe ']+data['lfl ']+data['lph ']+data['lct '])*data['cell_mass']*erg/(second*gram)
yt.add_field('Cool',function=_Cool,units="erg/s",sampling_type=st)

def _CoolH(field,data):
    return data['lh  ']*data['cell_mass']*erg/(second*gram)
yt.add_field('CoolH',function=_CoolH,units="erg/s",sampling_type=st)

def _CoolHe(field,data):
    return data['lhe ']*data['cell_mass']*erg/(second*gram)
yt.add_field('CoolHe',function=_CoolHe,units="erg/s",sampling_type=st)

def _CoolC(field,data):
    return data['lc  ']*data['cell_mass']*erg/(second*gram)
yt.add_field('CoolC',function=_CoolC,units="erg/s",sampling_type=st)

def _CoolN(field,data):
    return data['ln  ']*data['cell_mass']*erg/(second*gram)
yt.add_field('CoolN',function=_CoolN,units="erg/s",sampling_type=st)

def _CoolO(field,data):
    return data['lo  ']*data['cell_mass']*erg/(second*gram)
yt.add_field('CoolO',function=_CoolO,units="erg/s",sampling_type=st)

def _CoolNe(field,data):
    return data['lne ']*data['cell_mass']*erg/(second*gram)
yt.add_field('CoolNe',function=_CoolNe,units="erg/s",sampling_type=st)

def _CoolNa(field,data):
    return data['lna ']*data['cell_mass']*erg/(second*gram)
yt.add_field('CoolNa',function=_CoolNa,units="erg/s",sampling_type=st)

def _CoolMg(field,data):
    return data['lmg ']*data['cell_mass']*erg/(second*gram)
yt.add_field('CoolMg',function=_CoolMg,units="erg/s",sampling_type=st)

def _CoolSi(field,data):
    return data['lsi ']*data['cell_mass']*erg/(second*gram)
yt.add_field('CoolSi',function=_CoolSi,units="erg/s",sampling_type=st)

def _CoolS(field,data):
    return data['ls  ']*data['cell_mass']*erg/(second*gram)
yt.add_field('CoolS',function=_CoolS,units="erg/s",sampling_type=st)

def _CoolCa(field,data):
    return data['lca ']*data['cell_mass']*erg/(second*gram)
yt.add_field('CoolCa',function=_CoolCa,units="erg/s",sampling_type=st)

def _CoolFe(field,data):
    return data['lfe ']*data['cell_mass']*erg/(second*gram)
yt.add_field('CoolFe',function=_CoolFe,units="erg/s",sampling_type=st)

def _CoolFl(field,data):
    return data['lfl ']*data['cell_mass']*erg/(second*gram)
yt.add_field('CoolFl',function=_CoolFl,units="erg/s",sampling_type=st)

def _CoolPh(field,data):
    return data['lph ']*data['cell_mass']*erg/(second*gram)
yt.add_field('CoolPh',function=_CoolPh,units="erg/s",sampling_type=st)

def _CoolCt(field,data):
    return data['lct ']*data['cell_mass']*erg/(second*gram)
yt.add_field('CoolCt',function=_CoolCt,units="erg/s",sampling_type=st)

def _E_tot(field,data): #the E_tot in the sim is specific, needs to be multiplied by cell mass -> erg
    return data['ener']*data['cell_mass']
yt.add_field('E_tot',function=_E_tot,units="erg",sampling_type=st)

def _kinetic(field,data):
    v2 = data['velx']**2 + data['vely']**2 + data['velz']**2
    v = 0.5*data['cell_mass']*v2
    return v
yt.add_field("kinetic",function=_kinetic,units="erg",sampling_type=st)

def _internal(field,data): #the internal in the sim is specific, needs to be multiplied by cell mass -> erg
    return data['eint']*data['cell_mass']
yt.add_field('internal',function=_internal,units="erg",sampling_type=st)

#result = -G*integrate.quad(lambda x: mass**2/x**2,0,r200(mass))[0]
#result2 = -(G*mass/r200(mass)+4*np.pi*G*integrate.quad(lambda x: x*rho(mass,x),0,r200(mass))[0])

def _potential(field,data):
    rad = np.sqrt((data['x'].in_units('cm'))**2+(data['y'].in_units('cm'))**2+(data['z'].in_units('cm'))**2)
    R_s = 6.608e23*cm/10.0
    rho_0 = mass*g/(4*np.pi*R_s**3*(np.log(11)-10.0/11.0)) #3.667827169075327e-25 DM central density
    ans = 4*np.pi*G*rho_0*(R_s)**3/rad*np.log(1+rad/R_s)*data['density']*data['dx']*data['dy']*data['dz']
    return ans
yt.add_field('potential',function=_potential,units="erg",sampling_type=st)

def _magnetic(field,data): #the internal in the sim is specific, needs to be multiplied by cell mass -> erg
    return np.sqrt((data['magx'].in_units('G'))**2+(data['magy'].in_units('G'))**2+(data['magz'].in_units('G'))**2)
yt.add_field('magnetic',function=_magnetic,units="uG",sampling_type=st)

def _mag_energy(field,data): 
    return ((data['magx'].in_units('G'))**2+(data['magy'].in_units('G'))**2+(data['magz'].in_units('G'))**2)*data['dx'].in_units('cm')*data['dy'].in_units('cm')*data['dz'].in_units('cm')/(8.0*np.pi)
yt.add_field('mag_energy',function=_mag_energy,units="erg",sampling_type=st)

def _mag_density(field,data): 
    return (data['magx']**2+data['magy']**2+data['magz']**2)*gauss**2/(8.0*np.pi)*erg/(gauss*cm**3)
yt.add_field('mag_density',function=_mag_density,units="erg/cm**3",sampling_type=st)

def _plasma_beta(field,data): 
    return 8.0*np.pi*data['pressure']/(data['magx']**2+data['magy']**2+data['magz']**2)
yt.add_field('plasma_beta',function=_plasma_beta,sampling_type=st)

def _mag_z(field,data): #getting the magnitude
    return np.sqrt((data['magz'].in_units('G'))**2)
yt.add_field('mag_z',function=_mag_z,units="uG",sampling_type=st)


ylabels = ['N_{O VI} \ (cm^{-2})','N_{N V} \ (cm^{-2})','N_{Si IV} \ (cm^{-2})','N_{Si III} \ (cm^{-2})','N_{Si II} \ (cm^{-2})','N_{H I} \ (cm^{-2})','N_{C IV} \ (cm^{-2})','N_{Mg II} \ (cm^{-2})','N_{C II} \ (cm^{-2})','N_{C III} \ (cm^{-2})','N_{Fe II} \ (cm^{-2})','N_{Fe III} \ (cm^{-2})','N_{N II} \ (cm^{-2})','N_{O I} \ (cm^{-2})','Temperature\ (K)','Number\ Density\ (cm^{-3})','Entropy\ (cm^{2} \dot keV)','Pressure\ (dyn \dot cm^{-2})','Magnetic\ (\mu G)']
#OVI(0), NV(1), SiIV(2), SiIII(3), SiII(4), HI(5), CIV(6), MgII(7), CII(8), CIII(9), FeII(10), FeIII(11), NII(12), OI(13), temp(14), ndens(15),

path = "/Users/mariahjones/Desktop/Research/CGM/data" #my path name
for p in range(1):
    files = [53, 63, 74, 85, 96]
    for i in range(5):
        i = files[i]
        if i < 10:
            stuff = yt.load(path+"/ISM_hdf5_chk_000"+str(i))
        elif 10 <= i < 100:
            stuff = yt.load(path+"/ISM_hdf5_chk_00"+str(i))
        else:
            stuff = yt.load(path+"/ISM_hdf5_chk_0"+str(i))
        sim = stuff.sphere([0.5,0.5,0.5],(214,"kpc"))#all_data() 
        #Various fields
        #O5(0), N4(1), Si3(2), Si2(3), Si(4), H(5), C3(6), Mg(7), C(8), C2(9), 
        #Fe(10), Fe2(11), N(12), O(13), temp(14), #dens(15), pressure(16), 
        #magnetic(17), internal(18), plasma_beta(19), mag_energy(20), 
        #mag_density(21), velz(22)
        field =['o5pdens','n4pdens','si3pdens','si2pdens','sipdens','hdens','c3pdens','mgpdens','cpdens','c2pdens','fepdens','fe2pdens','npdens','odens','temperature','ndens','pressure','magnetic','internal','plasma_beta','mag_energy','mag_density','velz','cell_mass'] 
        # Limits for the aformentioned fields. 
        # If plotting elements, these limits are for projected column densities (ProjectionPlot)
        lim1 = [2e16,2e16,2e16,2e16,2e16,2e20,2e16,2e16,2e16,2e16,2e16,2e16,2e16,2e16,2e7,1e0,2e2,5e-13,2e1,1e55] #
        lim2 = [8e12,8e12,8e12,8e12,8e12,8e12,8e12,8e12,8e12,8e12,8e12,8e12,8e12,8e12,8e3,8e-6,8e-2,5e-15,5e-3,1e51]  #
        for j in range(2): 
            j = 22
           # slc = yt.SlicePlot(stuff,2,field[j],center=[0.5,0.5,0.5],width=(428,'kpc'))
            slc = yt.ProjectionPlot(stuff,1,field[j],center=[0.5,0.5,0.5],width=(428,'kpc'),weight_field=None)
          #  x = 14 ; y = 17 ; k = -1
            #slc = yt.PhasePlot(sim, field[x], field[y], field[k] ,weight_field=None)
            #slc.set_zlim(field[k],lim1[k],lim2[k])
            #slc.annotate_sphere([0.5, 0.5, 0.5], radius=(214, 'kpc'), circle_args={'color':'black'})
            slc.annotate_timestamp(draw_inset_box=True)
            #slc.set_cmap(field=field[k],cmap='YlGnBu') #for projections
            slc.set_cmap(field[j],"YlGnBu") #"RdBu_r")
            slc.set_font({'style':'normal','weight':'bold','size':'32'})
            #slc.set_xlabel(r"$\bf x\ (kpc)$") #if axis=2, x-x,y-y; if axis=1, x-z,y-x; if axis=0, x-y, y-z
            #slc.set_ylabel(r"$\bf y\ (kpc)$")
            #slc.set_colorbar_label(field[j],r"$\rm \bf "+ylabels[j]+"$")
            #slc.annotate_scale(coeff=30,unit='kpc')
            #slc.save() #figure out dpi
            

def phase(names,file_names):
    """
    Function that makes 2D phase plots of temperature and density along with pressure. 
    Original phase plot program (notice no j loops)
    """
    field = ["ndens","temperature","Cool","plasma_beta","cell_mass","Si4_O6","C4_O6","N5_O6","radius","y","ratio","entropy","v_rad","angular_mom","magnetic"]
    axis_font = {'size':'22'}
    for p in range(2):
        p = p+2
        if p==0 : #pretty_gal
            ii = np.arange(0,96,1)
            path = "/Users/mariahjones/Desktop/Research/CGM/data/"
        if p == 1: #high_rot_no_mag_z1
            ii = np.arange(0,101,1)
            path = "/Users/mariahjones/Desktop/Research/CGM/data/"
        if p==2: #high_rot_z1
            ii = [53,96]#np.arange(0,100,1)
            path = "/Users/mariahjones/Desktop/Research/CGM/data/"
       
            ''' 
            if p==3: #Z03 high_rot
                 ii = [76,116]#np.arange(0,118,1)
                 path = "/Volumes/research2/MAIHEM_Z03/"
             if p==4: #Z03 high_rot_no_mag
                 ii = np.arange(0,121,1)
                 path = "/Volumes/research2/MAIHEM_Z03/"
             if p==5: #Z03 high_rot_low_mag
                 ii = np.arange(0,103,1)
                 path = "/Volumes/research2/MAIHEM_Z03/"
            '''
        for i in range(len(ii)):
            i=i + 0
            if ii[i]< 10:
                stuff = yt.load(path+names[p]+"/ISM_hdf5_chk_000"+str(ii[i]))
            elif 10 < ii[i] < 100:
                stuff = yt.load(path+names[p]+"/ISM_hdf5_chk_00"+str(ii[i]))
            elif 100 <= ii[i]:
                stuff = yt.load(path+names[p]+"/ISM_hdf5_chk_0"+str(ii[i]))
            sim = stuff.sphere([0.5,0.5,0.5],(225,"kpc"))#all_data()
            #sim = simulation.cut_region(["obj['ratio'] != np.inf"])
            #ndens [0], temp[1], entropy[2], pressure[3], cell mass[4], Si4/O6[5], C4/O6[6], N5/O6[7], radius[8], y[9], ratio[10], entropy[11], v rad[12], ang y[13]
            x = 1; y = 5; k = 4
            plot = yt.PhasePlot(sim,field[x],field[y],field[k],weight_field=None)
            #plot.set_xlim(10**(-3),5e2)#(4.e-3,1.e1)
            #plot.set_log(field[y],False)
            plot.set_unit(field[k],"Msun")
            #plot.set_ylim(-225,225)#(1e-3,1e3) entropy
            #plot.set_log(field[x],False)
            plot.set_xlim(5e3,1e7) #(2e-5,500)#(8.e-30,2e-23) density
            plot.set_ylim(1e32,1e41)#(2e-18,4e-11)#
            plot.set_cmap(field=field[k],cmap='Spectral')
            plot.set_zlim(field[k],8e3,6e8)#temp,8e3,2e7)#pressure,5e-18,4e-11)#mass,8e3,6e8)
            plot.set_colorbar_label(field[k],r"$\rm Plasma\ \beta$")
            #plot.set_xlabel(r"$\rm \sqrt{X^2 + Z^2}\ (kpc)$")
            plot.set_font({'style':'normal','weight':'bold','size':'28'})
            plot.set_figure_size(10)
            plot.save("Phase_"+file_names[p]+"_"+field[x]+"_"+field[y]+"_"+field[k]+"_"+str(ii[i])+".png")
    return