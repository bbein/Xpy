import math as m
import xpy

l_min = 0.8
l_max = 1.05
l_step = 0.0002 #5000 points per order of l

wavelength = 10**-10
direct = 2.0E8
background = 1

STO = xpy.CrystalStructureCheck(path = xpy.scriptpath+"Perovskites/STOh.str")
substrate = xpy.ThickFilmCheck(STO)

SRO = xpy.CrystalStructureCheck(path = xpy.scriptpath+"Perovskites/SROs.str")
electrode = xpy.ThinFilmCheck(SRO,41)
SRO.lattice_parameters[0] = SRO.lattice_parameters[1] = 3.925*10**-10
SRO.lattice_parameters[2] = 3.999*10**-10
electrode.delta = 0.0

substrate = xpy.ThickFilmCheck(STO)
sample = xpy.Sample(substrate)

sub001 = xpy.l_scan(sample, STO, l_min = l_min, l_max = l_max, l_step = l_step, direct = direct, background = background, wavelength = wavelength)