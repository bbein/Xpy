import math as m
import xpy

l_min = 0.8
l_max = 1.05
l_step = 0.001

wavelength = 1.5409E-10
direct = 1.0E16
background = 1

STO =  xpy.SimpleCrystalStructure(path = xpy.scriptpath+"Perovskites/STO.str")
PTO =  xpy.SimpleCrystalStructure(path = xpy.scriptpath+"Perovskites/PTOs.str")
substrate =  xpy.Film(STO)
film = xpy.ThinFilm(PTO)

q = [0.0,0.0,0.0]
steps_l = (l_max-l_min)/l_step
#Output
data = []
for i in range(int(steps_l)):
    l=l_min+i*l_step
    q[2] = l*2*m.pi/STO.lattice_parameters[2]
    #lin = abs(film.get_reflection(q, wavelength) + substrate.get_reflection(q, wavelength) * film.get_phase_shift(q, wavelength))
    lin = abs(film.get_reflection(q, wavelength))
    #lin = abs(substrate.get_reflection(q, wavelength))
    lin *= lin
    lin *= direct
    lin += background;
    data.append([l, lin])
    print(l, lin, m.log10(lin))