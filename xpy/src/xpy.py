import math as m
import cmath as mc
#import numpy as np

words = __file__.split("/")
scriptpath = ""
for i in range(len(words)):
    if ( i == len(words)-3): break
    scriptpath += words[i]+"/"

####################
# FormFactor #
####################

class FormFactor(object): 
    
    def __init__(self,atom="",a=None,b=None,q=0.0,value=0.0,path=None):
        """initializes the class data and checks data types
        
            atom:  chemical Symbol                   default: ""
            q:     absolute value of the wave vector default: 0.0   
            value: last calculated value             default: 0.0   
            a: 5 element list of prefactors          default: [0.0 ...]  
            b: 4 element list of exponents           default: [0.0 ...]
        """
        self.q = q
        self.value = value
        self.atom = atom
        if (a):
            self.a = a
        else:
            self.a = [0.0,0.0,0.0,0.0,0.0]
        if(b): 
            self.b = b
        else:
            self.b = [0.0,0.0,0.0,0.0]
        #load data from path overwrites given data        
        if path:
            self.load_file(path)
        #check data types
        self._FormFactor__type_check__()
        
    def __str__(self):
        """Returns the data in a structured form.
        
            returns:
            element: 'atom'
            a1 = 'a[0]; a2 = 'a[1]; a3 = 'a[2]; a4 = 'a[3]; a5 = 'a[4];
            b1 = 'b[0]; b2 = 'b[1]; b3 = 'b[2]; b4 = 'b[3];
        """
        stri = "element:" + self.atom + "\n"
        for i in range(len(self.a)): 
            stri += 'a' + str(i+1) + " = " + str(self.a[i]) + "; "
        stri += "\n"
        for i in range(len(self.b)): 
            stri += 'b' + str(i+1) + " = " + str(self.b[i]) + "; "
        return stri
        
    def __type_check__(self):
        """checks the class data to have the right types.
        
            q -> float
            value -> float
            atom -> str
            a -> list
            a[i] -> float
            b -> list
            b[i] -> float
        """
        assert isinstance(self.q, float) , 'q needs to be a float'
        assert isinstance(self.value, float) , 'value needs to be a float'
        assert isinstance(self.atom, str) , 'atom needs to be a string'
        assert (isinstance(self.a, list) or (len(self.a) != 5)) , 'a needs to be a list with 5 elements'
        for i in range(len(self.a)):
            assert (isinstance(self.a[i], (float, int))) , 'a' + str(i) + 'needs to be a number'
        assert (isinstance(self.b, list) or (len(self.b) != 4)) , 'b needs to be a list with 4 elements' 
        for i in range(len(self.b)):
            assert (isinstance(self.b[i], (float, int))) , 'b' + str(i) + 'needs to be a number'

    _FormFactor__type_check__ = __type_check__ #private copy of the function to avoid overload
    
    def get_value(self, q):
        """returns the structure factor value for the given q
        
           input q: absolute value of the wave vector      
        """
        assert isinstance(q,self.q.__class__), 'q needs to be ' + str(self.q.__class__)
        if (self.q != q): 
            q2 = (q*q) / (16*m.pi*m.pi) * 1e-20 
            self.value =(self.a[0] * m.exp(self.b[0] * q2) + 
                         self.a[1] * m.exp(self.b[1] * q2) + 
                         self.a[2] * m.exp(self.b[2] * q2) + 
                         self.a[3] * m.exp(self.b[3] * q2) + 
                         self.a[4]
                        )
            self.q = q
        return self.value
        
    def load_file(self, path):
        """loads the structure factor the file in path"""
        f = open(path, 'r')
        line = f.readline()
        while (line):
            words = line.split()
            if words:
                if words[0] == "atom": self.atom = words[1]
                elif words[0] == "a1": self.a[0]= float(words[1])
                elif words[0] == "a2": self.a[1]= float(words[1])
                elif words[0] == "a3": self.a[2]= float(words[1])
                elif words[0] == "a4": self.a[3]= float(words[1])
                elif words[0] == "a5": self.a[4]= float(words[1])
                elif words[0] == "b1": self.b[0]= float(words[1])
                elif words[0] == "b2": self.b[1]= float(words[1])
                elif words[0] == "b3": self.b[2]= float(words[1])
                elif words[0] == "b4": self.b[3]= float(words[1])
            line = f.readline()
        f.close()

    
    
########
# Atom #
########

class Atom(object):
    
    def __init__(self, form, pos = None):
        """initializes the class data and checks data types
          
            pos: 3 element list of relative atom position default: [0.0 ...]  
            struc: Structure Factor of the Atom default: N/A
        """
        if (pos):
            self.pos = pos
        else:
            self.pos = [0.0,0.0,0.0] 
        self.form = form
        #check data types
        self._Atom__type_check__()    
        
    def __str__(self):
        """Returns the data in a structured form.
           
           returns:
           atom: atom x=pos[1] y=pos[2] z= pos[3] 
        """
        
        stri = ("atom: " + self.form.atom + " x=" + str(self.pos[0]) +
                " y=" + str(self.pos[1]) + " z=" + str(self.pos[2]) 
                )
        return stri
        
    def __type_check__(self):
        """checks the class data to have the right types
        
            pos    -> list
            pos[i] -> float
            struc  -> StructureFactor
        """
        assert isinstance(self.form, FormFactor) , 'struc needs to be a StructureFactor'
        assert (isinstance(self.pos, list) or (len(self.pos) != 3)) , 'a needs to be a list with 3 elements'
        for i in range(len(self.pos)):
            assert (isinstance(self.pos[i], (float, int))) , 'pos' + str(i) + 'needs to be a number'
            
    _Atom__type_check__ = __type_check__ #private copy of the function to avoid overload

##########################
# SimpleCrystalStructure #
##########################

class SimpleCrystalStructure(object):
    
    def __init__(self, lattice_parameters = None, atoms = None, q = None, path = None ):
        """initializes the class data and checks data types
          
            lattice_parameters: 3 element list default: [0.0 ...]  
            atoms: list of Atom instances default:[]
            q: scattering vector default: [0.0,0.0,0.0]
            path: path to a file with structure information default: None
        """
        if (atoms):
            self.atoms = atoms
        else:
            self.atoms = []
        if (q):
            self._q_ = q
        else:
            self._q_ = [0.0,0.0,0.0]
        if (lattice_parameters):
            self.lattice_parameters = lattice_parameters
        else:
            self.lattice_parameters = [1.0,1.0,1.0]
        self._formfactor_ = 0.0+0.0j
        #load data from path overwrites given data        
        if path:
            self.load_file(path)
        #check data types
        self._SimpleCrystalStructure__type_check__()
        
    def __str__(self):
        """Returns the data in a structured form.
        
        """
        stri = ("Lattice parameters: a=" + str(self.lattice_parameters[0]) +
                " b="+ str(self.lattice_parameters[1]) +
                " c="+ str(self.lattice_parameters[2]) + "\n"
                )
        for atom in self.atoms:
            stri += atom.__str__() + "\n"
        return stri
    
    def __type_check__(self):
        """checks the class data to have the right types
        
            _wavelength_ -> float
            lattice_parameters -> list
            lattice_parameters[i] -> float
            atoms -> list
            atoms[i] -> Atom
            _q_ -> list
            _q_[i[ -> float
            _formfactor_ -> complex
        """
        assert (isinstance(self.lattice_parameters, list) or (len(self.lattice_parameters) != 3)) , 'lattice_parameters needs to be a list with 3 elements' 
        for i in range(len(self.lattice_parameters)):
            assert (isinstance(self.lattice_parameters[i], float)) , 'lattice_parameters[' + str(i) + '] needs to be a number'
        assert (isinstance(self._q_, list) or (len(self._q_) != 3)) , 'q needs to be a list with 3 elements' 
        for i in range(len(self._q_)):
            assert (isinstance(self._q_[i], float)) , 'q[' + str(i) + ']n eeds to be a number'
        assert (isinstance(self.atoms, list)) , 'atoms needs to be a list ' 
        for i in range(len(self.atoms)):
            assert (isinstance(self.atoms[i], Atom)) , 'atoms[' + str(i) + '] needs to be a Atom'
        assert (isinstance(self._formfactor_, complex)) , 'formfactor needs to be a complex number '
    
    _SimpleCrystalStructure__type_check__ = __type_check__ #private copy of the function to avoid overload
            
    def _add_atom_(self, usedfile, atomtype):
        """adds the next atom in the file to the atom list """
        form = FormFactor(path=scriptpath+"/atoms/" + atomtype + ".at")
        pos = [0.0,0.0,0.0];
        line = usedfile.readline()
        count = 0
        while (line):
            words = line.split()
            if words:
                if words[0]   == "x": 
                    pos[0] = float(words[1]) * self.lattice_parameters[0]
                    count += 1
                elif words[0] == "y": 
                    pos[1] = float(words[1]) * self.lattice_parameters[1]
                    count += 1
                elif words[0] == "z": 
                    pos[2] = float(words[1]) * self.lattice_parameters[2]
                    count += 1
                if (count == 3): break 
            line = usedfile.readline()
        newatom = Atom(form, pos)
        self.add_atom(newatom)    

    def add_atom(self, atom):
        """adds atom to the atoms list"""
        assert isinstance(atom, Atom) , 'atom needs to be an Atom'
        for oldatom in self.atoms:
            if(atom.form.atom == oldatom.form.atom): 
                atom.form = oldatom.form
                break
        self.atoms.append(atom)

    def load_file(self, path):
        """loads the structure from the file in path"""
        f = open(path, 'r')
        line = f.readline()
        while (line):
            words = line.split()
            if words:
                if words[0]   == "a": self.lattice_parameters[0] = float(words[1])
                elif words[0] == "b": self.lattice_parameters[1] = float(words[1])
                elif words[0] == "c": self.lattice_parameters[2] = float(words[1])
                elif words[0] == "element": self._add_atom_(f, words[1])
            line = f.readline()
        f.close()
        
    def get_structure_factor(self, q):
        """calculates the structure factor for the structure 
        
            q: scattering vector
        """
        if (self.atoms == []):
            return 0.0
        assert isinstance(q,self._q_.__class__), 'q needs to be ' + str(self._q_.__class__)
        if ((self._q_ != q)): 
            self._formfactor_ = 0.0 + 0.0j
            mq = m.sqrt(sum([x**2 for x in q]))
            for atom in self.atoms:
                for i in range(len(q)):
                    multi = atom.pos[i]*q[i] 
                self._formfactor_ +=(atom.form.get_value(mq) *
                                   mc.exp(-1.0j * multi)
                                   )
            self._q_ = q[:]
        return self._formfactor_   

    
    
########
# Film #
########
 
class Film (object):

    def __init__(self, crystal, damping=5.98E+4):
        """initializes the class data and checks data types
          
           crystal SimpleCrystalStructure of the film
           damping for the film 
        """
        self.crystal= crystal
        self.damping = damping
        self._reflection_ = 0.0;
        self._wavelength_ = 0.0;
        self._icq_ = self._ibq_ = self._iaq_ = 0.0;
        self._amplitude_ = 0.0;
        self._q_ = [0.0, 0.0, 0.0]
        self._Film__type_check__()

    def __str__(self):
        """Returns the data in a structured form.
           
           returns its crystal structure 
        """
        return self.crystal.__str__()
    
    def __type_check__(self):
        """checks the class data to have the right types
        
           crystal -> SimpleCrystalstructure
           damping -> float 
        """
        assert isinstance(self.crystal, SimpleCrystalStructure) , 'crystal needs to be a SimpleCrystalstructure'
        assert isinstance(self.damping, float) , 'damping needs to be a float'
        
    _Film__type_check__ = __type_check__ #private copy of the function to avoid overload

    def _set_icq_(self):
        """recalculates the values for i*(latticeparameters*q)"""
        factor = (4 * m.pi / self._wavelength_ * 
                  self.damping / m.sqrt(sum([x**2 for x in self._q_]))
                  )
        a = self.crystal.lattice_parameters[0]
        b = self.crystal.lattice_parameters[1]
        c = self.crystal.lattice_parameters[2]
        self._iaq_ = (-1J * (a * self._q_[0])) + (-a * factor)   
        self._ibq_ = (-1J * (b * self._q_[1])) + (-b * factor) 
        self._icq_ = (-1J * (c * self._q_[2])) + (-c * factor)
    
    def _calc_amplitude_(self):
        """calculates the scattering amplitude"""
        #where is a^2 comming from?
        a = self.crystal.lattice_parameters[0]
        R = 2.1879e-15 #R= 2.1879e-15 is the classical electron radius
        amp1d = 4 * m.pi * R / m.sqrt(sum([x**2 for x in self._q_])) / a / a
        self._amplitude_ = 0.0 + 1j*(amp1d);

    def _calc_reflection_(self):
        """calculates the complex reflection value of the Film."""
        self._set_icq_();
        self._calc_amplitude_();
        r = (self._amplitude_ * self.crystal.get_structure_factor(self._q_) / 
             (#(1 - mc.exp(self._iaq_)) * (1 - mc.exp(self._ibq_)) * 
             (1 - mc.exp(self._icq_))
            ))
        self._reflection_ = 2 * r / (1 + m.sqrt(1 + 4 * abs(r**2)))
    
    def _check_q_l_(self, q, wavelength):
        """checks if q or lambda are changed.
        
           q: scattering vector
           wavelength: wavelength 
        """
        assert isinstance(q,self._q_.__class__), 'q needs to be ' + str(self._q_.__class__)
        assert isinstance(wavelength,self._wavelength_.__class__), 'wavelength needs to be ' + str(self._wavelength_.__class__)
        if ( (self._q_ != q) or (self._wavelength_ != wavelength) ):
            self._q_ = q[:];
            self._wavelength_ = wavelength;
            return True;
        return False;
    
    def get_reflection(self, q, wavelength):
        """returns the complex reflection value of the Film.
             
           q: wave vector
           wavelength: wavelength
        """
        if (self._check_q_l_(q, wavelength)):
            self._calc_reflection_()
        return self._reflection_;

############
# ThinFilm #
############

class ThinFilm(Film):

    def __init__(self, crystal, layers=1, damping=5.98E+4):
        """initializes the class data and checks data types
          
           crystal: SimpleCrystalStructure of the film
           damping: damping for the film
           layers:  number of layers in the film
        """
        super(ThinFilm, self).__init__(crystal, damping)
        self._layers_ = layers
        self._icnq_ = 0.0+0.0j
        self._phaseshift_ = 0.0+0.0j
        self._ThinFilm__type_check__()

    def __type_check__(self):
        """checks the class data to have the right types
        
           layers -> float 
        """
        assert isinstance(self._layers_, (float, int)) , 'layers needs to be a number'
        
    _ThinFilm__type_check__ = __type_check__ #private copy of the function to avoid overload
    
    def _set_icq_(self):
        """recalculates the values for i*(latticeparameters*q)"""
        super(ThinFilm, self)._set_icq_()
        c = self.crystal.lattice_parameters[2]
        factor = (4 * m.pi / self._wavelength_ * self.damping / 
                 m.sqrt(sum([x**2 for x in self._q_]))
                 )
        self._icnq_ = (-1.0j * self._layers_ * c * self._q_[2] - 
                     self._layers_ * c * factor
                     )  
    
    def _calc_reflection_(self):
        """calculates the reflection and phase shift of the Film."""
        self._set_icq_()
        self._phaseshift_ = mc.exp(self._icnq_)
        self._calc_amplitude_()
        self._reflection_ = (self._amplitude_ * 
                             self.crystal.get_structure_factor(self._q_) * 
                             (1-self._phaseshift_) / (#(1-mc.exp(self._iaq_))* 
                                                       #(1-mc.exp(self._ibq_))* 
                                                       (1-mc.exp(self._icq_))
                                                      )
                            )    

    def get_phase_shift(self,  q, wavelength):
        """returns the phase of the film.

           q: scattering vector
           wavelength: wavelength 
        """
        if (self._check_q_l_, q, wavelength):
            self._calc_reflection_()
        return self._phaseshift_
    
class SuperLattice(object):
    
    def __init__(self, film1, film2, bilayers=1):
        """initializes the class data and checks data types
          
           crystal: SimpleCrystalStructure of the film
           damping: damping for the film
           layers:  number of layers in the film
        """
        
        self.film1 = film1
        self.film2 = film2
        self.bilayers = bilayers
        self._SuperLattice__type_check__()
        
    def __type_check__(self):
        """checks the class data to have the right types
        
           bilayers -> float 
           film1 -> ThinFilm
           film2 -> ThinFilm
        """
        assert isinstance(self.bilayers, (float, int)) , 'bilayers needs to be a number'
        assert isinstance(self.film1, ThinFilm) , 'film1 needs to be a ThinFilm'
        assert isinstance(self.film2, ThinFilm) , 'film2 needs to be a ThinFilm'
    
    _SuperLattice__type_check__ = __type_check__ #private copy of the function to avoid overload
       
    def get_phase_shift(self,  q, wavelength):
        """returns the phase of the Superlattice.

           q: scattering vector
           wavelength: wavelength 
        """
        return ((self.film1.get_phase_shift(q, wavelength) * 
                self.film2.get_phase_shift(q, wavelength)) ** self.bilayers
               )
    
    def get_reflection(self, q, wavelength):
        """returns the complex reflection value of the Superlattice.
             
           q: wave vector
           wavelength: wavelength
        """
        r1 = self.film1.get_reflection(q, wavelength)
        r2 = self.film2.get_reflection(q, wavelength)
        p1 = self.film1.get_phase_shift(q, wavelength)
        p2 = self.film2.get_phase_shift(q, wavelength)
        r = 0.0 + 0.0j
        for i in range(self.bilayers):
            r += ((r1 * p2 + r2) * ((p1 * p2) ** i))
        return r
    

    def get_n(self):
        """returns the number of layers per bilayer"""
        return self.film1._layers_ + self.film2._layers_
    
    
    def get_c(self):
        """returns the average c lattice parameter"""
        c=((self.film1.crystal.lattice_parameters[2] * self.film1._layers_ +
            self.film2.crystal.lattice_parameters[2] * self.film2._layers_) / 
            self.get_n() 
          )
        return c
    
    def get_Lambda(self):
        return self.get_c()*self.get_n()