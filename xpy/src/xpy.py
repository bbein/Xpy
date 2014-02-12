import math as m
import cmath as mc
#import numpy as np

words = __file__.split("/")
scriptpath = ""
for i in range(len(words)):
    if ( i == len(words)-3): break
    scriptpath += words[i]+"/"

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

##############
# FormFactor #
##############

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

####################
# CrystalStructure #
####################

class CrystalStructure(object):
    
    def __init__(self, lattice_parameters = None, atoms = None, damping=5.98E+4, path = None ):
        """initializes the class data and checks data types
          
            lattice_parameters: 3 element list default: [0.0 ...]  
            atoms: list of Atom instances default:[]
            damping: damping coefficient of the structure
            q: scattering vector default: [0.0,0.0,0.0]
            path: path to a file with structure information default: None
        """
        self.damping = damping
        if (atoms):
            self.atoms = atoms
        else:
            self.atoms = []
        if (lattice_parameters):
            self.lattice_parameters = lattice_parameters
        else:
            self.lattice_parameters = [1.0,1.0,1.0]
        self._structure_factor_ = 0.0+0.0j
        self._q_ = [0.0,0.0,0.0]
        self._wavelength_ = 0.0
        #load data from path overwrites given data        
        if path:
            self.load_file(path)
        #check data types
        self._CrystalStructure__type_check__()
             
    def __str__(self):
        """Returns the data in a structured form.
        
        """
        stri = ("Lattice parameters: a=" + str(self.lattice_parameters[0]) +
                " b="+ str(self.lattice_parameters[1]) +
                " c="+ str(self.lattice_parameters[2]) + "\n"
                )
        stri += "damping: " + self.damping
        for atom in self.atoms:
            stri += atom.__str__() + "\n"
        return stri
    
    def __type_check__(self):
        """checks the class data to have the right types
                  
            lattice_parameters -> list
            lattice_parameters[i] -> float
            atoms -> list
            atoms[i] -> Atom
            damping -> float
            _wavelength_ -> float
            _q_ -> list
            _q_[i[ -> float
            _formfactor_ -> complex
        """
        assert isinstance(self.damping, float) , 'damping needs to be a float'
        assert (isinstance(self.lattice_parameters, list) or (len(self.lattice_parameters) != 3)) , 'lattice_parameters needs to be a list with 3 elements' 
        for i in range(len(self.lattice_parameters)):
            assert (isinstance(self.lattice_parameters[i], float)) , 'lattice_parameters[' + str(i) + '] needs to be a number'
        assert (isinstance(self.atoms, list)) , 'atoms needs to be a list ' 
        for i in range(len(self.atoms)):
            assert (isinstance(self.atoms[i], Atom)) , 'atoms[' + str(i) + '] needs to be a Atom'
    
    _CrystalStructure__type_check__ = __type_check__ #private copy of the function to avoid overload
            
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
                if words[0]   == "element": self._add_atom_(f, words[1])
                elif words[0] == "a": self.lattice_parameters[0] = float(words[1])
                elif words[0] == "b": self.lattice_parameters[1] = float(words[1])
                elif words[0] == "c": self.lattice_parameters[2] = float(words[1])
                elif words[0] == "damping": self.damping = float(words[1])
                
            line = f.readline()
        f.close()
        
    def _calc_structure_factor_(self, q, sinomega):
        """calculates the structure factor for the structure 
        
            q: scattering vector
            sinomega: sin(angle between incident beam and sample surface
        """
        if (self.atoms == []):
            return 0.0
        structure_factor = 0.0 + 0.0j
        mq = m.sqrt(sum([x**2 for x in q]))
        for atom in self.atoms:
            multi = 0.0
            for i in range(len(q)):
                multi += atom.pos[i]*q[i] 
            structure_factor +=(atom.form.get_value(mq) * 
                                 mc.exp(-1.0j * multi - 
                                        (self.damping * sinomega *
                                         (self.lattice_parameters[2] - 
                                          atom.pos[2]
                                         ) 
                                        ) 
                                       )     
                                )
        return structure_factor
                  
    def get_structure_factor(self, q, sinomega):
        """calculates the structure factor for the structure 
        
            q: scattering vector
            sinomega: sin(angle between incident beam and sample surface
        """
        return self._calc_structure_factor_(q, sinomega)

#########################
# CrystalStructureCheck #
#########################

class CrystalStructureCheck(CrystalStructure):
    
    def get_structure_factor(self, q, sinomega):
        """calculates the structure factor for the structure 
        
            q: scattering vector
            sinomega: sin(angle between incident beam and sample surface
        """
        #assert isinstance(q,self._q_.__class__), 'q needs to be ' + str(self._q_.__class__)
        if (self._q_ != q or self._sinomega_ != sinomega): 
            self._structure_factor_ = self._calc_structure_factor_(q, sinomega)
            self._q_ = q[:]
            self._sinomega_ = sinomega
        return self._structure_factor_

#########
# Layer #
#########

class Layer(object):
        
    def __init__(self, crystal):
        """initializes the class data and checks data types
          
           crystal SimpleCrystalStructure of the film 
        """
        self.crystal= crystal
        self._icq_ = self._ibq_ = self._iaq_ = 0.0
        self.delta = 0.0
        
        self._Layer__type_check__()
    
    def __str__(self):
        """Returns the data in a structured form.
           
           returns its crystal structure 
        """
        return self.crystal.__str__()
    
    def __type_check__(self):
        """checks the class data to have the right types
        
           crystal -> SimpleCrystalstructure 
        """
        assert isinstance(self.crystal, CrystalStructure) , 'crystal needs to be a SimpleCrystalstructure'
        
    _Layer__type_check__ = __type_check__ #private copy of the function to avoid overload
    
    def _set_icq_(self, q, sinomega):
        """recalculates the values for i*(latticeparameters*q)"""                   
        a = self.crystal.lattice_parameters[0]
        b = self.crystal.lattice_parameters[1]
        c = self.crystal.lattice_parameters[2]
        factor = sinomega * self.crystal.damping * c
        self._iaq_ = -1J * a * q[0] - factor 
        self._ibq_ = -1J * b * q[1] - factor 
        self._icq_ = -1J * c * q[2] - factor
    
    def _calc_amplitude_(self, q):
        """returns calculates the scattering amplitude"""
        a = self.crystal.lattice_parameters[0]
        b = self.crystal.lattice_parameters[1]
        R = 2.1879e-15 #R= 2.1879e-15 is the classical electron radius
        amp1d = 4 * m.pi * R / m.sqrt(sum([x**2 for x in q])) / a / b
        return (1j*(amp1d) * (1 - mc.exp(self._icq_.real))**2)
    
    def _calc_reflection_(self, q, sinomega):
        """calculates the complex reflection value of the Film."""
        self._set_icq_(q, sinomega)
        amp = self._calc_amplitude_(q)
        r = (amp * self.crystal.get_structure_factor(q, sinomega) / 
             ((1 - mc.exp(self._iaq_)) * (1 - mc.exp(self._ibq_)))             
            )
        return r
    
    def _calc_phase_shift_(self, q, sinomega):
        """returns the phase shift of the layer"""
        return mc.exp( (1 + self.delta) * self._icq_)
        
    def get_reflection(self, q, sinomega):
        """returns the complex reflection value of the Film.
             
           q: wave vector
           sinomega: sin(angle between incident beam and sample surface
        """
        return self._calc_reflection_(q, sinomega)
    
    def get_phase_shift(self, q, sinomega):
        """returns the phase of the film.

           q: scattering vector
           sinomega: sin angle between surface plane and incedent beam 
        """
        return self._calc_phase_shift_(q, sinomega)
    
    def get_thickness(self):
        """returns the thickness of the film."""
        return (1 + self.delta)*self.crystal.lattice_parameters[2]

##############
# LayerCheck #
##############

class LayerCheck(Layer):
    
    def __init__(self, crystal):
        Layer.__init__(self, crystal)
        self._q_ = [0,0,0]
        self._sinomega_ = 0.0
        self._reflection_ = 0.0 + 0.0j
        self._phase_shift_ = 0.0 + 0.0j
    
    def get_reflection(self, q, sinomega):
        """returns the complex reflection value of the Film.
             
           q: wave vector
           sinomega: sin(angle between incident beam and sample surface
        """
        self._check_q_s_(q, sinomega)
        return self._reflection_
    
    def get_phase_shift(self, q, sinomega):
        """returns the phase of the film.

           q: scattering vector
           sinomega: sin angle between surface plane and incedent beam 
        """
        self._check_q_s_(q, sinomega)
        return self._phase_shift_
    
    def _check_q_s_(self, q, sinomega):
        """checks if q or sinomega are changed.
        
           q: scattering vector
           sinomega: sin(angle between incident beam and sample surface 
        """
        #assert isinstance(q,self._q_.__class__), 'q needs to be ' + str(self._q_.__class__)
        #assert isinstance(sinomega,self._sinomega_.__class__), 'sinomega needs to be ' + str(self._sinomega_.__class__)
        if ( (self._q_ != q) or 
             (self._sinomega_ != sinomega) ):
            self._q_ = q[:]
            self._sinomega_ = sinomega
            self._reflection_ = self._calc_reflection_(q, sinomega)
            self._phase_shift_ = self._calc_phase_shift_(q, sinomega)
            return True
        return False

#############
# LayerSave #
#############

class LayerSave(Layer):
    
    def __init__(self, crystal):
        Layer.__init__(self, crystal)
        self._data_ = []
        self._reflection_ = []
        self._phase_shift_ = []
        
    def _check_q_s_(self, q, sinomega):
        """checks if q or sinomega are changed.
        
           q: scattering vector
           sinomega: sin(angle between incident beam and sample surface""" 
        values = [q[:], sinomega]
        try:
            return self._data_.index(values)
        except:
            self._data_.append([q[:], sinomega])
            self._reflection_.append(self._calc_reflection_(q, sinomega))
            self._phase_shift_.append(self._calc_phase_shift_(q, sinomega))
            return  (len(self._data_) - 1)
    
    def get_reflection(self, q, sinomega):
        """returns the complex reflection value of the Film.
             
           q: wave vector
           sinomega: sin(angle between incident beam and sample surface
        """
        return self._reflection_[self._check_q_s_(q, sinomega)]
    
    def get_phase_shift(self, q, sinomega):
        """returns the phase of the film.

           q: scattering vector
           sinomega: sin angle between surface plane and incedent beam 
        """
        return self._phase_shift_[self._check_q_s_(q, sinomega)]
           

#############
# ThickFilm #
#############
 
class ThickFilm (Layer):

    def _calc_reflection_(self, q, sinomega):
        """calculates the complex reflection value of the Film."""
        r = Layer._calc_reflection_(self, q, sinomega)
        r *= 1 / (1 - mc.exp(self._icq_))
        r = 2 * r / (1 + m.sqrt(1 + 4 * abs(r**2)))
        return r

##################
# ThickFilmCheck #
##################

class ThickFilmCheck (LayerCheck, ThickFilm):
    pass

#################
# ThickFilmSave #
#################

class ThickFilmSave (LayerSave, ThickFilm):
    pass

############
# ThinFilm #
############

class ThinFilm(Layer):

    def __init__(self, crystal, layers=1):
        """initializes the class data and checks data types
          
           crystal: SimpleCrystalStructure of the film
           layers:  number of layers in the film
        """
        super(ThinFilm, self).__init__(crystal)
        self._layers_ = layers
        self._icnq_ = 0.0+0.0j   
        self._ThinFilm__type_check__()
        
    def __type_check__(self):
        """checks the class data to have the right types
        
           layers -> int 
        """
        assert isinstance(self._layers_, (int)) , 'layers needs to be a int'
        
    _ThinFilm__type_check__ = __type_check__ #private copy of the function to avoid overload
    
    def _set_icq_(self,  q, sinomega):
        """recalculates the values for i*(latticeparameters*q*n)"""
        super(ThinFilm, self)._set_icq_(q, sinomega)
        self._icnq_ =  self._layers_ * self._icq_
    
    def _calc_reflection_(self, q, sinomega):
        """calculates the reflection and phase shift of the Film."""
        r = Layer._calc_reflection_(self, q, sinomega)
        r *= (1-mc.exp(self._icnq_)) / (1 - mc.exp(self._icq_))
        return r      

    def _calc_phase_shift_(self, q, sinomega):
        """returns the phase shift of the layer"""
        return mc.exp( (self._layers_ + self.delta) * self._icq_)
    
    def get_thickness(self):
        """returns the thickness of the film."""
        return (self._layers_ + self.delta)*self.crystal.lattice_parameters[2]

#################
# ThinFilmCheck #
#################

class ThinFilmCheck(LayerCheck, ThinFilm):    
    
    def __init__(self, crystal, layers=1):
        ThinFilm.__init__(self, crystal, layers)
        LayerCheck.__init__(self, crystal)
        
################
# ThinFilmSave #
################

class ThinFilmSave(LayerSave, ThinFilm):    
    
    def __init__(self, crystal, layers=1):
        ThinFilm.__init__(self, crystal, layers)
        LayerSave.__init__(self, crystal)  

################
# SuperLattice #
################

class SuperLattice(object):
        
    def __init__(self, film1, film2, bilayers=1):
        """initializes the class data and checks data types
          
           film1: first film of the superlattice 
           film2: second film of the superlattice 
           bilayers: number of times both films get repeted
        """
        
        self.film1 = film1
        self.film2 = film2
        self.bilayers = bilayers
        self._SuperLattice__type_check__()
        
    def __type_check__(self):
        """checks the class data to have the right types
        
           bilayers -> float 
           film1 -> Layer
           film2 -> Layer
        """
        assert isinstance(self.bilayers, (int)) , 'bilayers needs to be an int'
        assert isinstance(self.film1, Layer) , 'film1 needs to be a Layer'
        assert isinstance(self.film2, Layer) , 'film2 needs to be a Layer'
    
    _SuperLattice__type_check__ = __type_check__ #private copy of the function to avoid overload
       
    def get_phase_shift(self,  q, sinomega):
        """returns the phase of the Superlattice.

           q: scattering vector
           wavelength: wavelength 
        """
        return ((self.film1.get_phase_shift(q, sinomega) * 
                self.film2.get_phase_shift(q, sinomega)) ** self.bilayers
               )
    
    def get_reflection(self, q, sinomega):
        """returns the complex reflection value of the Superlattice.
             
           q: wave vector
           wavelength: wavelength
        """
        r1 = self.film1.get_reflection(q, sinomega)
        r2 = self.film2.get_reflection(q, sinomega)
        p1 = self.film1.get_phase_shift(q, sinomega)
        p2 = self.film2.get_phase_shift(q, sinomega)
        r = 0.0 + 0.0j
        for i in range(self.bilayers):
            r += ((r1 * p2 + r2) * ((p1 * p2) ** i))
        return r
    

    def get_n(self):
        """returns the number of layers per bilayer"""
        if (self.film1._layers_): 
            l1 = self.film1._layers_
        else:
            l1 = 1
        if (self.film2._layers_): 
            l2 = self.film2._layers_
        else:
            l2 = 1
        return l1 + l2
    
    def get_c(self):
        """returns the average c lattice parameter"""
        return self.get_Lambda() / self.get_n()
    
    def get_Lambda(self):
        """returns the Total thicknes of the superlatice"""
        return (self.film1.get_thickness() + self.film2.get_thickness())

#######################
# SuperLatticeComplex #
#######################

class SuperLatticeComplex(SuperLattice):
    
    precision = 1
    
    def __init__(self, crystal1, crystal2, layers1=1, layers2=1, bilayers=1, damping=5.98E+4):
        """initializes the class data and checks data types
          
           crystal1: SimpleCrystalStructure of the first film
           layers1:  number of layers in the first film
           crystal2: SimpleCrystalStructure of the second film
           layers2:  number of layers in the second film
           bilayers: number of times both films get repeted
           damping: damping for the film
           
        """
        
        self.crystal1 = crystal1
        self.layers1 = layers1
        self.crystal2 = crystal2
        self.layers2 = layers2
        self.bilayers = bilayers
        self.damping = damping
        self._SuperLatticeComplex__type_check__()
        
    def __type_check__(self):
        """checks the class data to have the right types
           crystal1 -> SimpleCrystalstructure
           layers1 -> float
           crystal2 -> SimpleCrystalstructure
           layers2 -> float
           damping -> float
        """
        assert isinstance(self.crystal1, CrystalStructure) , 'crystal1 needs to be a SimpleCrystalstructure'
        assert isinstance(self.layers1, (float, int)) , 'bilayers needs to be a number'
        assert isinstance(self.crystal2, CrystalStructure) , 'crystal2 needs to be a SimpleCrystalstructure'
        assert isinstance(self.layers2, (float, int)) , 'bilayers needs to be a number'
        assert isinstance(self.damping, float) , 'damping needs to be a float'
        
    _SuperLatticeComplex__type_check__ = __type_check__ #private copy of the function to avoid overload
    

    def get_n(self):
        """returns the number of layers per bilayer"""
        return self.layers1 + self.layers2
        
    def get_c(self):
        """returns the average c lattice parameter"""
        c=((self.crystal1.lattice_parameters[2] * self.layers1 +
            self.crystal2.lattice_parameters[2] * self.layers2) / 
            self.get_n() 
          )
        return c 

    