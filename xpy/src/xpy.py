import math as m
import cmath as mc
import fractions
import random

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
        assert isinstance(self.form, FormFactor) , 'form needs to be a FormFactor'
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
        stri += "damping: " + str(self.damping) + '\n'
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
                    pos[0] = float(words[1])
                    count += 1
                elif words[0] == "y": 
                    pos[1] = float(words[1])
                    count += 1
                elif words[0] == "z": 
                    pos[2] = float(words[1])
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
    
    def save_file(self, path):
        """saves the structure to the file path"""
        f = open(path, 'w')
        f.write('a ' + str(self.lattice_parameters[0]) + '\n')
        f.write('b ' + str(self.lattice_parameters[1]) + '\n')
        f.write('c ' + str(self.lattice_parameters[2]) + '\n')
        for atom in self.atoms:
            f.write('element ' + str(atom.form.atom) + '\n')
            f.write('x ' + str(atom.pos[0]) + '\n')
            f.write('y ' + str(atom.pos[1]) + '\n')
            f.write('z ' + str(atom.pos[2]) + '\n')
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
                multi += (1 - atom.pos[i])*self.lattice_parameters[i]*q[i] 
            structure_factor +=(atom.form.get_value(mq) * 
                                 mc.exp(-1.0j * multi - 
                                        (self.damping / sinomega *
                                         self.lattice_parameters[2] * 
                                          (1 - atom.pos[2])
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
          
           crystal CrystalStructure of the film 
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
        
           crystal -> Crystalstructure 
        """
        assert isinstance(self.crystal, CrystalStructure) , 'crystal needs to be a SimpleCrystalstructure'
        
    _Layer__type_check__ = __type_check__ #private copy of the function to avoid overload
    
    def _set_icq_(self, q, sinomega):
        """recalculates the values for i*(latticeparameters*q)"""                   
        a = self.crystal.lattice_parameters[0]
        b = self.crystal.lattice_parameters[1]
        c = self.crystal.lattice_parameters[2]
        factor = self.crystal.damping * c / sinomega
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
        """returns the phase shift of the layer
        
            q: wave vector
           sinomega: sin(angle between incident beam and sample surface
        """
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
        self._counter_ = 0
        
    def _check_q_s_(self, q, sinomega):
        """checks if q or sinomega are changed.
        
           q: scattering vector
           sinomega: sin(angle between incident beam and sample surface""" 
        values = [q[:], sinomega]
        self._counter_ += 1
        try:
            if (values == self._data_[self._counter_]): #exception if counter bigger than list
                pass
            elif (values == self._data_[self._counter_-1]):
                self._counter_ -= 1
            else:                
                self._counter_ = self._data_.index(values) #exception if values not in data
        except:
            self._data_.append([q[:], sinomega])
            self._reflection_.append(self._calc_reflection_(q, sinomega))
            self._phase_shift_.append(self._calc_phase_shift_(q, sinomega))
            self._counter_ = (len(self._data_) - 2) #makes sure that every data set gets checked
            return self._counter_ + 1
        return self._counter_
        
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

           q: wave vector
           sinomega: sin(angle between incident beam and sample surface
        """
        return ((self.film1.get_phase_shift(q, sinomega) * 
                self.film2.get_phase_shift(q, sinomega)) ** self.bilayers
               )
    
    def get_reflection(self, q, sinomega):
        """returns the complex reflection value of the Superlattice.
             
           q: wave vector
           sinomega: sin(angle between incident beam and sample surface
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
        """returns the Bilayer thickness of the superlatice"""
        return (self.film1.get_thickness() + self.film2.get_thickness())
    
    def get_thickness(self):
        """returns the total thickness of the superlatice"""
        return self.get_Lambda()*self.bilayers

##########
# Sample #
##########

class Sample(object):
        
    def __init__(self, substrate = None, electrode = None, film = None, layers = None):
        """initializes the class data and checks data types
          
           Substrate: substrate of the Sample 
           Electrode: Electrode of the Sample
           Film: Film of the Sample
           Layers: any number of Layers to add to the sample
        """
        
        self._Layers_ = []
        if (substrate):
            self._Layers_.append(substrate)
        if (electrode):
            self._Layers_.append(electrode)
        if (film):
            self._Layers_.append(film)
        if (layers):
            for layer in layers:
                self._Layers_.append(layer)
        self._Sample__type_check__()
        
    def __type_check__(self):
        """checks the class data to have the right types
        
           _Layers_ -> list 
           _Layers_[i] -> Layer
        """
        assert isinstance(self._Layers_, list) , '_Layers_ needs to be an list'
        for i in range(len(self._Layers_)):
            assert (isinstance(self._Layers_[i], (Layer, SuperLattice, SuperLatticeComplex, Sample))) , '_Layers_[' + str(i) + '] needs to be a Layer'
    
    _Sample__type_check__ = __type_check__ #private copy of the function to avoid overload
       
    def get_phase_shift(self,  q, sinomega):
        """returns the phase of the Superlattice.

           q: scattering vector
           sinomega: sin of angle between surface and incedent beam
        """
        phase = 0.0 + 0.0j
        for layer in self._Layers_:
            phase += layer.get_phase_shift(q, sinomega)
        return phase
    
    def get_reflection(self, q, sinomega):
        """returns the complex reflection value of the Superlattice.
             
           q: wave vector
           sinomega: sin of angle between surface and incedent beam
        """
        r = 0.0 + 0.0j
        for i in range(len(self._Layers_)):
            temp = self._Layers_[i].get_reflection(q, sinomega)
            for j in range(i+1, len(self._Layers_)):
                temp *= self._Layers_[j].get_phase_shift(q, sinomega)
            r += temp 
        return r
    
    def get_thickness(self):
        """returns the thickness of the sample."""
        thickness = 0.0 + 0.0j
        for layer in self._Layers_:
            thickness += layer.get_thickness()
        return thickness
    
    def add_Layer(self, layer):
        """adds the layer to the top of the sample."""
        assert (isinstance(layer, (Layer, SuperLattice))) , 'layer needs to be a Layer'
        self._Layers_.append(layer)

#######################
# SuperLatticeComplex #
#######################

class SuperLatticeComplex(Sample):
    
    def __init__(self, crystals = None, layers = None, repetitions=1):
        """initializes the class data and checks data types
          
           crystals: Crystal Structures of the superlattice 
           layers:  number of layers of each CrystalStructure
           repetions: number of repetions          
        """
        self.crystals = crystals
        if not self.crystals:
            self.crystals = []
        self.layers = layers
        if not self.layers:
            self.layers = []
        self.repetitions = repetitions
        self._films_ = {}
        self._Layers_ = []
        
        self._SuperLatticeComplex__type_check__()
        
        self._create_structures_()
        
    def __type_check__(self):
        """checks the class data to have the right types
           crystals -> list
           crystals[] -> CrystalStructure
           layers -> list
           layers[] -> float, int
           repetitions -> float, int

        """
        assert isinstance(self.repetitions, (float, int)) , 'repetitions needs to be a number'
        assert (isinstance(self.crystals, list)) , 'atoms needs to be a list ' 
        for i in range(len(self.crystals)):
            assert (isinstance(self.crystals[i], CrystalStructure)) , 'crystals[' + str(i) + '] needs to be a CrystalStructure'
        assert (isinstance(self.layers, list)) , 'layers needs to be a list ' 
        for i in range(len(self.layers)):
            assert (isinstance(self.layers[i], (float, int))) , 'layers[' + str(i) + '] needs to be a number'
        
    _SuperLatticeComplex__type_check__ = __type_check__ #private copy of the function to avoid overload
    
    def _create_structures_(self):
        """creates the crystal structures needed"""
        leftover = 0.0
        layers = sum(self.layers)*self.repetitions
        i = 0
        while(layers > 0):
            layer = self.layers[i]
            while(layer >= 1):
                if leftover != 0.0:
                    j2 = i-1
                    if j2 < 0:
                        j2 = len(self.layers)-1
                    key = '{}-{}_{}-{}'.format(i, round(1 - leftover,8), j2, leftover)
                    if key not in self._films_.keys():
                        crystal = create_layer(self.crystals[j2], self.crystals[i], leftover)
                        self._films_[key] = LayerSave(crystal)
                    self._Layers_.append(self._films_[key])
                    layer -= 1-leftover
                    leftover = 0.0
                else:
                    key = '{}-{}'.format(i, 1)
                    if key not in self._films_.keys():
                        self._films_[key] = LayerSave(self.crystals[i])
                    self._Layers_.append(self._films_[key])
                    layer -= 1
                
                layers -= 1
                if (layers < 0):
                    break
            leftover = round(layer, 8)
            i += 1
            if i > len(self.layers)-1:
                i = 0

    def get_n(self):
        """returns the number of layers per repetitions"""
        return sum(self.layers)
        
    def get_c(self):
        """returns the average c lattice parameter"""
        return self.get_Lambda() / self.get_n() 
    
    def get_Lambda(self):
        """returns the thicknes of one repetion"""
        l=0.0
        n = 0
        for crystal in self.crystals:
            l += crystal.lattice_parameters[2] * self.layers[n]
            n += 1
        return l
    
    def get_thickness(self):
        """returns the total thickness of the superlatice"""
        return self.get_Lambda()*self.repetitions

##################
# Help Functions #
##################

def l_scan(sample, base_struc, l_min = 0.8, l_max = 1.05, l_step = 0.0002,  h = 0, k = 0, sinomegain = 0, direct = 2.0E8, background = 1, wavelength = 1.5409E-10):

    q = [0.0,0.0,0.0]
    steps_l = (l_max-l_min)/l_step
    q[0] = h*2*m.pi/base_struc.lattice_parameters[0]
    q[1] = k*2*m.pi/base_struc.lattice_parameters[1]
    scan = []
    for i in range(int(steps_l)):
        l=l_min+i*l_step
        q[2] = l*2*m.pi/base_struc.lattice_parameters[2]
        absq = sum([x**2 for x in q])
        if (sinomegain == 0):
            sinomega = absq * wavelength / 4 / m.pi
        else:
            sinomega = sinomegain
        lin = abs(sample.get_reflection(q, sinomega))
        lin *= lin
        lin *= direct
        lin += background
        scan.append([h, k, l, lin])
    
    return scan

def lcm(a,b): 
    return abs(a * b) / fractions.gcd(a,b) if a and b else 0 

def get_multi(x):
    S = str(x).split(".")
    if (len(S) == 1):
        return 1
    L = 10**len(S[1])
    div = alldividers(L)
    for i in div:
        if(i*x == float(int(i*x))):
            return i
                      
def alldividers(number):
    multi = 10.0
    try: #only changes the multiplicator if number is a decimal number
        temp = str(number).split(".")
        multi = 10.0**len(temp[1])
    except:
        pass        
    result = []
    for i in range (1, number + 1):
        if number * multi / i % multi == 0:
            result.append(i)
    return result

def create_layer(struc1, struc2, amount1):
    L = get_multi(amount1) 
    a = max(struc1.lattice_parameters[0],struc2.lattice_parameters[0])
    b = max(struc1.lattice_parameters[1],struc2.lattice_parameters[1])
    c = max(struc1.lattice_parameters[2],struc2.lattice_parameters[2])
    crystal = CrystalStructure([a*L,b,c])
    for j in range(L):
        if (j < amount1*L):
            struc = struc1                
        else:
            struc = struc2 
        for atom in struc.atoms:
            addatom = Atom(atom.form, atom.pos[:] )
            addatom.pos[0] += j*a
            crystal.add_atom(addatom)
    return crystal

def create_structure(a_in, c_in, a11 ,a12 ,a21 ,a22, l1, l2, bl, size = 0.0, zoffset1 = 0.0, zoffset2 = 0.0):
    random.seed()
    L1 = get_multi(l1)
    L2 = get_multi(l2)
    L = lcm(L1,L2)
    layers = bl*(l1+l2)
    if (float(int(layers)) == layers):
        layers = int(layers)
    else:
        layers = int(layers)+1
    (a, b, c) = (a_in*L*10**(-10), a_in*10**(-10) ,c_in*layers*10**(-10) )
    crystal = CrystalStructureCheck()
    crystal.lattice_parameters[0] = a
    crystal.lattice_parameters[1] = b
    crystal.lattice_parameters[2] = c
    
    atom_kind = 1
    l = l1
    offset = 0
    zoffset = zoffset1
    for i in range(layers):
        for j in range(L):
            #check that the right atom is printed
            if (i*L+j - offset >= l*L):
                if (atom_kind == 1):
                    l = l2
                    atom_kind = 2
                    zoffset += zoffset2
                else:
                    l = l1
                    atom_kind = 1
                offset = i*L+j
                zoffset += zoffset1
            #create atoms
            #corner atom
            if (atom_kind == 1):
                atomtype = a11
            else:
                atomtype = a21
            form = FormFactor(path=scriptpath+"/atoms/" + atomtype + ".at")
            pos = [0.0,0.0,0.0]
            pos[0] =(j + 0.0 + random.randrange(-1,1) * size) / L
            pos[1] =(    0.0 + random.randrange(-1,1) * size) 
            pos[2] =(i + 1.0 + zoffset + random.randrange(-1,1) * size) / layers
            atom = Atom(form, pos)
            crystal.add_atom(atom)
            #center atom
            if (atom_kind == 1):
                atomtype = a12
            else:
                atomtype = a22
            form = FormFactor(path=scriptpath+"/atoms/" + atomtype + ".at")
            pos = [0.0,0.0,0.0]
            pos[0] =(j + 0.5 + random.randrange(-1,1) * size) / L
            pos[1] =(    0.5 + random.randrange(-1,1) * size)
            pos[2] =(i + 0.5 + zoffset + random.randrange(-1,1) * size) / layers
            atom = Atom(form, pos)
            crystal.add_atom(atom)
            #oxygen 1
            form = FormFactor(path=scriptpath+"/atoms/O.at")
            pos = [0.0,0.0,0.0]
            pos[0] =(j + 0.5 + random.randrange(-1,1) * size) / L
            pos[1] =(    0.5 + random.randrange(-1,1) * size)
            pos[2] =(i + 1.0 + zoffset + random.randrange(-1,1) * size) / layers
            atom = Atom(form, pos)
            crystal.add_atom(atom)
            #oxygen 2
            form = FormFactor(path=scriptpath+"/atoms/O.at")
            pos = [0.0,0.0,0.0]
            pos[0] =(j + 0.5 + random.randrange(-1,1) * size) / L
            pos[1] =(    0.0 + random.randrange(-1,1) * size)
            pos[2] =(i + 0.5 + zoffset + random.randrange(-1,1) * size) / layers
            atom = Atom(form, pos)
            crystal.add_atom(atom)
            #oxygen 3
            form = FormFactor(path=scriptpath+"/atoms/O.at")
            pos = [0.0,0.0,0.0]
            pos[0] =(j + 0.0 + random.randrange(-1,1) * size) / L
            pos[1] =(    0.5 + random.randrange(-1,1) * size)
            pos[2] =(i + 0.5 + zoffset + random.randrange(-1,1) * size) / layers
            atom = Atom(form, pos)
            crystal.add_atom(atom)
    return crystal

def showup(l1, l2):
    L1 = get_multi(l1)
    L2 = get_multi(l2)
    L = lcm(L1,L2)
    bl = get_multi((l1+l2))
    layers = float(bl)*(l1+l2)
    
    if (float(int(layers)) == layers):
        layers = int(layers)  
    
    atom = 1
    l = l1
    offset = 0
    result = ''
    for i in range(layers):
        for j in range(L):
            #check that the right atom is printed
            if (i*L+j - offset >= l*L):
                if (atom == 1):
                    l = l2
                    atom = 2
                else:
                    l = l1
                    atom = 1
                offset = i*L+j
            if (atom == 1):
                result += 'B'
            else:
                result += 'S'
        result += '\n'
    each = result.split('\n')
    temp = []
    for i in range(len(each)-1, -1 ,-1):
        temp.append(each[i])
    #for line in temp:
        #print line
        
#######################
# Load Data Functions #
#######################

def load_sim(path):
    data = []
    f = open(path, 'r')
    line=f.readline()
    line=f.readline()
    while (line):
        words = line.split(' ')
        if words:
            if (len(words) > 1):
                data.append([words[2],words[3]])
        line=f.readline()
    f.close()
    return data

def load_data(path):
    data = []
    f = open(path, 'r')
    line=f.readline()
    line=f.readline()
    while (line):
        words = line.split(' ')
        if words:
            if (len(words) > 1):
                data.append([words[0],words[1]])
        line=f.readline()
    f.close()
    return data

def load_scan(path, scan_number):
    data = []
    f = open(path, 'r')
    line = f.readline()
    line = f.readline()
    while (line):
        words = line.split(',')
        if words:
            if (len(words) > 1):
                if (words[0] == str(scan_number)):
                    data.append([float(words[1]),float(words[3])])
        line=f.readline()
    f.close()
    return data

def import_scan(path, row1 = 1, row2 = 3):
    data = []
    data.append([])
    f = open(path, 'r')
    line = f.readline()
    line = f.readline()
    while (line):
        words = line.split(',')
        if words:
            if (len(words) > 1):
                try:
                    data[int(words[0])].append([float(words[row1]),float(words[row2])])  
                except:
                    data.append([])
                    data[int(words[0])].append([float(words[row1]),float(words[row2])])
        line=f.readline()
    f.close()
    return data