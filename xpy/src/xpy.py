import math as m
import cmath as mc
import fractions
import random

words = __file__.split("/")
_SCRIPTPATH_ = ""
for i in range(len(words)):
    if ( i == len(words)-3): break
    _SCRIPTPATH_ += words[i]+"/"
    
_NUMBER_ = (int, float) # np is not working in 3.3 yet ,np.floating, np.integer)



########
# Atom #
########

class Atom(object):
    
    def __init__(self, form, pos = None , B = 0):
        """initializes the class data and checks data types
          
            pos: 3 element list of relative atom position default: [0.0 ...]  
            struc: Structure Factor of the Atom default: N/A
            B: is the Debye-Waller factor (Temperature dependent movement)
        """
        if (pos):
            self.pos = pos
        else:
            self.pos = [0.0,0.0,0.0] 
        self.form = form
        self.B = B
        #check data types
        self._Atom__type_check__()    
        
    def __str__(self):
        """Returns the data in a structured form.
           
           returns:
           atom: atom x=pos[1] y=pos[2] z= pos[3] 
        """
        
        stri = ("atom: " + self.form.atom + " x=" + str(self.pos[0]) +
                " y=" + str(self.pos[1]) + " z=" + str(self.pos[2]) + 
                " B=" + str(self.B)
                )
        return stri
        
    def __type_check__(self):
        """checks the class data to have the right types
        
            pos    -> list
            pos[i] -> float
            struc  -> StructureFactor
            B -> number
        """
        assert isinstance(self.form, FormFactor) , 'form needs to be a FormFactor'
        assert (isinstance(self.pos, list) or (len(self.pos) != 3)) , 'a needs to be a list with 3 elements'
        for i in range(len(self.pos)):
            assert (isinstance(self.pos[i], _NUMBER_)) , 'pos' + str(i) + 'needs to be a number'
        assert (isinstance(self.B , _NUMBER_)) , 'pos' + str(i) + 'needs to be a number'
            
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
            self.a = [0.0,0.0,0.0,0.0,0.0,0.0]
        if(b): 
            self.b = b
        else:
            self.b = [0.0,0.0,0.0,0.0,0.0]
        #load data from path overwrites given data        
        if path:
            self.load_file(path)
        #check data types
        self._FormFactor__type_check__()
        
    def __str__(self):
        """Returns the data in a structured form.
        
            returns:
            element: 'atom'
            a1 = 'a[0]; a2 = 'a[1]; a3 = 'a[2]; a4 = 'a[3]; a5 = 'a[4] a6 = a[5];
            b1 = 'b[0]; b2 = 'b[1]; b3 = 'b[2]; b4 = 'b[3] b5 = b[4];
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
        assert isinstance(self.q, _NUMBER_) , 'q needs to be a Number'
        assert isinstance(self.value, _NUMBER_) , 'value needs to be a Number'
        assert isinstance(self.atom, str) , 'atom needs to be a string'
        assert (isinstance(self.a, list) or (len(self.a) != 6)) , 'a needs to be a list with 5 elements'
        for i in range(len(self.a)):
            assert (isinstance(self.a[i], _NUMBER_)) , 'a' + str(i) + 'needs to be a Number'
        assert (isinstance(self.b, list) or (len(self.b) != 5)) , 'b needs to be a list with 4 elements' 
        for i in range(len(self.b)):
            assert (isinstance(self.b[i], _NUMBER_)) , 'b' + str(i) + 'needs to be a Number'

    _FormFactor__type_check__ = __type_check__ #private copy of the function to avoid overload
    
    def get_value(self, q):
        """returns the structure factor value for the given q
        
           input q: absolute value of the wave vector      
        """
        assert isinstance(q,self.q.__class__), 'q needs to be ' + str(self.q.__class__)
        if (self.q != q): 
            q2 = (q*q) / (16*m.pi*m.pi) * 1e-20 
            self.value =(self.a[0] * m.exp(-self.b[0] * q2) + 
                         self.a[1] * m.exp(-self.b[1] * q2) + 
                         self.a[2] * m.exp(-self.b[2] * q2) + 
                         self.a[3] * m.exp(-self.b[3] * q2) +
                         self.a[4] * m.exp(-self.b[4] * q2) + 
                         self.a[5]
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
                elif words[0] == "c": self.a[5]= float(words[1])
                elif words[0] == "b1": self.b[0]= float(words[1])
                elif words[0] == "b2": self.b[1]= float(words[1])
                elif words[0] == "b3": self.b[2]= float(words[1])
                elif words[0] == "b4": self.b[3]= float(words[1])
                elif words[0] == "b5": self.b[4]= float(words[1])
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
            self._lattice_parameters_ = lattice_parameters
        else:
            self._lattice_parameters_ = [1.0,1.0,1.0]
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
        stri = ("Lattice parameters: a=" + str(self._lattice_parameters_[0]) +
                " b="+ str(self._lattice_parameters_[1]) +
                " c="+ str(self._lattice_parameters_[2]) + "\n"
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
        assert isinstance(self.damping, _NUMBER_) , 'damping needs to be a Number'
        assert (isinstance(self._lattice_parameters_, list) or (len(self._lattice_parameters_) != 3)) , 'lattice_parameters needs to be a list with 3 elements' 
        for i in range(len(self._lattice_parameters_)):
            assert (isinstance(self._lattice_parameters_[i], _NUMBER_)) , 'lattice_parameters[' + str(i) + '] needs to be a number'
        assert (isinstance(self.atoms, list)) , 'atoms needs to be a list ' 
        for i in range(len(self.atoms)):
            assert (isinstance(self.atoms[i], Atom)) , 'atoms[' + str(i) + '] needs to be a Atom'
    
    _CrystalStructure__type_check__ = __type_check__ #private copy of the function to avoid overload
            
    def _add_atom_(self, usedfile, atomtype):
        """adds the next atom in the file to the atom list """
        form = FormFactor(path=_SCRIPTPATH_+"/atoms/" + atomtype + ".at")
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
                elif words[0] == "a": self._lattice_parameters_[0] = float(words[1])
                elif words[0] == "b": self._lattice_parameters_[1] = float(words[1])
                elif words[0] == "c": self._lattice_parameters_[2] = float(words[1])
                elif words[0] == "damping": self.damping = float(words[1])
                
            line = f.readline()
        f.close()
    
    def save_file(self, path):
        """saves the structure to the file path"""
        f = open(path, 'w')
        f.write('a ' + str(self._lattice_parameters_[0]) + '\n')
        f.write('b ' + str(self._lattice_parameters_[1]) + '\n')
        f.write('c ' + str(self._lattice_parameters_[2]) + '\n')
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
                multi += (1 - atom.pos[i])*self._lattice_parameters_[i]*q[i] 
            structure_factor +=(atom.form.get_value(mq) *
                                 mc.exp(-1.0j * multi - 
                                        (self.damping / sinomega *
                                         self._lattice_parameters_[2] * 
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
    
    @property
    def a(self):
        """a lattice parameter"""
        return self._lattice_parameters_[0]

    @a.setter
    def a(self, value):
        self._lattice_parameters_[0] = value
        
    @property
    def b(self):
        """b lattice parameter"""
        return self._lattice_parameters_[1]

    @b.setter
    def b(self, value):
        self._lattice_parameters_[1] = value
        
    @property
    def c(self):
        """c lattice parameter"""
        return self._lattice_parameters_[2]

    @c.setter
    def c(self, value):
        self._lattice_parameters_[2] = value

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


##############
# BasicLayer #
##############

class BasicLayer(object):
        
    def __init__(self):
        """initializes the class data and checks data types"""
        pass
        
    def _calc_reflection_(self, q, sinomega):
        """calculates the complex reflection value of the Basic Layer."""
        return 1
    
    def _calc_phase_shift_(self, q, sinomega):
        """returns the phase shift of the Basic layer
        
            q: wave vector
           sinomega: sin(angle between incident beam and sample surface
        """
        return 1
        
    def _calc_thickness_(self):
        """returns the thickness of the Basic Layer"""
        return 1
    
    def _calc_N_(self):
        """returns the total number of layers of the Basic layer"""
        return 1
    
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
        return self._calc_thickness_()
    
    def get_N(self):
        """returns the total number of layers of the Basic layer"""
        return self._calc_N_()

###################
# BasicLayerCheck #
###################

class BasicLayerCheck(BasicLayer):
    
    def __init__(self):
        self._q_ = [0,0,0]
        self._sinomega_ = 0.0
        self._reflection_ = 0.0 + 0.0j
        self._phase_shift_ = 0.0 + 0.0j
    
    def get_reflection(self, q, sinomega):
        """returns the complex reflection value of the Basic Layer.
             
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

##################
# BasicLayerSave #
##################

class BasicLayerSave(BasicLayer):
    
    def __init__(self):
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
        """returns the complex reflection value of the Basic Layer.
             
           q: wave vector
           sinomega: sin(angle between incident beam and sample surface
        """
        return self._reflection_[self._check_q_s_(q, sinomega)]
    
    def get_phase_shift(self, q, sinomega):
        """returns the phase of the Basic Layer.

           q: scattering vector
           sinomega: sin angle between surface plane and incedent beam 
        """
        return self._phase_shift_[self._check_q_s_(q, sinomega)]

#########
# Layer #
#########

class Layer(BasicLayer):
        
    def __init__(self, crystal):
        """initializes the class data and checks data types
          
           crystal CrystalStructure of the film 
        """
        BasicLayer.__init__(self)
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
        a = self.crystal.a
        b = self.crystal.b
        c = self.crystal.c
        factor = self.crystal.damping * c / sinomega
        self._iaq_ = -1J * a * q[0] - factor 
        self._ibq_ = -1J * b * q[1] - factor 
        self._icq_ = -1J * c * q[2] - factor
    
    def _calc_amplitude_(self, q):
        """returns calculates the scattering amplitude"""
        a = self.crystal.a
        b = self.crystal.b
        R = 2.1879e-15 #R= 2.1879e-15 is the classical electron radius
        amp1d = 4 * m.pi * R / a / b #/ m.sqrt(sum([x**2 for x in q]))
        return (1 * amp1d * (1-mc.exp(self._icq_.real))**2)
    
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
    
    def _calc_thickness_(self):
        """returns the thickness of the film."""
        return (1 + self.delta)*self.crystal.c

##############
# LayerCheck #
##############

class LayerCheck(Layer, BasicLayerCheck):
    
    def __init__(self, crystal):
        Layer.__init__(self, crystal)
        BasicLayerCheck.__init__(self)

#############
# LayerSave #
#############

class LayerSave(Layer, BasicLayerSave):
    
    def __init__(self, crystal):
        Layer.__init__(self, crystal)
        BasicLayerSave.__init__(self)          

#############
# ThickFilm #
#############
 
class ThickFilm (Layer):

    def _calc_reflection_(self, q, sinomega):
        """calculates the complex reflection value of the Thick Film."""
        r = Layer._calc_reflection_(self, q, sinomega)
        r *= 1 / (1 - mc.exp(self._icq_))
        #r = 2 * r / (1 + m.sqrt(1 + 4 * abs(r**2)))
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
        super(self.__class__, self).__init__(crystal)
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
    
    def _calc_thickness_(self):
        """returns the thickness of the film."""
        return (self._layers_ + self.delta)*self.crystal.c
    
    def _calc_N_(self):
        """returns the total number of layers of the Basic layer"""
        return self._layers_

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

##########
# Sample #
##########

class Sample(BasicLayer):
        
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
            assert (isinstance(self._Layers_[i], BasicLayer)) , '_Layers_[' + str(i) + '] needs to be a Layer'
    
    _Sample__type_check__ = __type_check__ #private copy of the function to avoid overload
       
    def _calc_phase_shift_(self,  q, sinomega):
        """returns the phase of the Superlattice.

           q: scattering vector
           sinomega: sin of angle between surface and incedent beam
        """
        phase = 1.0+0.0j
        for layer in self._Layers_:
            phase *= layer.get_phase_shift(q, sinomega)
        return phase
    
    def _calc_reflection_(self, q, sinomega):
        """returns the complex reflection value of the Superlattice.
             
           q: wave vector
           sinomega: sin of angle between surface and incedent beam
        """
        r = 0.0 + 0.0j        
        phase = 1.0+0.0j
        for layer in reversed(self._Layers_):
            r += layer.get_reflection(q, sinomega)*phase
            phase *= layer.get_phase_shift(q, sinomega)
        return r 
            
    
    def _calc_thickness_(self):
        """returns the thickness of the sample."""
        thickness = 0.0 + 0.0j
        for layer in self._Layers_:
            thickness += layer.get_thickness()
        return thickness
    
    def _calc_N_(self):
        """returns the total number of layers of the Sample"""
        t = 0
        for film in self._Layers_:
            t += film.get_N() 
        return t
    
    def add_Layer(self, layer):
        """adds the layer to the top of the sample."""
        assert (isinstance(layer, (Layer, SuperLattice))) , 'layer needs to be a Layer'
        self._Layers_.append(layer)

#############
# MixSample #
#############

class MixSample(BasicLayer):
    
    def __init__(self, films = None, probabileties = None):
        """initializes the class data and checks data types
          
           films: different films in the Sample
        """
        self._films_ = []
        self._P_ = []
        if (films):
            for layer in films:
                self._films_.append(layer)
        if (probabileties):
            for p in probabileties:
                self._P_.append(p)
        self._MixSample__type_check__()
        
    def __type_check__(self):
        """checks the class data to have the right types
        
           _films_ -> list 
           _films_[i] -> Layer
           _P_ -> list 
           _P_[i] -> number
        """
        assert isinstance(self._films_, list) , '_films_ needs to be an list'
        for i in range(len(self._films_)):
            assert (isinstance(self._films_[i], (BasicLayer))) , '_films_[' + str(i) + '] needs to be a Basic Layer'
        assert isinstance(self._P_, list) , '_films_ needs to be an list'
        for i in range(len(self._P_)):
            assert (isinstance(self._P_[i], _NUMBER_)) , '_films_[' + str(i) + '] needs to be a Number'
        self.__check_P__()
    
    _MixSample__type_check__ = __type_check__ #private copy of the function to avoid overload
    
    def __check_P__(self):
        """checks that the probabilities are still right"""
        eps = 0.0001
        assert (len(self._P_) == len(self._films_)) , 'each film needs a probability'
        if (self._P_):
            if (sum(self._P_) > 1.0+eps or sum(self._P_) < 1.0-eps):
                assert (sum(self._P_) == 1.0) , 'sum of all probabilities needs to be 1 but is ' + str(sum(self._P_))    
    
    def _calc_reflection_(self, q, sinomega):
        """returns the complex reflection value of the Sample.
             
           q: wave vector
           sinomega: sin of angle between surface and incedent beam
        """
        self.__check_P__()
        r = 0.0 + 0.0j
        i = 0        
        for layer in self._films_:
            r += layer.get_reflection(q, sinomega) * self._P_[i]
            i += 1
        return r
    
    def _calc_phase_shift_(self, q, sinomega):
        """returns the phase shift of the layer"""
        p = 0
        i = 0
        for film in self._films_:
            p += film.get_phase_shift(q, sinomega) * self._P_[i]
            i += 1
        return p
    
    def _calc_thickness_(self):
        """returns the thickness of the film."""
        t = 0
        i = 0
        for film in self._films_:
            t += film.get_thickness() * self._P_[i]
            i += 1
        return t
    
    def _calc_N_(self):
        """returns the total number of layers of the Basic layer"""
        t = 0
        i = 0
        for film in self._films_:
            t += film.get_N() * self._P_[i]
            i += 1
        return t
    
    def get_c(self):
        return self._calc_thickness_() / self._calc_N_()

###################
# ComplexThinFilm #
###################

class ComplexThinFilm(MixSample):
        
    
    def __init__(self, crystal, layers=1.0, width = 0.01, filmType=ThinFilm):
        """initializes the class data and checks data types
          
           crystal: SimpleCrystalStructure of the film
           layers:  number of layers in the film
           width: width of the gausian distribution around layers
           filmType: the used filmtype to calculate the thinfilms
        """
        MixSample.__init__(self)
        self._crystal_ = crystal
        self._layers_ = layers
        self._width_ = width
        self._filmType_ = filmType  
        self._ComplexThinFilm__type_check__()
        self._init_P_()        
        
    def __type_check__(self):
        """checks the class data to have the right types
        
           layers -> number
           width -> number 
        """
        assert isinstance(self._layers_, _NUMBER_) , 'layers needs to be a Number'
        assert isinstance(self._width_, _NUMBER_) , 'width needs to be a Number'
        assert isinstance(self._crystal_, CrystalStructure) , 'crystal needs to be a CrystalStructure'
        #assert isinstance(self.filmType, (ThinFilm)) , 'filmType needs to be a ThinFilm'
        
    _ComplexThinFilm__type_check__ = __type_check__ #private copy of the function to avoid overload
    
    def _init_P_(self):
        """initializes the probabilities for different Thicknesses"""
        self._P_ = []
        self._films_ = []
        minimum = int(self._layers_ - 3 * self._width_) #int always rounds down
        maximum = int(round(self._layers_ + 3 * self._width_ + 1.5)) #+0.5 to round up +1 for the range
        summe = 0
        for i in range(minimum, maximum):
            value = m.exp((i - self._layers_)**2/(-2 * self._width_**2))
            self._P_.append(value)
            self._films_.append(self._filmType_(self._crystal_,i))
            summe += value
        for i in range(len(self._P_)):
            self._P_[i] /= summe
            
    def _init_P_perfect(self):
        """initializes the probabilities for different Thicknesses"""
        self._P_ = []
        self._films_ = []
        p2 = self._layers_ - int(self._layers_) #probabilety for second film
        p1 = 1-p2
        self._P_.append(p1)
        self._films_.append(self._filmType_(self._crystal_,int(self._layers_)))
        self._P_.append(p2)
        self._films_.append(self._filmType_(self._crystal_,int(self._layers_)+1))

            
########################
# ComplexThinFilmCheck #
########################

class ComplexThinFilmCheck(BasicLayerCheck, ComplexThinFilm):    
    
    def __init__(self, crystal, layers=1.0, width = 0.01, filmType=ThinFilm):
        ComplexThinFilm.__init__(self, crystal, layers, width, ThinFilm)
        BasicLayerCheck.__init__(self) 
 
#######################
# ComplexThinFilmSave #
#######################

class ComplexThinFilmSave(BasicLayerSave, ComplexThinFilm):    
    
    def __init__(self, crystal, layers=1.0, width = 0.01, filmType=ThinFilm):
        ComplexThinFilm.__init__(self, crystal, layers, width, ThinFilm)
        BasicLayerSave.__init__(self)   

##############
# MultiLayer #
##############

class Multilayer(Sample):
    
    def __init__(self, film1 = None, film2 = None, film3 = None, films = None):
        Sample.__init__(self, substrate=film1, electrode=film2, film=film3, layers=films)
        
    def get_c(self):
        return self.get_thickness()/self.get_N()    
       
################
# SuperLattice #
################

class SuperLattice(Multilayer):
        
    def __init__(self, film1 = None, film2 = None, film3 = None, films = None, bilayers=1):
        
        """initializes the class data and checks data types
          
           film1: first film of the superlattice 
           film2: second film of the superlattice 
           bilayers: number of times both films get repeted
        """
        Sample.__init__(self, substrate=film1, electrode=film2, film=film3, layers=films)
        self.bilayers = bilayers
        self._SuperLattice__type_check__()
        
    def __type_check__(self):
        """checks the class data to have the right types
        
           bilayers -> int 

        """
        assert isinstance(self.bilayers, (int)) , 'bilayers needs to be an int'
    
    _SuperLattice__type_check__ = __type_check__ #private copy of the function to avoid overload
       
    def _calc_phase_shift_(self,  q, sinomega):
        """returns the phase of the Superlattice.

           q: wave vector
           sinomega: sin(angle between incident beam and sample surface
        """
        return (Multilayer._calc_phase_shift_(self, q, sinomega) ** self.bilayers)
    
    def _calc_reflection_(self, q, sinomega):
        """returns the complex reflection value of the Superlattice.
             
           q: wave vector
           sinomega: sin(angle between incident beam and sample surface
        """
        r = Multilayer._calc_reflection_(self, q, sinomega)
        p = Multilayer._calc_phase_shift_(self, q, sinomega)
        reflection = 0.0 + 0.0j
        phase = 1.0
        for i in range(self.bilayers):
            reflection += (r * phase)
            phase *= p
        return reflection
    
    def _calc_thickness_(self):
        """returns the total thickness of the superlatice"""
        return self.get_multilayer_thickness() *self.bilayers
    
    def _calc_N_(self):
        """returns the total number of layers in the Superlattice"""
        return self.get_n() * self.bilayers
    
    def get_c(self):
        """returns the average c lattice parameter"""
        return self.get_multilayer_thickness() / self.get_n()
    
    def get_multilayer_thickness(self):
        """returns the Bilayer thickness of the superlatice"""
        return Multilayer._calc_thickness_(self)
    
    def get_n(self):
        """returns the number of layers per multilayer"""
        return Multilayer._calc_N_(self)
    
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
            l += crystal.c * self.layers[n]
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
    q[0] = h*2*m.pi/base_struc.a
    q[1] = k*2*m.pi/base_struc.b
    scan = []
    for i in range(int(steps_l)):
        l=l_min+i*l_step
        q[2] = l*2*m.pi/base_struc.c
        absq = m.sqrt(sum([x**2 for x in q]))
        if (sinomegain == 0):
            sinomega = absq * wavelength / 4 / m.pi
        else:
            sinomega = sinomegain
        lin = abs(sample.get_reflection(q, sinomega))
        lin *= lin
        lin *= direct
        lin += background
        #lin *= 1/m.sin(2*m.asin(sinomega))*(1+m.cos(2*m.asin(sinomega))**2)/2/sinomega
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
    a = max(struc1.a,struc2.a)
    b = max(struc1.b,struc2.b)
    c = max(struc1.c,struc2.c)
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
    crystal.a = a
    crystal.b = b
    crystal.c = c
    
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
            form = FormFactor(path=_SCRIPTPATH_+"/atoms/" + atomtype + ".at")
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
            form = FormFactor(path=_SCRIPTPATH_+"/atoms/" + atomtype + ".at")
            pos = [0.0,0.0,0.0]
            pos[0] =(j + 0.5 + random.randrange(-1,1) * size) / L
            pos[1] =(    0.5 + random.randrange(-1,1) * size)
            pos[2] =(i + 0.5 + zoffset + random.randrange(-1,1) * size) / layers
            atom = Atom(form, pos)
            crystal.add_atom(atom)
            #oxygen 1
            form = FormFactor(path=_SCRIPTPATH_+"/atoms/O.at")
            pos = [0.0,0.0,0.0]
            pos[0] =(j + 0.5 + random.randrange(-1,1) * size) / L
            pos[1] =(    0.5 + random.randrange(-1,1) * size)
            pos[2] =(i + 1.0 + zoffset + random.randrange(-1,1) * size) / layers
            atom = Atom(form, pos)
            crystal.add_atom(atom)
            #oxygen 2
            form = FormFactor(path=_SCRIPTPATH_+"/atoms/O.at")
            pos = [0.0,0.0,0.0]
            pos[0] =(j + 0.5 + random.randrange(-1,1) * size) / L
            pos[1] =(    0.0 + random.randrange(-1,1) * size)
            pos[2] =(i + 0.5 + zoffset + random.randrange(-1,1) * size) / layers
            atom = Atom(form, pos)
            crystal.add_atom(atom)
            #oxygen 3
            form = FormFactor(path=_SCRIPTPATH_+"/atoms/O.at")
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