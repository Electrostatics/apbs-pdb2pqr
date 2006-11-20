# This file was created automatically by SWIG 1.3.29.
# Don't modify this file, modify the SWIG interface instead.
# This file is compatible with both classic and new-style classes.

import _apbslib
import new
new_instancemethod = new.instancemethod
def _swig_setattr_nondynamic(self,class_type,name,value,static=1):
    if (name == "thisown"): return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'PySwigObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    if (not static) or hasattr(self,name):
        self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)

def _swig_setattr(self,class_type,name,value):
    return _swig_setattr_nondynamic(self,class_type,name,value,0)

def _swig_getattr(self,class_type,name):
    if (name == "thisown"): return self.this.own()
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError,name

def _swig_repr(self):
    try: strthis = "proxy of " + self.this.__repr__()
    except: strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

import types
try:
    _object = types.ObjectType
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0
del types


class Valist(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Valist, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Valist, name)
    __repr__ = _swig_repr
    def __init__(self, *args): 
        this = _apbslib.new_Valist(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _apbslib.delete_Valist
    __del__ = lambda self : None;
    __swig_setmethods__["number"] = _apbslib.Valist_number_set
    __swig_getmethods__["number"] = _apbslib.Valist_number_get
    if _newclass:number = property(_apbslib.Valist_number_get, _apbslib.Valist_number_set)
Valist_swigregister = _apbslib.Valist_swigregister
Valist_swigregister(Valist)

Valist_getAtomList = _apbslib.Valist_getAtomList
Valist_getAtom = _apbslib.Valist_getAtom
class Vatom(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Vatom, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Vatom, name)
    __repr__ = _swig_repr
    def __init__(self, *args): 
        this = _apbslib.new_Vatom(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _apbslib.delete_Vatom
    __del__ = lambda self : None;
    __swig_setmethods__["id"] = _apbslib.Vatom_id_set
    __swig_getmethods__["id"] = _apbslib.Vatom_id_get
    if _newclass:id = property(_apbslib.Vatom_id_get, _apbslib.Vatom_id_set)
Vatom_swigregister = _apbslib.Vatom_swigregister
Vatom_swigregister(Vatom)

Vatom_getPosition = _apbslib.Vatom_getPosition
Vatom_setCharge = _apbslib.Vatom_setCharge
Vatom_getCharge = _apbslib.Vatom_getCharge
Vatom_getRadius = _apbslib.Vatom_getRadius
class MGparm(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, MGparm, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, MGparm, name)
    __repr__ = _swig_repr
    def __init__(self, *args): 
        this = _apbslib.new_MGparm(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _apbslib.delete_MGparm
    __del__ = lambda self : None;
    __swig_setmethods__["type"] = _apbslib.MGparm_type_set
    __swig_getmethods__["type"] = _apbslib.MGparm_type_get
    if _newclass:type = property(_apbslib.MGparm_type_get, _apbslib.MGparm_type_set)
MGparm_swigregister = _apbslib.MGparm_swigregister
MGparm_swigregister(MGparm)

MGparm_setCenterX = _apbslib.MGparm_setCenterX
MGparm_setCenterY = _apbslib.MGparm_setCenterY
MGparm_setCenterZ = _apbslib.MGparm_setCenterZ
class PBEparm(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, PBEparm, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, PBEparm, name)
    __repr__ = _swig_repr
    def __init__(self, *args): 
        this = _apbslib.new_PBEparm(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _apbslib.delete_PBEparm
    __del__ = lambda self : None;
    __swig_setmethods__["temp"] = _apbslib.PBEparm_temp_set
    __swig_getmethods__["temp"] = _apbslib.PBEparm_temp_get
    if _newclass:temp = property(_apbslib.PBEparm_temp_get, _apbslib.PBEparm_temp_set)
    __swig_setmethods__["pdie"] = _apbslib.PBEparm_pdie_set
    __swig_getmethods__["pdie"] = _apbslib.PBEparm_pdie_get
    if _newclass:pdie = property(_apbslib.PBEparm_pdie_get, _apbslib.PBEparm_pdie_set)
    __swig_setmethods__["sdie"] = _apbslib.PBEparm_sdie_set
    __swig_getmethods__["sdie"] = _apbslib.PBEparm_sdie_get
    if _newclass:sdie = property(_apbslib.PBEparm_sdie_get, _apbslib.PBEparm_sdie_set)
    __swig_setmethods__["molid"] = _apbslib.PBEparm_molid_set
    __swig_getmethods__["molid"] = _apbslib.PBEparm_molid_get
    if _newclass:molid = property(_apbslib.PBEparm_molid_get, _apbslib.PBEparm_molid_set)
PBEparm_swigregister = _apbslib.PBEparm_swigregister
PBEparm_swigregister(PBEparm)

class Vcom(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Vcom, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Vcom, name)
    __repr__ = _swig_repr
    def __init__(self, *args): 
        this = _apbslib.new_Vcom(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _apbslib.delete_Vcom
    __del__ = lambda self : None;
Vcom_swigregister = _apbslib.Vcom_swigregister
Vcom_swigregister(Vcom)

Vcom_ctor = _apbslib.Vcom_ctor
Vcom_size = _apbslib.Vcom_size
Vcom_rank = _apbslib.Vcom_rank
class Vmem(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Vmem, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Vmem, name)
    __repr__ = _swig_repr
    def __init__(self, *args): 
        this = _apbslib.new_Vmem(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _apbslib.delete_Vmem
    __del__ = lambda self : None;
Vmem_swigregister = _apbslib.Vmem_swigregister
Vmem_swigregister(Vmem)

Vmem_ctor = _apbslib.Vmem_ctor
class Vpmg(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Vpmg, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Vpmg, name)
    __repr__ = _swig_repr
    def __init__(self, *args): 
        this = _apbslib.new_Vpmg(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _apbslib.delete_Vpmg
    __del__ = lambda self : None;
Vpmg_swigregister = _apbslib.Vpmg_swigregister
Vpmg_swigregister(Vpmg)

class Vpbe(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Vpbe, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Vpbe, name)
    __repr__ = _swig_repr
    def __init__(self, *args): 
        this = _apbslib.new_Vpbe(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _apbslib.delete_Vpbe
    __del__ = lambda self : None;
    __swig_setmethods__["acc"] = _apbslib.Vpbe_acc_set
    __swig_getmethods__["acc"] = _apbslib.Vpbe_acc_get
    if _newclass:acc = property(_apbslib.Vpbe_acc_get, _apbslib.Vpbe_acc_set)
Vpbe_swigregister = _apbslib.Vpbe_swigregister
Vpbe_swigregister(Vpbe)

class NOsh_calc(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, NOsh_calc, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, NOsh_calc, name)
    __repr__ = _swig_repr
    def __init__(self, *args): 
        this = _apbslib.new_NOsh_calc(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _apbslib.delete_NOsh_calc
    __del__ = lambda self : None;
    __swig_setmethods__["mgparm"] = _apbslib.NOsh_calc_mgparm_set
    __swig_getmethods__["mgparm"] = _apbslib.NOsh_calc_mgparm_get
    if _newclass:mgparm = property(_apbslib.NOsh_calc_mgparm_get, _apbslib.NOsh_calc_mgparm_set)
    __swig_setmethods__["femparm"] = _apbslib.NOsh_calc_femparm_set
    __swig_getmethods__["femparm"] = _apbslib.NOsh_calc_femparm_get
    if _newclass:femparm = property(_apbslib.NOsh_calc_femparm_get, _apbslib.NOsh_calc_femparm_set)
    __swig_setmethods__["pbeparm"] = _apbslib.NOsh_calc_pbeparm_set
    __swig_getmethods__["pbeparm"] = _apbslib.NOsh_calc_pbeparm_get
    if _newclass:pbeparm = property(_apbslib.NOsh_calc_pbeparm_get, _apbslib.NOsh_calc_pbeparm_set)
    __swig_setmethods__["calctype"] = _apbslib.NOsh_calc_calctype_set
    __swig_getmethods__["calctype"] = _apbslib.NOsh_calc_calctype_get
    if _newclass:calctype = property(_apbslib.NOsh_calc_calctype_get, _apbslib.NOsh_calc_calctype_set)
NOsh_calc_swigregister = _apbslib.NOsh_calc_swigregister
NOsh_calc_swigregister(NOsh_calc)

class NOsh(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, NOsh, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, NOsh, name)
    __repr__ = _swig_repr
    def __init__(self, *args): 
        this = _apbslib.new_NOsh(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _apbslib.delete_NOsh
    __del__ = lambda self : None;
    __swig_setmethods__["ncalc"] = _apbslib.NOsh_ncalc_set
    __swig_getmethods__["ncalc"] = _apbslib.NOsh_ncalc_get
    if _newclass:ncalc = property(_apbslib.NOsh_ncalc_get, _apbslib.NOsh_ncalc_set)
    __swig_setmethods__["nprint"] = _apbslib.NOsh_nprint_set
    __swig_getmethods__["nprint"] = _apbslib.NOsh_nprint_get
    if _newclass:nprint = property(_apbslib.NOsh_nprint_get, _apbslib.NOsh_nprint_set)
    __swig_setmethods__["nelec"] = _apbslib.NOsh_nelec_set
    __swig_getmethods__["nelec"] = _apbslib.NOsh_nelec_get
    if _newclass:nelec = property(_apbslib.NOsh_nelec_get, _apbslib.NOsh_nelec_set)
    __swig_setmethods__["nmol"] = _apbslib.NOsh_nmol_set
    __swig_getmethods__["nmol"] = _apbslib.NOsh_nmol_get
    if _newclass:nmol = property(_apbslib.NOsh_nmol_get, _apbslib.NOsh_nmol_set)
    __swig_setmethods__["printwhat"] = _apbslib.NOsh_printwhat_set
    __swig_getmethods__["printwhat"] = _apbslib.NOsh_printwhat_get
    if _newclass:printwhat = property(_apbslib.NOsh_printwhat_get, _apbslib.NOsh_printwhat_set)
NOsh_swigregister = _apbslib.NOsh_swigregister
NOsh_swigregister(NOsh)

NPT_ENERGY = _apbslib.NPT_ENERGY
NPT_FORCE = _apbslib.NPT_FORCE
NOsh_getCalc = _apbslib.NOsh_getCalc
NOsh_elecname = _apbslib.NOsh_elecname
NOsh_elec2calc = _apbslib.NOsh_elec2calc
NOsh_printWhat = _apbslib.NOsh_printWhat
NOsh_parseInputFile = _apbslib.NOsh_parseInputFile
NOsh_ctor = _apbslib.NOsh_ctor
class AtomForce(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, AtomForce, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, AtomForce, name)
    __repr__ = _swig_repr
    def __init__(self, *args): 
        this = _apbslib.new_AtomForce(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _apbslib.delete_AtomForce
    __del__ = lambda self : None;
AtomForce_swigregister = _apbslib.AtomForce_swigregister
AtomForce_swigregister(AtomForce)

new_valist = _apbslib.new_valist
get_Valist = _apbslib.get_Valist
new_gridlist = _apbslib.new_gridlist
new_pmglist = _apbslib.new_pmglist
get_Vpmg = _apbslib.get_Vpmg
new_pmgplist = _apbslib.new_pmgplist
new_pbelist = _apbslib.new_pbelist
get_Vpbe = _apbslib.get_Vpbe
new_atomforcelist = _apbslib.new_atomforcelist
delete_atomforcelist = _apbslib.delete_atomforcelist
delete_valist = _apbslib.delete_valist
delete_gridlist = _apbslib.delete_gridlist
delete_pmglist = _apbslib.delete_pmglist
delete_pmgplist = _apbslib.delete_pmgplist
delete_pbelist = _apbslib.delete_pbelist
delete_Nosh = _apbslib.delete_Nosh
delete_Com = _apbslib.delete_Com
delete_Mem = _apbslib.delete_Mem
get_AtomForce = _apbslib.get_AtomForce
make_Valist = _apbslib.make_Valist
double_array = _apbslib.double_array
int_array = _apbslib.int_array
delete_double_array = _apbslib.delete_double_array
delete_int_array = _apbslib.delete_int_array
get_entry = _apbslib.get_entry
set_entry = _apbslib.set_entry
parseInputFromString = _apbslib.parseInputFromString
Valist_load = _apbslib.Valist_load
NOsh_setupElecCalc = _apbslib.NOsh_setupElecCalc
wrap_forceMG = _apbslib.wrap_forceMG
getAtomPosition = _apbslib.getAtomPosition
getPotentials = _apbslib.getPotentials
getEnergies = _apbslib.getEnergies
getForces = _apbslib.getForces
loadMolecules = _apbslib.loadMolecules
killMolecules = _apbslib.killMolecules
loadDielMaps = _apbslib.loadDielMaps
killDielMaps = _apbslib.killDielMaps
loadKappaMaps = _apbslib.loadKappaMaps
killKappaMaps = _apbslib.killKappaMaps
loadChargeMaps = _apbslib.loadChargeMaps
killChargeMaps = _apbslib.killChargeMaps
printPBEPARM = _apbslib.printPBEPARM
printMGPARM = _apbslib.printMGPARM
initMG = _apbslib.initMG
killMG = _apbslib.killMG
solveMG = _apbslib.solveMG
setPartMG = _apbslib.setPartMG
killEnergy = _apbslib.killEnergy
killForce = _apbslib.killForce
writedataMG = _apbslib.writedataMG
writematMG = _apbslib.writematMG
printForce = _apbslib.printForce
startVio = _apbslib.startVio
Vacc_molAcc = _apbslib.Vacc_molAcc
Vacc_vdwAcc = _apbslib.Vacc_vdwAcc
energyMG = _apbslib.energyMG
printEnergy = _apbslib.printEnergy
returnEnergy = _apbslib.returnEnergy


