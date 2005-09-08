# This file was created automatically by SWIG.
# Don't modify this file, modify the SWIG interface instead.
# This file is compatible with both classic and new-style classes.

import _apbslib

def _swig_setattr_nondynamic(self,class_type,name,value,static=1):
    if (name == "this"):
        if isinstance(value, class_type):
            self.__dict__[name] = value.this
            if hasattr(value,"thisown"): self.__dict__["thisown"] = value.thisown
            del value.thisown
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    if (not static) or hasattr(self,name) or (name == "thisown"):
        self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)

def _swig_setattr(self,class_type,name,value):
    return _swig_setattr_nondynamic(self,class_type,name,value,0)

def _swig_getattr(self,class_type,name):
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError,name

import types
try:
    _object = types.ObjectType
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0
del types


APBS_SWIG = _apbslib.APBS_SWIG
class MGparm(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, MGparm, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, MGparm, name)
    def __repr__(self):
        return "<%s.%s; proxy of C MGparm instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        _swig_setattr(self, MGparm, 'this', _apbslib.new_MGparm(*args))
        _swig_setattr(self, MGparm, 'thisown', 1)
    def __del__(self, destroy=_apbslib.delete_MGparm):
        try:
            if self.thisown: destroy(self)
        except: pass

    __swig_setmethods__["type"] = _apbslib.MGparm_type_set
    __swig_getmethods__["type"] = _apbslib.MGparm_type_get
    if _newclass:type = property(_apbslib.MGparm_type_get, _apbslib.MGparm_type_set)

class MGparmPtr(MGparm):
    def __init__(self, this):
        _swig_setattr(self, MGparm, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, MGparm, 'thisown', 0)
        _swig_setattr(self, MGparm,self.__class__,MGparm)
_apbslib.MGparm_swigregister(MGparmPtr)

class PBEparm(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, PBEparm, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, PBEparm, name)
    def __repr__(self):
        return "<%s.%s; proxy of C PBEparm instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        _swig_setattr(self, PBEparm, 'this', _apbslib.new_PBEparm(*args))
        _swig_setattr(self, PBEparm, 'thisown', 1)
    def __del__(self, destroy=_apbslib.delete_PBEparm):
        try:
            if self.thisown: destroy(self)
        except: pass

    __swig_setmethods__["temp"] = _apbslib.PBEparm_temp_set
    __swig_getmethods__["temp"] = _apbslib.PBEparm_temp_get
    if _newclass:temp = property(_apbslib.PBEparm_temp_get, _apbslib.PBEparm_temp_set)

class PBEparmPtr(PBEparm):
    def __init__(self, this):
        _swig_setattr(self, PBEparm, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, PBEparm, 'thisown', 0)
        _swig_setattr(self, PBEparm,self.__class__,PBEparm)
_apbslib.PBEparm_swigregister(PBEparmPtr)

class Vcom(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Vcom, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Vcom, name)
    def __repr__(self):
        return "<%s.%s; proxy of C Vcom instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        _swig_setattr(self, Vcom, 'this', _apbslib.new_Vcom(*args))
        _swig_setattr(self, Vcom, 'thisown', 1)
    def __del__(self, destroy=_apbslib.delete_Vcom):
        try:
            if self.thisown: destroy(self)
        except: pass


class VcomPtr(Vcom):
    def __init__(self, this):
        _swig_setattr(self, Vcom, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, Vcom, 'thisown', 0)
        _swig_setattr(self, Vcom,self.__class__,Vcom)
_apbslib.Vcom_swigregister(VcomPtr)


Vcom_ctor = _apbslib.Vcom_ctor

Vcom_size = _apbslib.Vcom_size

Vcom_rank = _apbslib.Vcom_rank
class Vmem(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Vmem, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Vmem, name)
    def __repr__(self):
        return "<%s.%s; proxy of C Vmem instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        _swig_setattr(self, Vmem, 'this', _apbslib.new_Vmem(*args))
        _swig_setattr(self, Vmem, 'thisown', 1)
    def __del__(self, destroy=_apbslib.delete_Vmem):
        try:
            if self.thisown: destroy(self)
        except: pass


class VmemPtr(Vmem):
    def __init__(self, this):
        _swig_setattr(self, Vmem, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, Vmem, 'thisown', 0)
        _swig_setattr(self, Vmem,self.__class__,Vmem)
_apbslib.Vmem_swigregister(VmemPtr)


Vmem_ctor = _apbslib.Vmem_ctor
class Vpmg(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Vpmg, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Vpmg, name)
    def __repr__(self):
        return "<%s.%s; proxy of C Vpmg instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        _swig_setattr(self, Vpmg, 'this', _apbslib.new_Vpmg(*args))
        _swig_setattr(self, Vpmg, 'thisown', 1)
    def __del__(self, destroy=_apbslib.delete_Vpmg):
        try:
            if self.thisown: destroy(self)
        except: pass


class VpmgPtr(Vpmg):
    def __init__(self, this):
        _swig_setattr(self, Vpmg, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, Vpmg, 'thisown', 0)
        _swig_setattr(self, Vpmg,self.__class__,Vpmg)
_apbslib.Vpmg_swigregister(VpmgPtr)

class Vpbe(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Vpbe, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Vpbe, name)
    def __repr__(self):
        return "<%s.%s; proxy of C Vpbe instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        _swig_setattr(self, Vpbe, 'this', _apbslib.new_Vpbe(*args))
        _swig_setattr(self, Vpbe, 'thisown', 1)
    def __del__(self, destroy=_apbslib.delete_Vpbe):
        try:
            if self.thisown: destroy(self)
        except: pass

    __swig_setmethods__["acc"] = _apbslib.Vpbe_acc_set
    __swig_getmethods__["acc"] = _apbslib.Vpbe_acc_get
    if _newclass:acc = property(_apbslib.Vpbe_acc_get, _apbslib.Vpbe_acc_set)

class VpbePtr(Vpbe):
    def __init__(self, this):
        _swig_setattr(self, Vpbe, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, Vpbe, 'thisown', 0)
        _swig_setattr(self, Vpbe,self.__class__,Vpbe)
_apbslib.Vpbe_swigregister(VpbePtr)

class NOsh_calc(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, NOsh_calc, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, NOsh_calc, name)
    def __repr__(self):
        return "<%s.%s; proxy of C NOsh_calc instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
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
    def __init__(self, *args):
        _swig_setattr(self, NOsh_calc, 'this', _apbslib.new_NOsh_calc(*args))
        _swig_setattr(self, NOsh_calc, 'thisown', 1)
    def __del__(self, destroy=_apbslib.delete_NOsh_calc):
        try:
            if self.thisown: destroy(self)
        except: pass


class NOsh_calcPtr(NOsh_calc):
    def __init__(self, this):
        _swig_setattr(self, NOsh_calc, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, NOsh_calc, 'thisown', 0)
        _swig_setattr(self, NOsh_calc,self.__class__,NOsh_calc)
_apbslib.NOsh_calc_swigregister(NOsh_calcPtr)

class NOsh(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, NOsh, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, NOsh, name)
    def __repr__(self):
        return "<%s.%s; proxy of C NOsh instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        _swig_setattr(self, NOsh, 'this', _apbslib.new_NOsh(*args))
        _swig_setattr(self, NOsh, 'thisown', 1)
    __swig_setmethods__["ncalc"] = _apbslib.NOsh_ncalc_set
    __swig_getmethods__["ncalc"] = _apbslib.NOsh_ncalc_get
    if _newclass:ncalc = property(_apbslib.NOsh_ncalc_get, _apbslib.NOsh_ncalc_set)
    __swig_setmethods__["nprint"] = _apbslib.NOsh_nprint_set
    __swig_getmethods__["nprint"] = _apbslib.NOsh_nprint_get
    if _newclass:nprint = property(_apbslib.NOsh_nprint_get, _apbslib.NOsh_nprint_set)
    __swig_setmethods__["nelec"] = _apbslib.NOsh_nelec_set
    __swig_getmethods__["nelec"] = _apbslib.NOsh_nelec_get
    if _newclass:nelec = property(_apbslib.NOsh_nelec_get, _apbslib.NOsh_nelec_set)
    __swig_setmethods__["printwhat"] = _apbslib.NOsh_printwhat_set
    __swig_getmethods__["printwhat"] = _apbslib.NOsh_printwhat_get
    if _newclass:printwhat = property(_apbslib.NOsh_printwhat_get, _apbslib.NOsh_printwhat_set)
    def __del__(self, destroy=_apbslib.delete_NOsh):
        try:
            if self.thisown: destroy(self)
        except: pass


class NOshPtr(NOsh):
    def __init__(self, this):
        _swig_setattr(self, NOsh, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, NOsh, 'thisown', 0)
        _swig_setattr(self, NOsh,self.__class__,NOsh)
_apbslib.NOsh_swigregister(NOshPtr)

NPT_ENERGY = _apbslib.NPT_ENERGY
NPT_FORCE = _apbslib.NPT_FORCE

NOsh_getCalc = _apbslib.NOsh_getCalc

NOsh_elecname = _apbslib.NOsh_elecname

NOsh_elec2calc = _apbslib.NOsh_elec2calc

NOsh_printWhat = _apbslib.NOsh_printWhat

NOsh_ctor2 = _apbslib.NOsh_ctor2

NOsh_parseFile = _apbslib.NOsh_parseFile
class AtomForce(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, AtomForce, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, AtomForce, name)
    def __repr__(self):
        return "<%s.%s; proxy of C AtomForce instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        _swig_setattr(self, AtomForce, 'this', _apbslib.new_AtomForce(*args))
        _swig_setattr(self, AtomForce, 'thisown', 1)
    def __del__(self, destroy=_apbslib.delete_AtomForce):
        try:
            if self.thisown: destroy(self)
        except: pass


class AtomForcePtr(AtomForce):
    def __init__(self, this):
        _swig_setattr(self, AtomForce, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, AtomForce, 'thisown', 0)
        _swig_setattr(self, AtomForce,self.__class__,AtomForce)
_apbslib.AtomForce_swigregister(AtomForcePtr)


new_valist = _apbslib.new_valist

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

wrap_forceMG = _apbslib.wrap_forceMG

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

npenergyMG = _apbslib.npenergyMG

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

