# This file was created automatically by SWIG.
import apbslibc
class MGparmPtr :
    def __init__(self,this):
        self.this = this
        self.thisown = 0
    def __del__(self):
        if self.thisown == 1 :
            apbslibc.delete_MGparm(self.this)
    def __setattr__(self,name,value):
        if name == "type" :
            apbslibc.MGparm_type_set(self.this,value)
            return
        self.__dict__[name] = value
    def __getattr__(self,name):
        if name == "type" : 
            return apbslibc.MGparm_type_get(self.this)
        raise AttributeError,name
    def __repr__(self):
        return "<C MGparm instance>"
class MGparm(MGparmPtr):
    def __init__(self) :
        self.this = apbslibc.new_MGparm()
        self.thisown = 1




class PBEparmPtr :
    def __init__(self,this):
        self.this = this
        self.thisown = 0
    def __del__(self):
        if self.thisown == 1 :
            apbslibc.delete_PBEparm(self.this)
    def __setattr__(self,name,value):
        if name == "temp" :
            apbslibc.PBEparm_temp_set(self.this,value)
            return
        self.__dict__[name] = value
    def __getattr__(self,name):
        if name == "temp" : 
            return apbslibc.PBEparm_temp_get(self.this)
        raise AttributeError,name
    def __repr__(self):
        return "<C PBEparm instance>"
class PBEparm(PBEparmPtr):
    def __init__(self) :
        self.this = apbslibc.new_PBEparm()
        self.thisown = 1




class VcomPtr :
    def __init__(self,this):
        self.this = this
        self.thisown = 0
    def __del__(self):
        if self.thisown == 1 :
            apbslibc.delete_Vcom(self.this)
    def __repr__(self):
        return "<C Vcom instance>"
class Vcom(VcomPtr):
    def __init__(self) :
        self.this = apbslibc.new_Vcom()
        self.thisown = 1




class VmemPtr :
    def __init__(self,this):
        self.this = this
        self.thisown = 0
    def __del__(self):
        if self.thisown == 1 :
            apbslibc.delete_Vmem(self.this)
    def __repr__(self):
        return "<C Vmem instance>"
class Vmem(VmemPtr):
    def __init__(self) :
        self.this = apbslibc.new_Vmem()
        self.thisown = 1




class VpmgPtr :
    def __init__(self,this):
        self.this = this
        self.thisown = 0
    def __del__(self):
        if self.thisown == 1 :
            apbslibc.delete_Vpmg(self.this)
    def __repr__(self):
        return "<C Vpmg instance>"
class Vpmg(VpmgPtr):
    def __init__(self) :
        self.this = apbslibc.new_Vpmg()
        self.thisown = 1




class NOsh_calcPtr :
    def __init__(self,this):
        self.this = this
        self.thisown = 0
    def __setattr__(self,name,value):
        if name == "mgparm" :
            apbslibc.NOsh_calc_mgparm_set(self.this,value.this)
            return
        if name == "femparm" :
            apbslibc.NOsh_calc_femparm_set(self.this,value)
            return
        if name == "pbeparm" :
            apbslibc.NOsh_calc_pbeparm_set(self.this,value.this)
            return
        if name == "calctype" :
            apbslibc.NOsh_calc_calctype_set(self.this,value)
            return
        self.__dict__[name] = value
    def __getattr__(self,name):
        if name == "mgparm" : 
            return MGparmPtr(apbslibc.NOsh_calc_mgparm_get(self.this))
        if name == "femparm" : 
            return apbslibc.NOsh_calc_femparm_get(self.this)
        if name == "pbeparm" : 
            return PBEparmPtr(apbslibc.NOsh_calc_pbeparm_get(self.this))
        if name == "calctype" : 
            return apbslibc.NOsh_calc_calctype_get(self.this)
        raise AttributeError,name
    def __repr__(self):
        return "<C NOsh_calc instance>"
class NOsh_calc(NOsh_calcPtr):
    def __init__(self,this):
        self.this = this




class NOshPtr :
    def __init__(self,this):
        self.this = this
        self.thisown = 0
    def __del__(self):
        if self.thisown == 1 :
            apbslibc.delete_NOsh(self.this)
    def __setattr__(self,name,value):
        if name == "ncalc" :
            apbslibc.NOsh_ncalc_set(self.this,value)
            return
        if name == "nprint" :
            apbslibc.NOsh_nprint_set(self.this,value)
            return
        if name == "nelec" :
            apbslibc.NOsh_nelec_set(self.this,value)
            return
        if name == "printwhat" :
            apbslibc.NOsh_printwhat_set(self.this,value)
            return
        self.__dict__[name] = value
    def __getattr__(self,name):
        if name == "ncalc" : 
            return apbslibc.NOsh_ncalc_get(self.this)
        if name == "nprint" : 
            return apbslibc.NOsh_nprint_get(self.this)
        if name == "nelec" : 
            return apbslibc.NOsh_nelec_get(self.this)
        if name == "printwhat" : 
            return apbslibc.NOsh_printwhat_get(self.this)
        raise AttributeError,name
    def __repr__(self):
        return "<C NOsh instance>"
class NOsh(NOshPtr):
    def __init__(self) :
        self.this = apbslibc.new_NOsh()
        self.thisown = 1




class AtomForcePtr :
    def __init__(self,this):
        self.this = this
        self.thisown = 0
    def __del__(self):
        if self.thisown == 1 :
            apbslibc.delete_AtomForce(self.this)
    def __repr__(self):
        return "<C AtomForce instance>"
class AtomForce(AtomForcePtr):
    def __init__(self) :
        self.this = apbslibc.new_AtomForce()
        self.thisown = 1






#-------------- FUNCTION WRAPPERS ------------------

def Vcom_ctor(arg0):
    val = apbslibc.Vcom_ctor(arg0)
    val = VcomPtr(val)
    return val

Vcom_dtor = apbslibc.Vcom_dtor

def Vcom_size(arg0):
    val = apbslibc.Vcom_size(arg0.this)
    return val

def Vcom_rank(arg0):
    val = apbslibc.Vcom_rank(arg0.this)
    return val

def Vmem_ctor(arg0):
    val = apbslibc.Vmem_ctor(arg0)
    val = VmemPtr(val)
    return val

Vmem_dtor = apbslibc.Vmem_dtor

def NOsh_getCalc(arg0,arg1):
    val = apbslibc.NOsh_getCalc(arg0.this,arg1)
    val = NOsh_calcPtr(val)
    return val

def NOsh_elecname(arg0,arg1):
    val = apbslibc.NOsh_elecname(arg0.this,arg1)
    return val

def NOsh_elec2calc(arg0,arg1):
    val = apbslibc.NOsh_elec2calc(arg0.this,arg1)
    return val

def NOsh_printWhat(arg0,arg1):
    val = apbslibc.NOsh_printWhat(arg0.this,arg1)
    return val

def NOsh_ctor2(arg0,arg1,arg2):
    val = apbslibc.NOsh_ctor2(arg0.this,arg1,arg2)
    return val

NOsh_dtor = apbslibc.NOsh_dtor

def NOsh_parseFile(arg0,arg1):
    val = apbslibc.NOsh_parseFile(arg0.this,arg1)
    return val

ptrcast = apbslibc.ptrcast

ptrvalue = apbslibc.ptrvalue

ptrset = apbslibc.ptrset

ptrcreate = apbslibc.ptrcreate

ptrfree = apbslibc.ptrfree

ptradd = apbslibc.ptradd

ptrmap = apbslibc.ptrmap

new_valist = apbslibc.new_valist

new_gridlist = apbslibc.new_gridlist

new_pmglist = apbslibc.new_pmglist

def get_Vpmg(arg0,arg1):
    val = apbslibc.get_Vpmg(arg0,arg1)
    val = VpmgPtr(val)
    return val

new_pmgplist = apbslibc.new_pmgplist

new_pbelist = apbslibc.new_pbelist

new_atomforcelist = apbslibc.new_atomforcelist

double_array = apbslibc.double_array

int_array = apbslibc.int_array

Valist_load = apbslibc.Valist_load

def getPotentials(arg0,arg1,arg2,arg3):
    val = apbslibc.getPotentials(arg0.this,arg1.this,arg2.this,arg3)
    return val

get_entry = apbslibc.get_entry

set_entry = apbslibc.set_entry

make_Valist = apbslibc.make_Valist

def loadMolecules(arg0,arg1):
    val = apbslibc.loadMolecules(arg0.this,arg1)
    return val

def killMolecules(arg0,arg1):
    val = apbslibc.killMolecules(arg0.this,arg1)
    return val

def loadDielMaps(arg0,arg1,arg2,arg3):
    val = apbslibc.loadDielMaps(arg0.this,arg1,arg2,arg3)
    return val

def killDielMaps(arg0,arg1,arg2,arg3):
    val = apbslibc.killDielMaps(arg0.this,arg1,arg2,arg3)
    return val

def loadKappaMaps(arg0,arg1):
    val = apbslibc.loadKappaMaps(arg0.this,arg1)
    return val

def killKappaMaps(arg0,arg1):
    val = apbslibc.killKappaMaps(arg0.this,arg1)
    return val

def loadChargeMaps(arg0,arg1):
    val = apbslibc.loadChargeMaps(arg0.this,arg1)
    return val

def killChargeMaps(arg0,arg1):
    val = apbslibc.killChargeMaps(arg0.this,arg1)
    return val

def printPBEPARM(arg0):
    val = apbslibc.printPBEPARM(arg0.this)
    return val

def printMGPARM(arg0,arg1):
    val = apbslibc.printMGPARM(arg0.this,arg1)
    return val

def initMG(arg0,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13):
    val = apbslibc.initMG(arg0,arg1.this,arg2.this,arg3.this,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13)
    return val

def killMG(arg0,arg1,arg2,arg3):
    val = apbslibc.killMG(arg0.this,arg1,arg2,arg3)
    return val

def solveMG(arg0,arg1,arg2):
    val = apbslibc.solveMG(arg0.this,arg1.this,arg2)
    return val

def setPartMG(arg0,arg1,arg2):
    val = apbslibc.setPartMG(arg0.this,arg1.this,arg2.this)
    return val

def energyMG(arg0,arg1,arg2,arg3,arg4,arg5,arg6,arg7):
    val = apbslibc.energyMG(arg0.this,arg1,arg2.this,arg3,arg4,arg5,arg6,arg7)
    return val

def npenergyMG(arg0,arg1,arg2,arg3,arg4):
    val = apbslibc.npenergyMG(arg0.this,arg1,arg2.this,arg3,arg4)
    return val

killEnergy = apbslibc.killEnergy

def forceMG(arg0,arg1,arg2,arg3,arg4,arg5,arg6,arg7):
    val = apbslibc.forceMG(arg0.this,arg1.this,arg2.this,arg3.this,arg4.this,arg5,arg6,arg7)
    return val

def killForce(arg0,arg1,arg2,arg3):
    val = apbslibc.killForce(arg0.this,arg1.this,arg2,arg3)
    return val

def writedataMG(arg0,arg1,arg2,arg3):
    val = apbslibc.writedataMG(arg0,arg1.this,arg2.this,arg3.this)
    return val

def writematMG(arg0,arg1,arg2,arg3):
    val = apbslibc.writematMG(arg0,arg1.this,arg2.this,arg3.this)
    return val

def printEnergy(arg0,arg1,arg2,arg3):
    val = apbslibc.printEnergy(arg0.this,arg1.this,arg2,arg3)
    return val

def printForce(arg0,arg1,arg2,arg3,arg4):
    val = apbslibc.printForce(arg0.this,arg1.this,arg2,arg3,arg4)
    return val

startVio = apbslibc.startVio



#-------------- VARIABLE WRAPPERS ------------------

APBS_SWIG = apbslibc.APBS_SWIG
NPT_ENERGY = apbslibc.NPT_ENERGY
NPT_FORCE = apbslibc.NPT_FORCE
