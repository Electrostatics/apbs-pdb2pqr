# This file was created automatically by SWIG.
import vgridc
class VgridPtr :
    def __init__(self,this):
        self.this = this
        self.thisown = 0
    def __del__(self):
        if self.thisown == 1 :
            vgridc.delete_Vgrid(self.this)
    def __setattr__(self,name,value):
        if name == "nx" :
            vgridc.Vgrid_nx_set(self.this,value)
            return
        if name == "ny" :
            vgridc.Vgrid_ny_set(self.this,value)
            return
        if name == "nz" :
            vgridc.Vgrid_nz_set(self.this,value)
            return
        if name == "hx" :
            vgridc.Vgrid_hx_set(self.this,value)
            return
        if name == "hy" :
            vgridc.Vgrid_hy_set(self.this,value)
            return
        if name == "hzed" :
            vgridc.Vgrid_hzed_set(self.this,value)
            return
        if name == "xmin" :
            vgridc.Vgrid_xmin_set(self.this,value)
            return
        if name == "ymin" :
            vgridc.Vgrid_ymin_set(self.this,value)
            return
        if name == "zmin" :
            vgridc.Vgrid_zmin_set(self.this,value)
            return
        if name == "data" :
            vgridc.Vgrid_data_set(self.this,value)
            return
        self.__dict__[name] = value
    def __getattr__(self,name):
        if name == "nx" : 
            return vgridc.Vgrid_nx_get(self.this)
        if name == "ny" : 
            return vgridc.Vgrid_ny_get(self.this)
        if name == "nz" : 
            return vgridc.Vgrid_nz_get(self.this)
        if name == "hx" : 
            return vgridc.Vgrid_hx_get(self.this)
        if name == "hy" : 
            return vgridc.Vgrid_hy_get(self.this)
        if name == "hzed" : 
            return vgridc.Vgrid_hzed_get(self.this)
        if name == "xmin" : 
            return vgridc.Vgrid_xmin_get(self.this)
        if name == "ymin" : 
            return vgridc.Vgrid_ymin_get(self.this)
        if name == "zmin" : 
            return vgridc.Vgrid_zmin_get(self.this)
        if name == "data" : 
            return vgridc.Vgrid_data_get(self.this)
        raise AttributeError,name
    def __repr__(self):
        return "<C Vgrid instance>"
class Vgrid(VgridPtr):
    def __init__(self) :
        self.this = vgridc.new_Vgrid()
        self.thisown = 1






#-------------- FUNCTION WRAPPERS ------------------

ptrcast = vgridc.ptrcast

ptrvalue = vgridc.ptrvalue

ptrset = vgridc.ptrset

ptrcreate = vgridc.ptrcreate

ptrfree = vgridc.ptrfree

ptradd = vgridc.ptradd

ptrmap = vgridc.ptrmap

null_array = vgridc.null_array

double_array = vgridc.double_array

get_entry = vgridc.get_entry

set_entry = vgridc.set_entry

def Vgrid_ctor(arg0,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9):
    val = vgridc.Vgrid_ctor(arg0,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9)
    val = VgridPtr(val)
    return val

def Vgrid_ctor2(arg0,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10):
    val = vgridc.Vgrid_ctor2(arg0.this,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10)
    return val

def Vgrid_value(arg0,arg1,arg2):
    val = vgridc.Vgrid_value(arg0.this,arg1,arg2)
    return val

Vgrid_dtor = vgridc.Vgrid_dtor

def Vgrid_dtor2(arg0):
    val = vgridc.Vgrid_dtor2(arg0.this)
    return val

def Vgrid_curvature(arg0,arg1,arg2,arg3):
    val = vgridc.Vgrid_curvature(arg0.this,arg1,arg2,arg3)
    return val

def Vgrid_gradient(arg0,arg1,arg2):
    val = vgridc.Vgrid_gradient(arg0.this,arg1,arg2)
    return val

def Vgrid_writeUHBD(arg0,arg1,arg2,arg3,arg4,arg5,arg6):
    val = vgridc.Vgrid_writeUHBD(arg0.this,arg1,arg2,arg3,arg4,arg5,arg6)
    return val

def Vgrid_writeDX(arg0,arg1,arg2,arg3,arg4,arg5,arg6):
    val = vgridc.Vgrid_writeDX(arg0.this,arg1,arg2,arg3,arg4,arg5,arg6)
    return val

def Vgrid_readDX(arg0,arg1,arg2,arg3,arg4):
    val = vgridc.Vgrid_readDX(arg0.this,arg1,arg2,arg3,arg4)
    return val

startVio = vgridc.startVio



#-------------- VARIABLE WRAPPERS ------------------

