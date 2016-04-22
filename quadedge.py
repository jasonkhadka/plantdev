# This file was automatically generated by SWIG (http://www.swig.org).
# Version 3.0.7
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.





from sys import version_info
if version_info >= (2, 6, 0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_quadedge', [dirname(__file__)])
        except ImportError:
            import _quadedge
            return _quadedge
        if fp is not None:
            try:
                _mod = imp.load_module('_quadedge', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _quadedge = swig_import_helper()
    del swig_import_helper
else:
    import _quadedge
del version_info
try:
    _swig_property = property
except NameError:
    pass  # Python < 2.2 doesn't have 'property'.


def _swig_setattr_nondynamic(self, class_type, name, value, static=1):
    if (name == "thisown"):
        return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'SwigPyObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name, None)
    if method:
        return method(self, value)
    if (not static):
        if _newclass:
            object.__setattr__(self, name, value)
        else:
            self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)


def _swig_setattr(self, class_type, name, value):
    return _swig_setattr_nondynamic(self, class_type, name, value, 0)


def _swig_getattr_nondynamic(self, class_type, name, static=1):
    if (name == "thisown"):
        return self.this.own()
    method = class_type.__swig_getmethods__.get(name, None)
    if method:
        return method(self)
    if (not static):
        return object.__getattr__(self, name)
    else:
        raise AttributeError(name)

def _swig_getattr(self, class_type, name):
    return _swig_getattr_nondynamic(self, class_type, name, 0)


def _swig_repr(self):
    try:
        strthis = "proxy of " + self.this.__repr__()
    except:
        strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

try:
    _object = object
    _newclass = 1
except AttributeError:
    class _object:
        pass
    _newclass = 0


class Cell(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Cell, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Cell, name)

    def __init__(self, *args, **kwargs):
        raise AttributeError("No constructor defined")
    __repr__ = _swig_repr
    __swig_getmethods__["make"] = lambda x: _quadedge.Cell_make
    if _newclass:
        make = staticmethod(_quadedge.Cell_make)
    __swig_getmethods__["makeTetrahedron"] = lambda x: _quadedge.Cell_makeTetrahedron
    if _newclass:
        makeTetrahedron = staticmethod(_quadedge.Cell_makeTetrahedron)
    __swig_getmethods__["kill"] = lambda x: _quadedge.Cell_kill
    if _newclass:
        kill = staticmethod(_quadedge.Cell_kill)

    def makeVertexEdge(self, vertex, left, right):
        return _quadedge.Cell_makeVertexEdge(self, vertex, left, right)

    def killVertexEdge(self, edge):
        return _quadedge.Cell_killVertexEdge(self, edge)

    def makeFaceEdge(self, face, org, dest):
        return _quadedge.Cell_makeFaceEdge(self, face, org, dest)

    def killFaceEdge(self, edge):
        return _quadedge.Cell_killFaceEdge(self, edge)

    def countVertices(self):
        return _quadedge.Cell_countVertices(self)

    def addVertex(self, vertex):
        return _quadedge.Cell_addVertex(self, vertex)

    def removeVertex(self, vertex):
        return _quadedge.Cell_removeVertex(self, vertex)

    def totalVertexSize(self):
        return _quadedge.Cell_totalVertexSize(self)

    def makeVertexID(self):
        return _quadedge.Cell_makeVertexID(self)

    def countFaces(self):
        return _quadedge.Cell_countFaces(self)

    def addFace(self, face):
        return _quadedge.Cell_addFace(self, face)

    def removeFace(self, face):
        return _quadedge.Cell_removeFace(self, face)

    def makeFaceID(self):
        return _quadedge.Cell_makeFaceID(self)
Cell_swigregister = _quadedge.Cell_swigregister
Cell_swigregister(Cell)

def Cell_make():
    return _quadedge.Cell_make()
Cell_make = _quadedge.Cell_make

def Cell_makeTetrahedron():
    return _quadedge.Cell_makeTetrahedron()
Cell_makeTetrahedron = _quadedge.Cell_makeTetrahedron

def Cell_kill(cell):
    return _quadedge.Cell_kill(cell)
Cell_kill = _quadedge.Cell_kill

class CellVertexIterator(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, CellVertexIterator, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, CellVertexIterator, name)
    __repr__ = _swig_repr

    def __init__(self, cell):
        this = _quadedge.new_CellVertexIterator(cell)
        try:
            self.this.append(this)
        except:
            self.this = this
    __swig_destroy__ = _quadedge.delete_CellVertexIterator
    __del__ = lambda self: None

    def next(self):
        return _quadedge.CellVertexIterator_next(self)
CellVertexIterator_swigregister = _quadedge.CellVertexIterator_swigregister
CellVertexIterator_swigregister(CellVertexIterator)

class CellFaceIterator(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, CellFaceIterator, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, CellFaceIterator, name)
    __repr__ = _swig_repr

    def __init__(self, cell):
        this = _quadedge.new_CellFaceIterator(cell)
        try:
            self.this.append(this)
        except:
            self.this = this
    __swig_destroy__ = _quadedge.delete_CellFaceIterator
    __del__ = lambda self: None

    def next(self):
        return _quadedge.CellFaceIterator_next(self)
CellFaceIterator_swigregister = _quadedge.CellFaceIterator_swigregister
CellFaceIterator_swigregister(CellFaceIterator)

class Face(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Face, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Face, name)

    def __init__(self, *args, **kwargs):
        raise AttributeError("No constructor defined")
    __repr__ = _swig_repr
    __swig_getmethods__["make"] = lambda x: _quadedge.Face_make
    if _newclass:
        make = staticmethod(_quadedge.Face_make)
    __swig_getmethods__["kill"] = lambda x: _quadedge.Face_kill
    if _newclass:
        kill = staticmethod(_quadedge.Face_kill)
    __swig_setmethods__["data"] = _quadedge.Face_data_set
    __swig_getmethods__["data"] = _quadedge.Face_data_get
    if _newclass:
        data = _swig_property(_quadedge.Face_data_get, _quadedge.Face_data_set)

    def getCell(self):
        return _quadedge.Face_getCell(self)

    def getID(self):
        return _quadedge.Face_getID(self)

    def setID(self, id):
        return _quadedge.Face_setID(self, id)

    def getEdge(self):
        return _quadedge.Face_getEdge(self)

    def addEdge(self, edge):
        return _quadedge.Face_addEdge(self, edge)

    def removeEdge(self, edge):
        return _quadedge.Face_removeEdge(self, edge)

    def getVertex(self):
        return _quadedge.Face_getVertex(self)

    def getIthVertex(self, position):
        return _quadedge.Face_getIthVertex(self, position)

    def addVertex(self, vertex):
        return _quadedge.Face_addVertex(self, vertex)

    def removeVertex(self, vertex):
        return _quadedge.Face_removeVertex(self, vertex)

    def countVertices(self):
        return _quadedge.Face_countVertices(self)
Face_swigregister = _quadedge.Face_swigregister
Face_swigregister(Face)

def Face_make(cell):
    return _quadedge.Face_make(cell)
Face_make = _quadedge.Face_make

def Face_kill(face):
    return _quadedge.Face_kill(face)
Face_kill = _quadedge.Face_kill

class FaceEdgeIterator(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, FaceEdgeIterator, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, FaceEdgeIterator, name)
    __repr__ = _swig_repr

    def __init__(self, face):
        this = _quadedge.new_FaceEdgeIterator(face)
        try:
            self.this.append(this)
        except:
            self.this = this
    __swig_destroy__ = _quadedge.delete_FaceEdgeIterator
    __del__ = lambda self: None

    def next(self):
        return _quadedge.FaceEdgeIterator_next(self)
FaceEdgeIterator_swigregister = _quadedge.FaceEdgeIterator_swigregister
FaceEdgeIterator_swigregister(FaceEdgeIterator)

class Edge(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Edge, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Edge, name)

    def __init__(self, *args, **kwargs):
        raise AttributeError("No constructor defined")
    __repr__ = _swig_repr
    __swig_getmethods__["make"] = lambda x: _quadedge.Edge_make
    if _newclass:
        make = staticmethod(_quadedge.Edge_make)
    __swig_getmethods__["kill"] = lambda x: _quadedge.Edge_kill
    if _newclass:
        kill = staticmethod(_quadedge.Edge_kill)
    __swig_getmethods__["splice"] = lambda x: _quadedge.Edge_splice
    if _newclass:
        splice = staticmethod(_quadedge.Edge_splice)
    __swig_setmethods__["data"] = _quadedge.Edge_data_set
    __swig_getmethods__["data"] = _quadedge.Edge_data_get
    if _newclass:
        data = _swig_property(_quadedge.Edge_data_get, _quadedge.Edge_data_set)

    def getID(self):
        return _quadedge.Edge_getID(self)

    def setID(self, id):
        return _quadedge.Edge_setID(self, id)

    def Org(self):
        return _quadedge.Edge_Org(self)

    def Dest(self):
        return _quadedge.Edge_Dest(self)

    def setOrg(self, org):
        return _quadedge.Edge_setOrg(self, org)

    def setDest(self, dest):
        return _quadedge.Edge_setDest(self, dest)

    def Left(self):
        return _quadedge.Edge_Left(self)

    def Right(self):
        return _quadedge.Edge_Right(self)

    def setLeft(self, left):
        return _quadedge.Edge_setLeft(self, left)

    def setRight(self, right):
        return _quadedge.Edge_setRight(self, right)

    def Rot(self):
        return _quadedge.Edge_Rot(self)

    def InvRot(self):
        return _quadedge.Edge_InvRot(self)

    def Sym(self):
        return _quadedge.Edge_Sym(self)

    def Onext(self):
        return _quadedge.Edge_Onext(self)

    def Oprev(self):
        return _quadedge.Edge_Oprev(self)

    def Dnext(self):
        return _quadedge.Edge_Dnext(self)

    def Dprev(self):
        return _quadedge.Edge_Dprev(self)

    def Lnext(self):
        return _quadedge.Edge_Lnext(self)

    def Lprev(self):
        return _quadedge.Edge_Lprev(self)

    def Rnext(self):
        return _quadedge.Edge_Rnext(self)

    def Rprev(self):
        return _quadedge.Edge_Rprev(self)
    __swig_setmethods__["Function1"] = _quadedge.Edge_Function1_set
    __swig_getmethods__["Function1"] = _quadedge.Edge_Function1_get
    if _newclass:
        Function1 = _swig_property(_quadedge.Edge_Function1_get, _quadedge.Edge_Function1_set)
    __swig_setmethods__["Function2"] = _quadedge.Edge_Function2_set
    __swig_getmethods__["Function2"] = _quadedge.Edge_Function2_get
    if _newclass:
        Function2 = _swig_property(_quadedge.Edge_Function2_get, _quadedge.Edge_Function2_set)
    __swig_setmethods__["Function3"] = _quadedge.Edge_Function3_set
    __swig_getmethods__["Function3"] = _quadedge.Edge_Function3_get
    if _newclass:
        Function3 = _swig_property(_quadedge.Edge_Function3_get, _quadedge.Edge_Function3_set)
    __swig_setmethods__["Ak"] = _quadedge.Edge_Ak_set
    __swig_getmethods__["Ak"] = _quadedge.Edge_Ak_get
    if _newclass:
        Ak = _swig_property(_quadedge.Edge_Ak_get, _quadedge.Edge_Ak_set)

    def setFunction1(self, arg2, arg3):
        return _quadedge.Edge_setFunction1(self, arg2, arg3)

    def setFunction2(self, arg2, arg3):
        return _quadedge.Edge_setFunction2(self, arg2, arg3)

    def setFunction3(self, arg2, arg3):
        return _quadedge.Edge_setFunction3(self, arg2, arg3)

    def setAk(self, arg2, arg3):
        return _quadedge.Edge_setAk(self, arg2, arg3)
Edge_swigregister = _quadedge.Edge_swigregister
Edge_swigregister(Edge)

def Edge_make():
    return _quadedge.Edge_make()
Edge_make = _quadedge.Edge_make

def Edge_kill(edge):
    return _quadedge.Edge_kill(edge)
Edge_kill = _quadedge.Edge_kill

def Edge_splice(a, b):
    return _quadedge.Edge_splice(a, b)
Edge_splice = _quadedge.Edge_splice


def objReadCell(name):
    return _quadedge.objReadCell(name)
objReadCell = _quadedge.objReadCell

def objWriteCell(cell, name):
    return _quadedge.objWriteCell(cell, name)
objWriteCell = _quadedge.objWriteCell

def objCloneCell(cell):
    return _quadedge.objCloneCell(cell)
objCloneCell = _quadedge.objCloneCell
class Vertex(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Vertex, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Vertex, name)

    def __init__(self, *args, **kwargs):
        raise AttributeError("No constructor defined")
    __repr__ = _swig_repr
    __swig_getmethods__["make"] = lambda x: _quadedge.Vertex_make
    if _newclass:
        make = staticmethod(_quadedge.Vertex_make)
    __swig_getmethods__["kill"] = lambda x: _quadedge.Vertex_kill
    if _newclass:
        kill = staticmethod(_quadedge.Vertex_kill)
    __swig_setmethods__["pos"] = _quadedge.Vertex_pos_set
    __swig_getmethods__["pos"] = _quadedge.Vertex_pos_get
    if _newclass:
        pos = _swig_property(_quadedge.Vertex_pos_get, _quadedge.Vertex_pos_set)
    __swig_setmethods__["data"] = _quadedge.Vertex_data_set
    __swig_getmethods__["data"] = _quadedge.Vertex_data_get
    if _newclass:
        data = _swig_property(_quadedge.Vertex_data_get, _quadedge.Vertex_data_set)

    def getCell(self):
        return _quadedge.Vertex_getCell(self)

    def getXcoordinate(self):
        return _quadedge.Vertex_getXcoordinate(self)

    def getYcoordinate(self):
        return _quadedge.Vertex_getYcoordinate(self)

    def getZcoordinate(self):
        return _quadedge.Vertex_getZcoordinate(self)

    def setXcoordinate(self, xcoord):
        return _quadedge.Vertex_setXcoordinate(self, xcoord)

    def setYcoordinate(self, ycoord):
        return _quadedge.Vertex_setYcoordinate(self, ycoord)

    def setZcoordinate(self, ycoord):
        return _quadedge.Vertex_setZcoordinate(self, ycoord)

    def getID(self):
        return _quadedge.Vertex_getID(self)

    def setID(self, id):
        return _quadedge.Vertex_setID(self, id)

    def insertProjectedXcoordinate(self, faceid, xcood):
        return _quadedge.Vertex_insertProjectedXcoordinate(self, faceid, xcood)

    def insertProjectedYcoordinate(self, faceid, ycood):
        return _quadedge.Vertex_insertProjectedYcoordinate(self, faceid, ycood)

    def getProjectedXcoordinate(self, faceid):
        return _quadedge.Vertex_getProjectedXcoordinate(self, faceid)

    def getProjectedYcoordinate(self, faceid):
        return _quadedge.Vertex_getProjectedYcoordinate(self, faceid)

    def getEdge(self):
        return _quadedge.Vertex_getEdge(self)

    def addEdge(self, edge):
        return _quadedge.Vertex_addEdge(self, edge)

    def removeEdge(self, edge):
        return _quadedge.Vertex_removeEdge(self, edge)
Vertex_swigregister = _quadedge.Vertex_swigregister
Vertex_swigregister(Vertex)

def Vertex_make(cell):
    return _quadedge.Vertex_make(cell)
Vertex_make = _quadedge.Vertex_make

def Vertex_kill(vertex):
    return _quadedge.Vertex_kill(vertex)
Vertex_kill = _quadedge.Vertex_kill

class VertexEdgeIterator(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, VertexEdgeIterator, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, VertexEdgeIterator, name)
    __repr__ = _swig_repr

    def __init__(self, vertex):
        this = _quadedge.new_VertexEdgeIterator(vertex)
        try:
            self.this.append(this)
        except:
            self.this = this
    __swig_destroy__ = _quadedge.delete_VertexEdgeIterator
    __del__ = lambda self: None

    def next(self):
        return _quadedge.VertexEdgeIterator_next(self)
VertexEdgeIterator_swigregister = _quadedge.VertexEdgeIterator_swigregister
VertexEdgeIterator_swigregister(VertexEdgeIterator)

# This file is compatible with both classic and new-style classes.


