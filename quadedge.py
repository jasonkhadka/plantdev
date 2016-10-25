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

    def getIthVertex(self, position):
        return _quadedge.Cell_getIthVertex(self, position)

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

    def setPressure(self, arg2):
        return _quadedge.Cell_setPressure(self, arg2)

    def getPressure(self):
        return _quadedge.Cell_getPressure(self)

    def setAlpha(self, arg2):
        return _quadedge.Cell_setAlpha(self, arg2)

    def getAlpha(self):
        return _quadedge.Cell_getAlpha(self)

    def setBeta(self, arg2):
        return _quadedge.Cell_setBeta(self, arg2)

    def getBeta(self):
        return _quadedge.Cell_getBeta(self)
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

    def setNormal(self, arg2):
        return _quadedge.Face_setNormal(self, arg2)

    def getNormal(self):
        return _quadedge.Face_getNormal(self)

    def setNormalTilde(self, arg2):
        return _quadedge.Face_setNormalTilde(self, arg2)

    def getNormalTilde(self):
        return _quadedge.Face_getNormalTilde(self)

    def setPivector(self, arg2):
        return _quadedge.Face_setPivector(self, arg2)

    def getPivector(self):
        return _quadedge.Face_getPivector(self)

    def setUnitx(self, arg2):
        return _quadedge.Face_setUnitx(self, arg2)

    def getUnitx(self):
        return _quadedge.Face_getUnitx(self)

    def setUnity(self, arg2):
        return _quadedge.Face_setUnity(self, arg2)

    def getUnity(self):
        return _quadedge.Face_getUnity(self)

    def setUnitz(self, arg2):
        return _quadedge.Face_setUnitz(self, arg2)

    def getUnitz(self):
        return _quadedge.Face_getUnitz(self)
    __swig_setmethods__["targetFormMatrix"] = _quadedge.Face_targetFormMatrix_set
    __swig_getmethods__["targetFormMatrix"] = _quadedge.Face_targetFormMatrix_get
    if _newclass:
        targetFormMatrix = _swig_property(_quadedge.Face_targetFormMatrix_get, _quadedge.Face_targetFormMatrix_set)
    __swig_setmethods__["traceSquaredTargetFormMatrix"] = _quadedge.Face_traceSquaredTargetFormMatrix_set
    __swig_getmethods__["traceSquaredTargetFormMatrix"] = _quadedge.Face_traceSquaredTargetFormMatrix_get
    if _newclass:
        traceSquaredTargetFormMatrix = _swig_property(_quadedge.Face_traceSquaredTargetFormMatrix_get, _quadedge.Face_traceSquaredTargetFormMatrix_set)
    __swig_setmethods__["mu1"] = _quadedge.Face_mu1_set
    __swig_getmethods__["mu1"] = _quadedge.Face_mu1_get
    if _newclass:
        mu1 = _swig_property(_quadedge.Face_mu1_get, _quadedge.Face_mu1_set)
    __swig_setmethods__["mu2"] = _quadedge.Face_mu2_set
    __swig_getmethods__["mu2"] = _quadedge.Face_mu2_get
    if _newclass:
        mu2 = _swig_property(_quadedge.Face_mu2_get, _quadedge.Face_mu2_set)
    __swig_setmethods__["mu3"] = _quadedge.Face_mu3_set
    __swig_getmethods__["mu3"] = _quadedge.Face_mu3_get
    if _newclass:
        mu3 = _swig_property(_quadedge.Face_mu3_get, _quadedge.Face_mu3_set)
    __swig_setmethods__["mu4"] = _quadedge.Face_mu4_set
    __swig_getmethods__["mu4"] = _quadedge.Face_mu4_get
    if _newclass:
        mu4 = _swig_property(_quadedge.Face_mu4_get, _quadedge.Face_mu4_set)

    def setMu(self):
        return _quadedge.Face_setMu(self)

    def getMu1(self):
        return _quadedge.Face_getMu1(self)

    def getMu2(self):
        return _quadedge.Face_getMu2(self)

    def getMu3(self):
        return _quadedge.Face_getMu3(self)

    def getMu4(self):
        return _quadedge.Face_getMu4(self)

    def setAngleOfTilt(self, arg2):
        return _quadedge.Face_setAngleOfTilt(self, arg2)

    def getAngleOfTilt(self):
        return _quadedge.Face_getAngleOfTilt(self)

    def setTargetFormMatrix(self):
        return _quadedge.Face_setTargetFormMatrix(self)

    def setTraceSquaredTargetFormMatrix(self):
        return _quadedge.Face_setTraceSquaredTargetFormMatrix(self)

    def getTraceSquaredTargetFormMatrix(self):
        return _quadedge.Face_getTraceSquaredTargetFormMatrix(self)

    def setCentralisedCoordinate(self, arg2, arg3, arg4):
        return _quadedge.Face_setCentralisedCoordinate(self, arg2, arg3, arg4)

    def setProjectedCoordinate(self):
        return _quadedge.Face_setProjectedCoordinate(self)

    def getXCentralised(self):
        return _quadedge.Face_getXCentralised(self)

    def getYCentralised(self):
        return _quadedge.Face_getYCentralised(self)

    def getZCentralised(self):
        return _quadedge.Face_getZCentralised(self)

    def setAreaOfFace(self):
        return _quadedge.Face_setAreaOfFace(self)

    def getAreaOfFace(self):
        return _quadedge.Face_getAreaOfFace(self)

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

    def setZcoordinate(self, zcoord):
        return _quadedge.Vertex_setZcoordinate(self, zcoord)

    def setAlphaBetaGamma(self):
        return _quadedge.Vertex_setAlphaBetaGamma(self)

    def setAlpha(self, arg2, arg3):
        return _quadedge.Vertex_setAlpha(self, arg2, arg3)

    def setBeta(self, arg2, arg3):
        return _quadedge.Vertex_setBeta(self, arg2, arg3)

    def setGamma(self, arg2, arg3):
        return _quadedge.Vertex_setGamma(self, arg2, arg3)

    def getAlpha(self, id):
        return _quadedge.Vertex_getAlpha(self, id)

    def getBeta(self, id):
        return _quadedge.Vertex_getBeta(self, id)

    def getGamma(self, id):
        return _quadedge.Vertex_getGamma(self, id)

    def getID(self):
        return _quadedge.Vertex_getID(self)

    def setID(self, id):
        return _quadedge.Vertex_setID(self, id)

    def insertProjectedXcoordinate(self, faceid, xcood):
        return _quadedge.Vertex_insertProjectedXcoordinate(self, faceid, xcood)

    def insertProjectedYcoordinate(self, faceid, ycood):
        return _quadedge.Vertex_insertProjectedYcoordinate(self, faceid, ycood)

    def insertProjectedZcoordinate(self, faceid, zcood):
        return _quadedge.Vertex_insertProjectedZcoordinate(self, faceid, zcood)

    def insertNonCentralisedProjectedXcoordinate(self, faceid, xcood):
        return _quadedge.Vertex_insertNonCentralisedProjectedXcoordinate(self, faceid, xcood)

    def insertNonCentralisedProjectedYcoordinate(self, faceid, ycood):
        return _quadedge.Vertex_insertNonCentralisedProjectedYcoordinate(self, faceid, ycood)

    def insertNonCentralisedProjectedZcoordinate(self, faceid, zcood):
        return _quadedge.Vertex_insertNonCentralisedProjectedZcoordinate(self, faceid, zcood)

    def getProjectedXcoordinate(self, faceid):
        return _quadedge.Vertex_getProjectedXcoordinate(self, faceid)

    def getProjectedYcoordinate(self, faceid):
        return _quadedge.Vertex_getProjectedYcoordinate(self, faceid)

    def getProjectedZcoordinate(self, faceid):
        return _quadedge.Vertex_getProjectedZcoordinate(self, faceid)

    def getNonCentralisedProjectedXcoordinate(self, faceid):
        return _quadedge.Vertex_getNonCentralisedProjectedXcoordinate(self, faceid)

    def getNonCentralisedProjectedYcoordinate(self, faceid):
        return _quadedge.Vertex_getNonCentralisedProjectedYcoordinate(self, faceid)

    def getNonCentralisedProjectedZcoordinate(self, faceid):
        return _quadedge.Vertex_getNonCentralisedProjectedZcoordinate(self, faceid)

    def setAreaDerivative(self):
        return _quadedge.Vertex_setAreaDerivative(self)

    def getAreaXDerivative(self, faceid):
        return _quadedge.Vertex_getAreaXDerivative(self, faceid)

    def getAreaYDerivative(self, faceid):
        return _quadedge.Vertex_getAreaYDerivative(self, faceid)

    def setAkDerivative(self):
        return _quadedge.Vertex_setAkDerivative(self)

    def getAkXDerivative(self, faceid):
        return _quadedge.Vertex_getAkXDerivative(self, faceid)

    def getAkYDerivative(self, faceid):
        return _quadedge.Vertex_getAkYDerivative(self, faceid)

    def setMuXDerivative(self):
        return _quadedge.Vertex_setMuXDerivative(self)

    def setMuSquaredXDerivative(self):
        return _quadedge.Vertex_setMuSquaredXDerivative(self)

    def setFirstTermXDerivative(self):
        return _quadedge.Vertex_setFirstTermXDerivative(self)

    def setSecondTermXDerivative(self):
        return _quadedge.Vertex_setSecondTermXDerivative(self)

    def setThirdTermXDerivative(self):
        return _quadedge.Vertex_setThirdTermXDerivative(self)

    def getMu1XDerivative(self, arg2):
        return _quadedge.Vertex_getMu1XDerivative(self, arg2)

    def getMu4XDerivative(self, arg2):
        return _quadedge.Vertex_getMu4XDerivative(self, arg2)

    def getFirstTermXDerivative(self):
        return _quadedge.Vertex_getFirstTermXDerivative(self)

    def getSecondTermXDerivative(self):
        return _quadedge.Vertex_getSecondTermXDerivative(self)

    def getThirdTermXDerivative(self):
        return _quadedge.Vertex_getThirdTermXDerivative(self)

    def getMu1SquaredXDerivative(self, facid):
        return _quadedge.Vertex_getMu1SquaredXDerivative(self, facid)

    def getMu2SquaredXDerivative(self, faceid):
        return _quadedge.Vertex_getMu2SquaredXDerivative(self, faceid)

    def getMu3SquaredXDerivative(self, faceid):
        return _quadedge.Vertex_getMu3SquaredXDerivative(self, faceid)

    def getMu4SquaredXDerivative(self, faceid):
        return _quadedge.Vertex_getMu4SquaredXDerivative(self, faceid)

    def setMuYDerivative(self):
        return _quadedge.Vertex_setMuYDerivative(self)

    def setMuSquaredYDerivative(self):
        return _quadedge.Vertex_setMuSquaredYDerivative(self)

    def setFirstTermYDerivative(self):
        return _quadedge.Vertex_setFirstTermYDerivative(self)

    def setSecondTermYDerivative(self):
        return _quadedge.Vertex_setSecondTermYDerivative(self)

    def setThirdTermYDerivative(self):
        return _quadedge.Vertex_setThirdTermYDerivative(self)

    def getFirstTermYDerivative(self):
        return _quadedge.Vertex_getFirstTermYDerivative(self)

    def getSecondTermYDerivative(self):
        return _quadedge.Vertex_getSecondTermYDerivative(self)

    def getThirdTermYDerivative(self):
        return _quadedge.Vertex_getThirdTermYDerivative(self)

    def getMu1YDerivative(self, arg2):
        return _quadedge.Vertex_getMu1YDerivative(self, arg2)

    def getMu4YDerivative(self, arg2):
        return _quadedge.Vertex_getMu4YDerivative(self, arg2)

    def getMu1SquaredYDerivative(self, facid):
        return _quadedge.Vertex_getMu1SquaredYDerivative(self, facid)

    def getMu2SquaredYDerivative(self, faceid):
        return _quadedge.Vertex_getMu2SquaredYDerivative(self, faceid)

    def getMu3SquaredYDerivative(self, faceid):
        return _quadedge.Vertex_getMu3SquaredYDerivative(self, faceid)

    def getMu4SquaredYDerivative(self, faceid):
        return _quadedge.Vertex_getMu4SquaredYDerivative(self, faceid)

    def setFunctions(self):
        return _quadedge.Vertex_setFunctions(self)

    def setFunction1(self):
        return _quadedge.Vertex_setFunction1(self)

    def setFunction2(self):
        return _quadedge.Vertex_setFunction2(self)

    def setFunction3(self):
        return _quadedge.Vertex_setFunction3(self)

    def setAk(self):
        return _quadedge.Vertex_setAk(self)

    def getAk(self, faceid):
        return _quadedge.Vertex_getAk(self, faceid)

    def getFunction1(self, faceid):
        return _quadedge.Vertex_getFunction1(self, faceid)

    def getFunction2(self, faceid):
        return _quadedge.Vertex_getFunction2(self, faceid)

    def getFunction3(self, faceid):
        return _quadedge.Vertex_getFunction3(self, faceid)

    def getMu1(self, faceid):
        return _quadedge.Vertex_getMu1(self, faceid)

    def getMu2(self, faceid):
        return _quadedge.Vertex_getMu2(self, faceid)

    def getMu3(self, faceid):
        return _quadedge.Vertex_getMu3(self, faceid)

    def getMu4(self, faceid):
        return _quadedge.Vertex_getMu4(self, faceid)

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


def energyOfCell(cell):
    return _quadedge.energyOfCell(cell)
energyOfCell = _quadedge.energyOfCell

def energyOfFace(face):
    return _quadedge.energyOfFace(face)
energyOfFace = _quadedge.energyOfFace

def jacobianOfCell(cell):
    return _quadedge.jacobianOfCell(cell)
jacobianOfCell = _quadedge.jacobianOfCell

def setProjectedCoordinateDerivative(face):
    return _quadedge.setProjectedCoordinateDerivative(face)
setProjectedCoordinateDerivative = _quadedge.setProjectedCoordinateDerivative
class CentralisedDerivative(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, CentralisedDerivative, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, CentralisedDerivative, name)
    __repr__ = _swig_repr

    def xtildeXDerivative(self, first, second, face):
        return _quadedge.CentralisedDerivative_xtildeXDerivative(self, first, second, face)

    def ytildeXDerivative(self, first, second, face):
        return _quadedge.CentralisedDerivative_ytildeXDerivative(self, first, second, face)

    def xiXDerivative(self, first, second, face):
        return _quadedge.CentralisedDerivative_xiXDerivative(self, first, second, face)

    def ex1XDerivative(self, first, second, face):
        return _quadedge.CentralisedDerivative_ex1XDerivative(self, first, second, face)

    def ex2XDerivative(self, first, second, face):
        return _quadedge.CentralisedDerivative_ex2XDerivative(self, first, second, face)

    def ex3XDerivative(self, first, second, face):
        return _quadedge.CentralisedDerivative_ex3XDerivative(self, first, second, face)

    def ey1XDerivative(self, first, second, face):
        return _quadedge.CentralisedDerivative_ey1XDerivative(self, first, second, face)

    def ey2XDerivative(self, first, second, face):
        return _quadedge.CentralisedDerivative_ey2XDerivative(self, first, second, face)

    def ey3XDerivative(self, first, second, face):
        return _quadedge.CentralisedDerivative_ey3XDerivative(self, first, second, face)

    def deltafunction(self, arg2, arg3):
        return _quadedge.CentralisedDerivative_deltafunction(self, arg2, arg3)

    def alphaXDerivative(self, arg2, arg3, arg4):
        return _quadedge.CentralisedDerivative_alphaXDerivative(self, arg2, arg3, arg4)

    def betaXDerivative(self, arg2, arg3, arg4):
        return _quadedge.CentralisedDerivative_betaXDerivative(self, arg2, arg3, arg4)

    def gammaXDerivative(self, arg2, arg3, arg4):
        return _quadedge.CentralisedDerivative_gammaXDerivative(self, arg2, arg3, arg4)

    def cpinormXDerivative(self, first, second, face):
        return _quadedge.CentralisedDerivative_cpinormXDerivative(self, first, second, face)

    def pixXDerivative(self, arg2, arg3, arg4):
        return _quadedge.CentralisedDerivative_pixXDerivative(self, arg2, arg3, arg4)

    def piyXDerivative(self, arg2, arg3, arg4):
        return _quadedge.CentralisedDerivative_piyXDerivative(self, arg2, arg3, arg4)

    def pizXDerivative(self, arg2, arg3, arg4):
        return _quadedge.CentralisedDerivative_pizXDerivative(self, arg2, arg3, arg4)

    def ncxXDerivative(self, arg2, arg3, arg4):
        return _quadedge.CentralisedDerivative_ncxXDerivative(self, arg2, arg3, arg4)

    def ncyXDerivative(self, arg2, arg3, arg4):
        return _quadedge.CentralisedDerivative_ncyXDerivative(self, arg2, arg3, arg4)

    def nczXDerivative(self, arg2, arg3, arg4):
        return _quadedge.CentralisedDerivative_nczXDerivative(self, arg2, arg3, arg4)

    def areatotalXDerivative(self, arg2, arg3, arg4):
        return _quadedge.CentralisedDerivative_areatotalXDerivative(self, arg2, arg3, arg4)

    def nctildenormXDerivative(self, arg2, arg3, arg4):
        return _quadedge.CentralisedDerivative_nctildenormXDerivative(self, arg2, arg3, arg4)

    def ncxtildeXDerivative(self, arg2, arg3, arg4):
        return _quadedge.CentralisedDerivative_ncxtildeXDerivative(self, arg2, arg3, arg4)

    def ncytildeXDerivative(self, arg2, arg3, arg4):
        return _quadedge.CentralisedDerivative_ncytildeXDerivative(self, arg2, arg3, arg4)

    def ncztildeXDerivative(self, arg2, arg3, arg4):
        return _quadedge.CentralisedDerivative_ncztildeXDerivative(self, arg2, arg3, arg4)

    def numericalXtildeXDerivative(self, arg2, arg3, arg4, arg5):
        return _quadedge.CentralisedDerivative_numericalXtildeXDerivative(self, arg2, arg3, arg4, arg5)

    def numericalYtildeXDerivative(self, arg2, arg3, arg4, arg5):
        return _quadedge.CentralisedDerivative_numericalYtildeXDerivative(self, arg2, arg3, arg4, arg5)

    def __init__(self):
        this = _quadedge.new_CentralisedDerivative()
        try:
            self.this.append(this)
        except:
            self.this = this
    __swig_destroy__ = _quadedge.delete_CentralisedDerivative
    __del__ = lambda self: None
CentralisedDerivative_swigregister = _quadedge.CentralisedDerivative_swigregister
CentralisedDerivative_swigregister(CentralisedDerivative)


def new_doublearray(nelements):
    return _quadedge.new_doublearray(nelements)
new_doublearray = _quadedge.new_doublearray

def delete_doublearray(ary):
    return _quadedge.delete_doublearray(ary)
delete_doublearray = _quadedge.delete_doublearray

def doublearray_getitem(ary, index):
    return _quadedge.doublearray_getitem(ary, index)
doublearray_getitem = _quadedge.doublearray_getitem

def doublearray_setitem(ary, index, value):
    return _quadedge.doublearray_setitem(ary, index, value)
doublearray_setitem = _quadedge.doublearray_setitem
# This file is compatible with both classic and new-style classes.


