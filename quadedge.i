// File : quadedge.i to hold all the interface files together
%module quadedge
%include cell.i
%include list.i
%include face.i
%include edge.i
%include obj.i
%include array.i
%include vertex.i
%include jacobian.i
%include derivative.i
%include "carrays.i"
%array_functions(double, doublearray)