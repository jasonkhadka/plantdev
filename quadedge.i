// File : quadedge.i to hold all the interface files together
%module quadedge

%include "numpy.i"

%include cell.i
%include list.i
%include face.i
%include edge.i
%include vertex.i
%include jacobian.i
%include derivative.i
%include "carrays.i"

%array_functions(double, doublearray)

//function to hopefully get the element values from a give array 
//swig documentation : http://www.iram.fr/~roche/code/python/files/SWIGDocumentation.pdf
%inline %{
double getCurrentFormMatrix(Face * face, int i, int j){
	return face->currentFormMatrix[i][j];
}	
%}