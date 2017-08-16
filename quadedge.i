// File : quadedge.i to hold all the interface files together
%module quadedge

//%include "numpy.i"

%include cell.i
%include list.i
%include face.i
%include edge.i
%include vertex.i
%include jacobian.i
%include derivative.i
%include "carrays.i"

%array_functions(double, doublearray)
//%apply (double* IN_ARRAY2, int DIM1, int DIM2) {(double* newValue, int i, int j)}
//function to hopefully get the element values from a give array 
//swig documentation : http://www.iram.fr/~roche/code/python/files/SWIGDocumentation.pdf
%inline %{
double getCurrentFormMatrix(Face * face, int i, int j){
	return face->currentFormMatrix[i][j];
}	
%}
//to get targetformmatrix
%inline %{
double getTargetFormMatrix(Face * face, int i, int j){
	return face->targetFormMatrix[i][j];
}	
%}
//to set targetformmatrix
%inline %{
double setTargetFormMatrix(Face * face, int i, int j,double newvalue){
	face->targetFormMatrix[i][j] = newvalue;
}	
%}
