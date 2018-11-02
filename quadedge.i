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
%include obj.i
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

%inline %{
double getEigenMatrix(Eigen::Matrix2d matrix, int i, int j){
	return matrix(i,j);
}	
%}

// Type map double *  (pointer) to array 

%typemap(out) double* {
	$result = PyList_New(3);
	for (int i = 0; i<3 ; ++ i){
		PyList_SetItem($result, i , PyFloat_FromDouble($1[i]));
	}
}

// Type map for double** to list of list in python
%typemap(out) double** {
	$result = PyList_New(2);
	for (int i = 0; i<2 ; ++ i){
		for (int j = 0; i<2 ; ++ j){
			l = PyList_New(2);
			PyList_SetItem(l, j, PyFloat_FromDouble($1[i][j]));
		}
		PyList_SetItem($result, i , l);
	}
}


