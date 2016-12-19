#ifndef jacobianINCLUDED
#define jacobianINCLUDED

#include "cell.hh"
#include "edge.hh"
#include "face.hh"
#include "vertex.hh"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>//for math operations sqrt and others
//forward declaring the classes
class Cell;
class Edge;
class Vector;
class Face;
/**----------------------------------------------------------------------------------------------
 * Jacobian Class : 
 * This calculates Energy and Jacobian of a given cell (tissue)
 * ----------------------------------------------------------------------------------------------
 */
/**
 * Task : give the total energy of the faces in the cell
 * 		  or in biological terms : total energy of the tissue ( sum over all cells)
 * Input : cell that contains all the information of faces and vertex of tissue
 * Return : a double that is the total energy of the cell (tissue)
 */
class Jacobian
{
public: 
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
	// %%%%%%%%%%%%%%%%%% MAIN FUNCTIONS %%%%%%%%%%%%%%%%%%%% //
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
	double energyOfCell(Cell *cell);
	/**
 * Task : Calculate the energy of the given face
 * Input : pointer to the face
 * Result : a double that is the energy of this face
 */
	double energyOfFace(Face *face);
/**
 * Task : total jacobian of the tissue
 * @details [long description]
 * 
 * @param cell [description]
 */
	double jacobianOfCell(Cell *cell);
/**
 * Numerical Derivative Calculation of setProjectedCoordinate Fucntion
 * Calcualted at for each vertex
 * 
 */
	void setProjectedCoordinateDerivative(Face *face);
};


#endif 
/* #ifndef edgeINCLUDED */