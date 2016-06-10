#include "cell.hh"
#include "edge.hh"
#include "face.hh"
#include "vertex.hh"

#include <algorithm>//for std::copy
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>//for math operations sqrt and others

/**----------------------------------------------------------------------------------------------
 * Jacobian library : 
 * This calculates Energy and Jacobian of a given cell (tissue)
 * ----------------------------------------------------------------------------------------------
 */


double energyOfCell(Cell *cell){
	return 0.0;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
double energyOfFace(Face *face){
	/*
	//getting need properties from this face
	double area = face->getAreaOfFace();//area of face
	Cell *cell = face->getCell();//getting the cell on which this face lies
	double alpha = cell->getAlpha();
	double beta = cell->getBeta();
	double pressure = cell->getPressure();
	unsigned int faceid = face->getID();//getting the face id
	// ****************************************************************************************** //
	//setting the target form matrix of this face
	face->setTargetFormMatrix();
	double targetFormMatrix[2][2];
	std::copy(&face->targetFormMatrix[0][0],&face->targetFormMatrix+2*2,&targetFormMatrix[0][0]);
	double targetFormMatrixTraceSquared = face->getTraceSquaredTargetFormMatrix();
	// ****************************************************************************************** //
	//temporary variables to calculate centralised Coordinates
	double firstterm(0), secondterm(0), thirdterm(0);// to store the values of each term of energy
	// ****************************************************************************************** //
	// calculating the Form matrix of this face
    FaceEdgeIterator faceIterator(face);
    Edge *currentEdge;//to keep the current edge
    Vertex *currentVertex;//to keep the current vertex
    double ak,f1,f2,f3;//to store the values of functions
    double term11(0), term12(0), term22(0); // to store the values of matrix element
    double currentFormMatrix[2][2];// current Form Matrix
    //iterating through all the edges
    while ((currentEdge = faceIterator.next())!=0){
          currentVertex = currentEdge->Org();
          ak = currentVertex->getAk(faceid);
          f1 = currentVertex->getFunction1(faceid);
          f2 = currentVertex->getFunction2(faceid);
          f3 = currentVertex->getFunction3(faceid);
          //calculating the terms of the form matrix
          term11 += (1./12.)*ak*f1;
          term12 += (1./24.)*ak*f2;
          term22 += (1./12.)*ak*f3;
    }
    currentFormMatrix[0][0] = term11;
    currentFormMatrix[0][1] = term12;
    currentFormMatrix[1][0] = term12;
    currentFormMatrix[1][1] = term22;
    // ****************************************************************************************** //
    //First Term
    firstterm = (1./targetFormMatrixTraceSquared)*alpha*area*(pow(currentFormMatrix[0][0]-targetFormMatrix[0][0],2)+
    				  pow(currentFormMatrix[0][1]-targetFormMatrix[0][1],2)+
    				  pow(currentFormMatrix[1][0]-targetFormMatrix[1][0],2)+
    				  pow(currentFormMatrix[1][1]-targetFormMatrix[1][1],2));
    //Second Term
    secondterm = (1./targetFormMatrixTraceSquared)*beta*area*(pow(currentFormMatrix[0][0]-targetFormMatrix[0][0]+
    						currentFormMatrix[1][1]-targetFormMatrix[1][1],2));
    //Third term
    thirdterm = pressure*area;
    // ****************************************************************************************** //
    double energy = firstterm + secondterm - thirdterm;
    return energy;
    */
    return 0.;
}		


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
double jacobianOfCell(Cell *cell){
 		printf("number of faces : %u",cell->countFaces());
 		return 0.0;
 }
