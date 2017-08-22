
#include "cell.hh"
#include "edge.hh"
#include "vertex.hh"
#include "face.hh"
#include "derivative.hh"
#include <math.h>

// ************************************************************************************ //
// ************************************************************************************ //

// -----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%//
// -----%-----%-----%----Public Members Only-----%-----%-----%-----%-----%-----%-----%//
// -----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%//
/*
This is a temprorary work and earlier work of 
double CentralisedDerivative::energyXDerivative(Cell * cell, Vertex * first){
	// iterating the edges arround the vertex first
	VertexEdgeIterator vertedgeiterator(first);
	Edge * currentEdge;
	Face * currentFace;
	Edge * tempEdge;
	Vertex * tempvertex;
	double alpha = cell->getAlpha();
	double beta = cell->getBeta();
	double pressure = cell->getPressure();
	double energyxtildederivative(0.), energyytildederivative(0.);
	double xtildexderivative(0.), ytildexderivative(0.);
	double energyxderivativeFace(0.), totalenergyxderivative(0.);
	// iterating through the vertex
	while((currentEdge = vertedgeiterator.next())!=0){
				currentFace = currentEdge->Left();
				// iterating the face to get the vertex arround
				FaceEdgeIterator faceedgeiterator(currentFace);
				while ((tempEdge = faceedgeiterator.next())!= 0 ){
						tempvertex = tempEdge->Org();//grabbing origin of the edge
						//now calculating the derivative
						energyxtildederivative = alpha*tempvertex->getFirstTermXDerivative()+
												beta*tempvertex->getSecondTermXDerivative()+
												pressure*tempvertex->getThirdTermXDerivative();

						energyytildederivative = alpha*tempvertex->getFirstTermYDerivative()+
												beta*tempvertex->getSecondTermYDerivative()+
												pressure*tempvertex->getThirdTermYDerivative();
						xtildexderivative = xtildeXDerivative(tempvertex,first,currentFace);
						ytildexderivative = ytildeXDerivative(tempvertex,first,currentFace);
						//summing up the derivatives 
						energyxderivativeFace+= energyxtildederivative*xtildexderivative + 
												 energyytildederivative*ytildexderivative;
						}
				totalenergyxderivative += energyxderivativeFace;
		}
	return totalenergyxderivative;
}
*/
//****************************************************************************************** //
//   Numerical Energy X derivative with respect to change in x coordinate of this vertex
double CentralisedDerivative::numericalEnergyXDerivative(Vertex* vertex){
	// forward difference
	Cell * cell = vertex->getCell();
	double currentEnergy = cell->getDerivativeEnergyCartesianVolume();
	//*** Now forward propagating the x vertex and the energy *** //
	double currentXCoordinate = vertex->getXcoordinate();
	// defining the stepsize for finite difference calculation
	double stepsize;
	if (currentXCoordinate != 0){
		stepsize = (cell->getSqrtEpsilon())*currentXCoordinate;
	}
	else{
		stepsize = cell->getSqrtEpsilon();
	}
	//now caculating the new energy by increasing the X-coordinate
	vertex->setXcoordinate(currentXCoordinate+stepsize);
	cell->setParameters();
	double plusEnergy = cell->getDerivativeEnergyCartesianVolume();
	//now caculating the new energy by decreasing the X-coordinate
	vertex->setXcoordinate(currentXCoordinate-stepsize);
	cell->setParameters();
	double minusEnergy = cell->getDerivativeEnergyCartesianVolume();
	//Resetting the coordinate and other values
	vertex->setXcoordinate(currentXCoordinate);
	cell->setParameters();
	//calculating the derivative
	double derivative = (plusEnergy-minusEnergy)/(2*stepsize);
	//to check if the values are reset to previous values
	/*std::cout<< "Vertex ID : "<< vertex->getID() 
	<<"\n Start Energy" << currentEnergy 
	<<"\n Energy Resting back to Normal : " << cell->getDerivativeEnergyCartesianVolume() 
	<< "\n plus Energy : " << plusEnergy
	<< "\n minus Energy : "<< minusEnergy
	<< "\n   difference : " <<plusEnergy-minusEnergy
	<<"\n Step size : "<< stepsize 
	<<"\n derivative :" <<derivative
	<<"\n ========================================"<< std::endl;
	*/
	// returning the derivative value
	return derivative;
}
//****************************************************************************************** //
//   Numerical Energy Y derivative with respect to change in y coordinate of this vertex
double CentralisedDerivative::numericalEnergyYDerivative(Vertex* vertex){
	// forward difference
	Cell * cell = vertex->getCell();
	double currentEnergy = cell->getDerivativeEnergyCartesianVolume();
	//*** Now forward propagating the x vertex and the energy *** //
	double currentYCoordinate = vertex->getYcoordinate();
	// defining the stepsize for finite difference calculation
	double stepsize;
	if (currentYCoordinate != 0){
		stepsize = (cell->getSqrtEpsilon())*currentYCoordinate;
	}
	else{
		stepsize = cell->getSqrtEpsilon();
	}
	//now caculating the new energy by increasing the Y-coordinate
	vertex->setYcoordinate(currentYCoordinate+stepsize);
	cell->setParameters();
	double plusEnergy = cell->getDerivativeEnergyCartesianVolume();
	//now caculating the new energy by decreasing the Y-coordinate
	vertex->setYcoordinate(currentYCoordinate-stepsize);
	cell->setParameters();
	double minusEnergy = cell->getDerivativeEnergyCartesianVolume();
	//Resetting the coordinate and other values
	vertex->setYcoordinate(currentYCoordinate);
	cell->setParameters();
	//calculating the derivative
	double derivative = (plusEnergy-minusEnergy)/(2*stepsize);
	//to check if the values are reset to previous values
	/*
	std::cout
	<<"Y - Derivative : "
	<<"\n Vertex ID : "<< vertex->getID() 
	<<"\n Start Energy" << currentEnergy 
	<<"\n Energy Resting back to Normal : " << cell->getDerivativeEnergyCartesianVolume() 
	<< "\n plus Energy : " << plusEnergy
	<< "\n minus Energy : "<< minusEnergy
	<< "\n   difference : " <<plusEnergy-minusEnergy
	<<"\n Step size : "<< stepsize 
	<<"\n derivative :" <<derivative
	<<"\n ========================================"<< std::endl;
	*/
	// returning the derivative value
	return derivative;
}
//****************************************************************************************** //
//   Numerical Energy Z derivative with respect to change in z coordinate of this vertex
double CentralisedDerivative::numericalEnergyZDerivative(Vertex* vertex){
	// forward difference
	Cell * cell = vertex->getCell();
	double currentEnergy = cell->getDerivativeEnergyCartesianVolume();
	//*** Now forward propagating the x vertex and the energy *** //
	double currentZCoordinate = vertex->getZcoordinate();
	// defining the stepsize for finite difference calculation
	double stepsize;
	if (currentZCoordinate != 0){
		stepsize = (cell->getSqrtEpsilon())*currentZCoordinate;
	}
	else{
		stepsize = cell->getSqrtEpsilon();
	}
	//now caculating the new energy by increasing the Z-coordinate
	vertex->setZcoordinate(currentZCoordinate+stepsize);
	cell->setParameters();
	double plusEnergy = cell->getDerivativeEnergyCartesianVolume();
	//now caculating the new energy by decreasing the Z-coordinate
	vertex->setZcoordinate(currentZCoordinate-stepsize);
	cell->setParameters();
	double minusEnergy = cell->getDerivativeEnergyCartesianVolume();
	//Resetting the coordinate and other values
	vertex->setZcoordinate(currentZCoordinate);
	cell->setParameters();
	//calculating the derivative
	double derivative = (plusEnergy-minusEnergy)/(2*stepsize);
	//to check if the values are reset to previous values
	/*std::cout<< "Vertex ID : "<< vertex->getID() 
	<<"\n Start Energy" << currentEnergy 
	<<"\n Energy Resting back to Normal : " << cell->getDerivativeEnergyCartesianVolume() 
	<< "\n plus Energy : " << plusEnergy
	<< "\n minus Energy : "<< minusEnergy
	<< "\n   difference : " <<plusEnergy-minusEnergy
	<<"\n Step size : "<< stepsize 
	<<"\n derivative :" <<derivative
	<<"\n ========================================"<< std::endl;
	*/
	// returning the derivative value
	return derivative;
}
//****************************************************************************************** //
//****************************************************************************************** //
double energyXDerivative(Cell * cell, Vertex * vertex){
	// iterating through the edges of the vertex
	VertexEdgeIterator vertedgeiterator(vertex);
	Edge * currentEdge;
	Face * currentFace;
	Edge * tempEdge;
	Vertex * tempVertex;
	while((currentEdge = vertedgeiterator.next())!=0){
			currentFace = currentEdge->Left(); // grabing left face of the edge
			FaceEdgeIterator faceedgeiterator(currentFace);
			while((tempEdge = faceedgeiterator.next())!= 0){
					tempVertex = tempEdge->Dest();
					
			}
	}
	return 0.;
}
// --------------------------------------------------------------------------------------------------- //
double CentralisedDerivative::xtildeXDerivative(Vertex *first, Vertex *second, Face *face){
	// gathering terms //
	if (face->getID() == 1){
		return 0.0;
	}
	else{
	double* unitxvector = face->getUnitx();
	double coordinates[3] = {first->getXcoordinate(),
							first->getYcoordinate(),
							first->getZcoordinate()};
	double centroid[3] ={face->getXCentralised(),
						face->getYCentralised(),
						face->getZCentralised()};
	//double * unitx = face->getUnitx();
	// calculating derivatives //
	double xiderivative = xiXDerivative(first, second, face);
	double ex1derivative =  ex1XDerivative(first, second, face);
	double ex2derivative = ex2XDerivative(first, second, face);
	double ex3derivative =  ex3XDerivative(first, second, face);
	// final calcuation //
	double derivativeValue = unitxvector[0]*xiderivative + coordinates[0]*ex1derivative +
								coordinates[1]*ex2derivative +
								coordinates[2]*ex3derivative
								- unitxvector[0]*1./6. -centroid[0]*ex1derivative- centroid[1]*ex2derivative
								- centroid[2]*ex3derivative;
	// returning value //
	return derivativeValue;
	}
	
}

// --------------------------------------------------------------------------------------------------- //
double CentralisedDerivative::ytildeXDerivative(Vertex *first, Vertex *second, Face *face){
	// gathering terms //
	if (face->getID() == 1){
		return 0.0;
	}
	else{
	//double* unitxvector = face->getUnitx();
	double coordinates[3] = {first->getXcoordinate(),
							first->getYcoordinate(),
							first->getZcoordinate()};
	double centroid[3] ={face->getXCentralised(),
						face->getYCentralised(),
						face->getZCentralised()};
	double * unityvector = face->getUnity();
	// calculating derivatives //
	double xiderivative = xiXDerivative(first, second, face);
	double ey1derivative = ey1XDerivative(first, second, face);
	double ey2derivative = ey2XDerivative(first, second, face);
	double ey3derivative = ey3XDerivative(first, second, face);
	// final calcuation //
	double derivativeValue = (coordinates[0]-centroid[0])*ey1derivative
							 +(coordinates[1]-centroid[1])*ey2derivative
							 +(coordinates[2]-centroid[2])*ey3derivative
							 +unityvector[0]*(xiderivative-1./6.);
	// returning value //
	//std::cout<<"Unit Y : "<<unityvector[0]<<unityvector[1]<<unityvector[3]<<std::endl;
	return derivativeValue;
	}
	
}

// -----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%//
// -----%-----%-----%----Private Members Only-----%-----%-----%-----%-----%-----%-----%//
// -----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%//




// --------------------------------------------------------------------------------- //
double CentralisedDerivative::xiXDerivative(Vertex *first, Vertex *second, Face *face){
	return deltafunction(first->getID(),second->getID());
}
// --------------------------------------------------------------------------------- //

double CentralisedDerivative::ex1XDerivative(Vertex *first, Vertex *second, Face *face){
	//gathering terms // 
	double *pntpivector = face->getPivector();
	double *pntnormal = face->getNormal();
	double cPiNorm = sqrt(pow(pntpivector[0]-face->getXCentralised(),2)+
							pow(pntpivector[1]-face->getYCentralised(),2) + 
							pow(pntpivector[2]-face->getZCentralised(),2));
	double ncxderivative = ncxXDerivative(first, second, face); //calculating the ncxderivative
	double cPiNormxderivative = cpinormXDerivative(first,second,face); //calculating the cpinorm X derivative
	//calculating the derivative value // 
	double derivativeValue = (cPiNorm*(-2.*pntnormal[0]*ncxderivative)-(1-pow(pntnormal[0],2))*cPiNormxderivative)/pow(cPiNorm,2.);
	return derivativeValue;
}
// --------------------------------------------------------------------------------- //
double CentralisedDerivative::ex2XDerivative(Vertex* first, Vertex* second, Face* face){
	// gathering the terms //
	double *pntpivector = face->getPivector();
	double *pntnormal = face->getNormal();
	double cpivector[3] = {pntpivector[0]-face->getXCentralised(),pntpivector[1]-face->getYCentralised(),pntpivector[2]-face->getZCentralised()};
	double cpinorm = sqrt(pow(cpivector[0],2)+
							pow(cpivector[1],2) + 
							pow(cpivector[2],2));
	double piyderivative = piyXDerivative(first, second, face);
	double cpinormderivative = cpinormXDerivative(first, second, face);
	// calculating derivative 
	double derivativeValue = (1./pow(cpinorm,2))*(cpinorm*piyderivative - cpivector[1]*cpinormderivative);
	return derivativeValue;
}
// --------------------------------------------------------------------------------- //
double CentralisedDerivative::ex3XDerivative(Vertex* first, Vertex* second, Face* face){
	// gathering the terms //
	double *pntpivector = face->getPivector();
	double *pntnormal = face->getNormal();
	double cpivector[3] = {pntpivector[0]-face->getXCentralised(),pntpivector[1]-face->getYCentralised(),pntpivector[2]-face->getZCentralised()};
	double cpinorm = sqrt(pow(cpivector[0],2)+
							pow(cpivector[1],2) + 
							pow(cpivector[2],2));
	double pizderivative = pizXDerivative(first, second, face);
	double cpinormderivative = cpinormXDerivative(first, second, face);
	// calculating derivative 
	double derivativeValue = (1./pow(cpinorm,2))*(cpinorm*pizderivative - cpivector[2]*cpinormderivative);
	return derivativeValue;
}



// --------------------------------------------------------------------------------- //
double CentralisedDerivative::alphaXDerivative(Vertex* first, Vertex* second, Face* face){
	//gathering terms//
	// calculating derivative//
	return 0.0;
}
// --------------------------------------------------------------------------------- //
double CentralisedDerivative::betaXDerivative(Vertex* first, Vertex* second, Face* face){
	//gathering terms//
	//iterating for second vertex in the face
	double derivativeValue; //variable to store derivative
	unsigned int faceid = face->getID();
	Vertex * currentVertex;
	Edge * currentEdge;
	VertexEdgeIterator vertEdges(first);// iterating through all the edges coming out from vertex First
	unsigned int tempfaceid;
	Face * tempface;
	while((currentEdge = vertEdges.next())!=0){
		tempface = currentEdge->Left();
		tempfaceid = tempface->getID(); 
		if (tempfaceid == faceid){ // if the left face is same as Face* face, then we found the outgoing edge from First vert on the Face* face
			// currentEdge :  Org=firstvertex(j)----->-----Dest=vertex(j+1)
			Vertex* vertj1 = currentEdge->Dest();//the j+1 vertex
			double zj1 = vertj1->getZcoordinate();// z coordinate of j+1
			double zj = first->getZcoordinate(); // z coordinate of j
			derivativeValue = (zj-(face->getZCentralised()))*deltafunction(vertj1->getID(),second->getID())-
							  1./6.*(zj - zj1) - 
							  (zj1 - face->getZCentralised())*deltafunction(first->getID(),second->getID());
			return derivativeValue;
		}
	}
	// calculating derivative//
	return 0.0;
}
// --------------------------------------------------------------------------------- //
double CentralisedDerivative::gammaXDerivative(Vertex* first, Vertex* second, Face* face){
	//gathering terms//
	//iterating for second vertex in the face
	double derivativeValue; //variable to store derivative
	unsigned int faceid = face->getID();
	Vertex * currentVertex;
	Edge * currentEdge;
	VertexEdgeIterator vertEdges(first);// iterating through all the edges coming out from vertex First
	unsigned int tempfaceid;
	Face * tempface;
	while((currentEdge = vertEdges.next())!=0){
		tempface = currentEdge->Left();
		tempfaceid = tempface->getID(); 
		if (tempfaceid == faceid){ // if the left face is same as Face* face, then we found the outgoing edge from First vert on the Face* face
			// currentEdge :  Org=firstvertex(j)----->-----Dest=vertex(j+1)
			Vertex* vertj1 = currentEdge->Dest();//the j+1 vertex
			double yj1 = vertj1->getYcoordinate();// z coordinate of j+1
			double yj = first->getYcoordinate(); // z coordinate of j
			derivativeValue = (yj1-(face->getYCentralised()))*deltafunction(first->getID(),second->getID())-
							  1./6.*(yj1 - yj) - 
							  (yj - face->getYCentralised())*deltafunction(vertj1->getID(),second->getID());
			return derivativeValue;
		}
	}
	return 0.0;
}

// --------------------------------------------------------------------------------- //

double CentralisedDerivative::cpinormXDerivative(Vertex *first, Vertex *second, Face *face){
	// gathering terms
	double *pntpivector = face->getPivector();
	double *pntnormal = face->getNormal();
	double cPivector[3] = {pntpivector[0]-face->getXCentralised(),pntpivector[1]-face->getYCentralised(),pntpivector[2]-face->getZCentralised()};
	double cPiNorm = sqrt(pow(cPivector[0],2)+
							pow(cPivector[1],2) + 
							pow(cPivector[2],2));
	double pixderivative = pixXDerivative(first, second, face);
	double piyderivative = piyXDerivative(first, second, face);
	double pizderivative = pizXDerivative(first, second, face);
	// calculating the derivative value //
	double derivativeValue = (0.5/cPiNorm)*(2*(cPivector[0])*(pixderivative - 1./6.) +
											2*(cPivector[1])*(piyderivative) +
											2*(cPivector[2])*(pizderivative));
	return derivativeValue;
}


// --------------------------------------------------------------------------------- //
double CentralisedDerivative::ncxXDerivative(Vertex* first, Vertex* second, Face* face){
	// gathering terms //
	double *pntnormaltilde = face->getNormalTilde();
	double normalTildeNorm = sqrt(pow(pntnormaltilde[0],2)+pow(pntnormaltilde[1],2)+pow(pntnormaltilde[2],2));
	double ncxtildederivative = ncxtildeXDerivative(first, second, face);
	double nctildenormderivative = nctildenormXDerivative(first, second, face);
	// calculating derivative //
	double derivativeValue = 1./(pow(normalTildeNorm,2))*(normalTildeNorm*ncxtildederivative
														 - pntnormaltilde[0]*nctildenormderivative);
	return derivativeValue;	
}
// --------------------------------------------------------------------------------- //

double CentralisedDerivative::ncxtildeXDerivative(Vertex* first, Vertex* second, Face* face){
	// gathering terms // 
	double areatotal = face->getAreaOfFace();
	double areatotalderivative = areatotalXDerivative(first, second, face); 
	// iterating the faces to get the all the alpha X derivatives
	double alphaderivative(0);
	double alphasum(0);
	unsigned int faceid = face->getID();//ID of the face
	Vertex* currentVertex; // a keeper for current Vertex
	Edge* currentEdge; // a keeper for current Edge
	FaceEdgeIterator faceEdges(face);//iterator to iterate through the vertex for outgoign edge
    while((currentEdge = faceEdges.next())!= 0){//runnign through the edges in the face again
     		currentVertex = currentEdge->Dest();//vertex to be derivatived by second vertex
     		// summing up alpha X derivatives of all vertices
     		alphaderivative += alphaXDerivative(currentVertex, second, face);
     		//summing up all the alpha values of the vertex
     		alphasum += currentVertex->getAlpha(faceid);
     	}
     // ******** Evaluating the value of derivative *********** //
     double derivativeValue = 0.5*(1./pow(areatotal,2))*(areatotal*alphaderivative - 
     													alphasum*areatotalderivative);
     return derivativeValue;
}
// --------------------------------------------------------------------------------- //

double CentralisedDerivative::ncytildeXDerivative(Vertex* first, Vertex* second, Face* face){
	// gathering terms // 
	double areatotal = face->getAreaOfFace();
	double areatotalderivative = areatotalXDerivative(first, second, face); 
	// iterating the faces to get the all the alpha X derivatives
	double betaderivative(0);
	double betasum(0);
	unsigned int faceid = face->getID();//ID of the face
	Vertex* currentVertex; // a keeper for current Vertex
	Edge* currentEdge; // a keeper for current Edge
	FaceEdgeIterator faceEdges(face);//iterator to iterate through the vertex for outgoign edge
    while((currentEdge = faceEdges.next())!= 0){//runnign through the edges in the face again
     		currentVertex = currentEdge->Dest();//vertex to be derivatived by second vertex
     		// summing up alpha X derivatives of all vertices
     		betaderivative += betaXDerivative(currentVertex, second, face);
     		//summing up all the alpha values of the vertex
     		betasum += currentVertex->getBeta(faceid);
     	}
     // ******** Evaluating the value of derivative *********** //
     double derivativeValue = 0.5*(1./pow(areatotal,2))*(areatotal*betaderivative - 
     													betasum*areatotalderivative);
     return derivativeValue;
}

// --------------------------------------------------------------------------------- //

double CentralisedDerivative::ncztildeXDerivative(Vertex* first, Vertex* second, Face* face){
	// gathering terms // 
	double areatotal = face->getAreaOfFace();
	double areatotalderivative = areatotalXDerivative(first, second, face); 
	// iterating the faces to get the all the alpha X derivatives
	double gammaderivative(0);
	double gammasum(0);
	unsigned int faceid = face->getID();//ID of the face
	Vertex* currentVertex; // a keeper for current Vertex
	Edge* currentEdge; // a keeper for current Edge
	FaceEdgeIterator faceEdges(face);//iterator to iterate through the vertex for outgoign edge
    while((currentEdge = faceEdges.next())!= 0){//runnign through the edges in the face again
     		currentVertex = currentEdge->Dest();//vertex to be derivatived by second vertex
     		// summing up alpha X derivatives of all vertices
     		gammaderivative += gammaXDerivative(currentVertex, second, face);
     		//summing up all the alpha values of the vertex
     		gammasum += currentVertex->getGamma(faceid);
     	}
     // ******** Evaluating the value of derivative *********** //
     double derivativeValue = 0.5*(1./pow(areatotal,2))*(areatotal*gammaderivative - 
     													gammasum*areatotalderivative);
     return derivativeValue;
}


// --------------------------------------------------------------------------------- //

double CentralisedDerivative::areatotalXDerivative(Vertex* first, Vertex* second, Face* face){
	// gathering terms //
	double derivativeValue(0), alphatemp(0), betatemp(0), gammatemp(0), normfactor(0);
	double betaderivative(0), gammaderivative(0);
	unsigned int faceid = face->getID();
	FaceEdgeIterator faceEdges(face); // iterating the edges of this face
	Edge * currentEdge;
	Vertex* currentVertex;
	while ((currentEdge = faceEdges.next())!= 0){
			currentVertex = currentEdge->Dest(); //grabbing the destination vertex
			// gathering all terms //
			alphatemp = currentVertex->getAlpha(faceid);
			betatemp = currentVertex->getBeta(faceid);
			gammatemp = currentVertex->getGamma(faceid);
			normfactor = sqrt(pow(alphatemp,2)+pow(betatemp,2)+pow(gammatemp,2));
			betaderivative = betaXDerivative(currentVertex, second, face);
			gammaderivative = gammaXDerivative(currentVertex, second, face);
			//calculating derivative sum //
			derivativeValue += (1./normfactor)*(betatemp*betaderivative+gammatemp*gammaderivative);
	}
	derivativeValue = (1./(2))*derivativeValue;
	return derivativeValue;
}

// --------------------------------------------------------------------------------- //

double CentralisedDerivative::nctildenormXDerivative(Vertex* first, Vertex* second, Face* face){
	// gathering terms //
	double * normaltilde = face->getNormalTilde(); //getting the normal tilde vector of tha face
	double normaltildenorm = sqrt(pow(normaltilde[0],2)+pow(normaltilde[1],2)+pow(normaltilde[2],2));
	double ncxtildederivative = ncxtildeXDerivative(first, second, face);
	double ncytildederivative = ncytildeXDerivative(first, second, face);
	double ncztildederivative = ncztildeXDerivative(first, second, face);
	// calculating the derivative //
	double derivativeValue = (1./normaltildenorm) * (normaltilde[0]*ncxtildederivative +
													 normaltilde[1]*ncytildederivative +
													 normaltilde[2]*ncztildederivative);
	return derivativeValue;
}
// --------------------------------------------------------------------------------- //


double CentralisedDerivative::pixXDerivative(Vertex* first, Vertex* second, Face* face){
	// gathering terms //
	double *pntnormal = face->getNormal();//grabing the normal of the face
	double ncxderivative = ncxXDerivative(first, second, face);
	//calcuating the derivative //
	double derivativeValue = 1./6. - 2*pntnormal[0]*ncxderivative;
	return derivativeValue;
}
// --------------------------------------------------------------------------------- //


double CentralisedDerivative::piyXDerivative(Vertex* first, Vertex* second, Face* face){
	// gathering terms //
	double * normal = face->getNormal(); //grabbing the normal of the face
	double ncyderivative = ncyXDerivative(first, second, face);
	double ncxderivative = ncxXDerivative(first, second, face);
	// calculating the derivative // 
	double derivativeValue = - normal[0]*ncyderivative - normal[1]*ncxderivative;
	return derivativeValue;
}

// --------------------------------------------------------------------------------- //

double CentralisedDerivative::pizXDerivative(Vertex* first, Vertex* second, Face* face){
	// gathering terms //
	double * normal = face->getNormal(); //grabbing the normal of the face
	double nczderivative = nczXDerivative(first, second, face);
	double ncxderivative = ncxXDerivative(first, second, face);
	// calculating the derivative // 
	double derivativeValue = - normal[0]*nczderivative - normal[2]*ncxderivative;
	return derivativeValue;
}
// --------------------------------------------------------------------------------- //

double CentralisedDerivative::ncyXDerivative(Vertex* first, Vertex* second, Face* face){
	// gathering terms //
	double * normaltilde = face->getNormalTilde(); //getting the normal tilde vector of tha face
	double normaltildenorm = sqrt(pow(normaltilde[0],2)+pow(normaltilde[1],2)+pow(normaltilde[2],2)); // norm of normaltilde
	double ncytildederivative = ncytildeXDerivative(first, second, face);
	double normaltildenormderivative = nctildenormXDerivative(first, second, face);
	// calculating the derivative 
	double derivativeValue = (1./pow(normaltildenorm,2))*(normaltildenorm*ncytildederivative - normaltilde[1]*normaltildenormderivative);
	return derivativeValue;
}



// --------------------------------------------------------------------------------- //

double CentralisedDerivative::nczXDerivative(Vertex* first, Vertex* second, Face* face){
	// gathering terms //
	double * normaltilde = face->getNormalTilde(); //getting the normal tilde vector of tha face
	double normaltildenorm = sqrt(pow(normaltilde[0],2)+pow(normaltilde[1],2)+pow(normaltilde[2],2)); // norm of normaltilde
	double ncztildederivative = ncztildeXDerivative(first, second, face);
	double normaltildenormderivative = nctildenormXDerivative(first, second, face);
	// calculating the derivative 
	double derivativeValue = (1./pow(normaltildenorm,2))*(normaltildenorm*ncztildederivative - normaltilde[2]*normaltildenormderivative);
	return derivativeValue;
}
// --------------------------------------------------------------------------------- //
double CentralisedDerivative::ey1XDerivative(Vertex* first, Vertex* second, Face* face){
	//gathering terms // 
	double * unitxvector = face->getUnitx();
	double * unitzvector = face->getUnitz();
	double ncyderivative = ncyXDerivative(first, second, face);
	double nczderivative = nczXDerivative(first, second, face);
	double ex2derivative = ex2XDerivative(first, second, face);
	double ex3derivative = ex3XDerivative(first, second, face);
	// calculating terms
	double derivativevalue = unitxvector[1]*nczderivative+unitzvector[2]*ex2derivative-
							 unitxvector[2]*ncyderivative-unitzvector[1]*ex3derivative;
	//returning 
	return derivativevalue;
}
// --------------------------------------------------------------------------------- //
double CentralisedDerivative::ey2XDerivative(Vertex* first, Vertex* second, Face* face){
	//gathering terms // 
	double * unitxvector = face->getUnitx();
	double * unitzvector = face->getUnitz();
	double ncxderivative = ncxXDerivative(first, second, face);
	double nczderivative = nczXDerivative(first, second, face);
	double ex1derivative = ex1XDerivative(first, second, face);
	double ex3derivative = ex3XDerivative(first, second, face);
	// calculating terms
	double derivativevalue = unitxvector[2]*ncxderivative+unitzvector[0]*ex3derivative-
							 unitxvector[0]*nczderivative-unitzvector[2]*ex1derivative;
	//returning 
	return derivativevalue;
}
// --------------------------------------------------------------------------------- //
double CentralisedDerivative::ey3XDerivative(Vertex* first, Vertex* second, Face* face){
	//gathering terms // 
	double * unitxvector = face->getUnitx();
	double * unitzvector = face->getUnitz();
	double ncxderivative = ncxXDerivative(first, second, face);
	double ncyderivative = ncyXDerivative(first, second, face);
	double ex1derivative = ex1XDerivative(first, second, face);
	double ex2derivative = ex2XDerivative(first, second, face);
	// calculating terms
	double derivativevalue = unitxvector[0]*ncyderivative+unitzvector[1]*ex1derivative-
							 unitxvector[1]*ncxderivative-unitzvector[0]*ex2derivative;
	//returning 
	return derivativevalue;
}
/*
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
double CentralisedDerivative::numericalXtildeXDerivative(Vertex* first, Vertex* second, Face* face, double stepsizeH = pow(10.,-8)){
  //setting the projected coordinates of all vertices in this face
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
    // getting the mean center position of this face //
    unsigned int faceid = face->getID();// first grabbingat id of  current face
   	unsigned int vertexid1 = first->getID();//geting id of first vertex
   	unsigned int vertexid2 = second->getID();//getting id of second vertex
   	// step size for derivative 
   	//double epsilon = std::numeric_limits<double>::epsilon();// smallest precision for calculation, ~10^-16
   	//double stepsizeH = sqrt(epsilon);// ~10^-8
   	// excluding face 1 because calcualtions are not need for it
    if (faceid == 1){//if faceid == 1 then face is external face, that doesnot need projection or calculation
        return 0.0;
    } 
    //starting the calculation for xtilde for vertex 1 with x_k+h of vertex2 
    FaceEdgeIterator faceEdges(face);//iterator to iterate through the vertex for outgoign edge
    double nx,ny;//normal coordinate
    double xCentroid(0), yCentroid(0), zCentroid(0); // coordinate of the cnetroid
    // array of vertices
    double xcood[6], ycood[6], zcood[6];// coordinate of the vertices of this face
    double xTriCen[6], yTriCen[6], zTriCen[6], areaTri;// coordinate of center of triangluated triangles of face
    //areaTri : area of triangles
    double xmean(0), ymean(0), zmean(0);//mean center of the face
    int counter(0);//counter to keep track of vertices
    Edge *currentEdge;//pointer to keep track of a edge that is iterated 
    Vertex *currentVertex; // vertex to store the dest() of currentEdge
    while ((currentEdge = faceEdges.next())!=0){//iterating the edges in the face
          currentVertex = currentEdge->Dest();
          if (currentVertex->getID() ==vertexid2){
          		//changing the value of x coordinate of second vertex by 'h' or stepsize
          		xcood[counter] = currentVertex->getXcoordinate() + stepsizeH;
          		ycood[counter] = currentVertex->getYcoordinate();
          		zcood[counter] = currentVertex->getZcoordinate();
          		xmean += xcood[counter];
          		ymean += ycood[counter];
          		zmean += zcood[counter];
          		counter += 1;//increasing the counter value
          }
          else
          {
         		// getting the vertices of the current face and storring it in array
          		xcood[counter] = currentVertex->getXcoordinate();
		        ycood[counter] = currentVertex->getYcoordinate();
		        zcood[counter] = currentVertex->getZcoordinate();
		        xmean += xcood[counter];
		        ymean += ycood[counter];
		        zmean += zcood[counter];
	            counter += 1;//increasing the counter value
      	  }
    }
    //printf(" face id : %u \n", faceid);
    //printf(" calculated means : X = %F ; Y = %F ; Z = %F \n", xmean, ymean,zmean );
    //printf("counter %h \n",counter);
    // divinding by the number of vertices to get the mean
    xmean = (1./counter)*xmean;
    ymean = (1./counter)*ymean;
    zmean = (1./counter)*zmean;
    //printf(">>>>new output start<<<< face id = %u \n", faceid);
    //printf("calculated means : X = %F ; Y = %F ; Z = %F \n", xmean, ymean,zmean );
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
    // triangulating the face and getting the center of each triangle
    //vectors of the triangles 
    //now calculating weighted center
    xCentroid = xmean;
    yCentroid = ymean;
    zCentroid = zmean;
    //printf("%s %u %F \n ", "calculating total area of Face :",faceid, totalarea);
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
    //Calculating normal of the triangles of the shape and then calculating weighted normal
    double normalTempX[6],normalTempY[6],normalTempZ[6],totalarea;//vertices of normal of triangles
    double vector1[3], vector2[3], crossProductVector[3],crossProductMagnitude;
    totalarea = 0;
    double normalX(0), normalY(0), normalZ(0); // vertices of weighted normal X, Y, Z 
    // iterating the faces also
    FaceEdgeIterator faceEdges1(face);//iterator to iterate through the vertex for outgoign edge
    // curentVertex and currentEdge have been already declared above
    currentEdge = faceEdges1.next();// grabbing the one edge of face = this
    for (int counter = 0; counter<6; counter++){//triangulating the vertices with centroid of this face
          //currentVertex = currentEdge->Dest();//getting the destination of the edge
          //getting two vectors of this triangle
          vector1[0] = xcood[counter]-xCentroid;
          vector1[1] = ycood[counter]-yCentroid;
          vector1[2] = zcood[counter]-zCentroid;
          vector2[0] = xcood[(counter+1)%6]-xCentroid;
          vector2[1] = ycood[(counter+1)%6]-yCentroid;             
          vector2[2] = zcood[(counter+1)%6]-zCentroid;
          //cross product of the two vectors
          crossProductVector[0] = vector1[1]*vector2[2] - vector1[2]*vector2[1];//also definition of alpha
          crossProductVector[1] = vector1[2]*vector2[0] - vector1[0]*vector2[2];//also definition of beta
          crossProductVector[2] = vector1[0]*vector2[1] - vector1[1]*vector2[0];//also definition of gamma
          //maginitude of cross product
          crossProductMagnitude = sqrt(crossProductVector[0]*crossProductVector[0] + crossProductVector[1]*crossProductVector[1]+
                                  crossProductVector[2]*crossProductVector[2]);
          //normalising the crossProductVector
          crossProductVector[0] = crossProductVector[0]/crossProductMagnitude;
          crossProductVector[1] = crossProductVector[1]/crossProductMagnitude;
          crossProductVector[2] = crossProductVector[2]/crossProductMagnitude;
          areaTri = 0.5*crossProductMagnitude;//area of this triangle
          //adding the weigthed centroid 
          normalX += crossProductVector[0]*areaTri;
          normalY += crossProductVector[1]*areaTri;
          normalZ += crossProductVector[2]*areaTri; 
          totalarea += areaTri;
          //total area is same as before calculated for this face
          currentEdge = faceEdges1.next();//keep on iterating through the edges in the face
      }
    //printf("%s %u %F \n ", "calculating total area of Face :",faceid, totalarea);
    //weighted normal to this face
    normalX = normalX/totalarea;
    normalY = normalY/totalarea;
    normalZ = normalZ/totalarea;
    //normalising the normal vector 
    double normalMagnitude = sqrt(pow(normalX,2.0)+pow(normalY,2.0)+pow(normalZ,2.0));
    normalX = normalX/normalMagnitude;
    normalY = normalY/normalMagnitude;
    normalZ = normalZ/normalMagnitude;
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
    //now need to get this projected 3d plane coordinates onto rotated 2d coordinate system
    // getting the unit vector of 2d Plane
    double unitx[3], unity[3], vectorVertex[3];//new unit vector of x, y on the plane
    double dotproduct, xvertex, yvertex, zvertex;
    //to get the unit vector in X direction, lets project X axis (1,0,0) on to the plane
    //dot product of vectorVertex and normal, which is equal to normalX; this is n^c_x = \phi_P
    dotproduct = normalX;//vectorVertex[0]*normalX+vectorVertex[1]*normalY+vectorVertex[2]*normalZ;
    //now calculating the projected vertices --at this stage unitx is storing pi_vector
    unitx[0] = 1.0 - dotproduct*normalX;
    unitx[1] = 0.0 - dotproduct*normalY;
    unitx[2] = 0.0 - dotproduct*normalZ;
    // normalising unitx
    double normUnitx   = sqrt(pow(unitx[0],2)+pow(unitx[1],2)+pow(unitx[2],2));
    unitx[0] = unitx[0]/normUnitx;
    unitx[1] = unitx[1]/normUnitx;
    unitx[2] = unitx[2]/normUnitx;
    // now getting unity : it is the cross product of  normal and unitx to the plane. 
    unity[0] = normalY*unitx[2]-normalZ*unitx[1];
    unity[1] = normalZ*unitx[0]-normalX*unitx[2];
    unity[2] = normalX*unitx[1]-normalY*unitx[0];
    //normalising unity
    double normUnity = sqrt(pow(unity[0],2)+pow(unity[1],2)+pow(unity[2],2));
    unity[0] = unity[0]/normUnity;
    unity[1] = unity[1]/normUnity;
    unity[2] = unity[2]/normUnity;
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
    //%%%%%%%%%%%%%%%%%%%Assigning new coordinates%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
    // now projecting the 3d coordinates to the 2D coordinates and assinging the vertices
    double xprojection2; // this is f(x+h)
    FaceEdgeIterator faceEdges2(face);//iterator to iterate through the vertex for outgoign edge
    while((currentEdge = faceEdges2.next())!= 0){//runnign through the edges in the face again
          currentVertex = currentEdge->Dest();// the current vertex, iterated
          if (currentVertex->getID() == vertexid1){
          xvertex = currentVertex->getXcoordinate();
          yvertex = currentVertex->getYcoordinate();
          zvertex = currentVertex->getZcoordinate();
          //getting the vector form the centroid or new origin
          // CENTRALISED PorjectedCoordinate
          vectorVertex[0] = xvertex-xCentroid;
          vectorVertex[1] = yvertex-yCentroid;
          vectorVertex[2] = zvertex-zCentroid;
          //  !!! CHANGE HERE !!!
          // NON -- CENTRALISED PorjectedCoordinate
          vectorVertex[0] = xvertex;
          vectorVertex[1] = yvertex;
          vectorVertex[2] = zvertex;
          //now getting the new x, y coordinates, dot product of unit vector with the vectorVertex
          xprojection2 = unitx[0]*vectorVertex[0]+unitx[1]*vectorVertex[1]+unitx[2]*vectorVertex[2];
          //yprojection2 = unity[0]*vectorVertex[0]+unity[1]*vectorVertex[1]+unity[2]*vectorVertex[2];
          //zprojection2 = normalX*vectorVertex[0]+normalY*vectorVertex[1]+normalZ*vectorVertex[2];
          }
    }
    double xprojection1; // this is initial value or f(x)
    xprojection1 = first->getProjectedXcoordinate(faceid);
    // *********************************************************************************** //
    //calculating derivative value, f' =(f(x+h)-f(x))/h
   double derivativeValue =  (xprojection2 - xprojection1)/stepsizeH;
   return derivativeValue;
}
*/


// ================================================================================================= //
// ========================= Y DERIVATIVE ========================================================== //
// ================================================================================================= //

double CentralisedDerivative::xtildeYDerivative(Vertex *first, Vertex *second, Face *face){
	// gathering terms //
	if (face->getID() == 1){
		return 0.0;
	}
	else{
	double* unitxvector = face->getUnitx();
	double coordinates[3] = {first->getXcoordinate(),
							first->getYcoordinate(),
							first->getZcoordinate()};
	double centroid[3] ={face->getXCentralised(),
						face->getYCentralised(),
						face->getZCentralised()};
	//double * unitx = face->getUnitx();
	// calculating derivatives //
	double yiderivative = yiYDerivative(first, second, face);
	double ex1derivative =  ex1YDerivative(first, second, face);
	double ex2derivative = ex2YDerivative(first, second, face);
	double ex3derivative =  ex3YDerivative(first, second, face);
	// final calcuation //
	double derivativeValue = (coordinates[0]-centroid[0])*ex1derivative
							 +(coordinates[1]-centroid[1])*ex2derivative
							 +(coordinates[2]-centroid[2])*ex3derivative
							 +unitxvector[1]*(yiderivative-1./6.);
	// returning value //
	return derivativeValue;
	}
	
}
// -------------------------------------------------------------------------------------------------- //

double CentralisedDerivative::ytildeYDerivative(Vertex *first, Vertex *second, Face *face){
	// gathering terms //
	if (face->getID() == 1){
		return 0.0;
	}
	else{
	//double* unitxvector = face->getUnitx();
	double coordinates[3] = {first->getXcoordinate(),
							first->getYcoordinate(),
							first->getZcoordinate()};
	double centroid[3] ={face->getXCentralised(),
						face->getYCentralised(),
						face->getZCentralised()};
	double * unityvector = face->getUnity();
	// calculating derivatives //
	double yiderivative = yiYDerivative(first, second, face);
	double ey1derivative = ey1YDerivative(first, second, face);
	double ey2derivative = ey2YDerivative(first, second, face);
	double ey3derivative = ey3YDerivative(first, second, face);
	// final calcuation //
	double derivativeValue = (coordinates[0]-centroid[0])*ey1derivative
							 +(coordinates[1]-centroid[1])*ey2derivative
							 +(coordinates[2]-centroid[2])*ey3derivative
							 +unityvector[1]*(yiderivative-1./6.);
	// returning value //
	//std::cout<<"Unit Y : "<<unityvector[0]<<unityvector[1]<<unityvector[3]<<std::endl;
	return derivativeValue;
	}
	
}

// -------------------------------------------------------------------------------------------------- //

double CentralisedDerivative::alphaYDerivative(Vertex* first, Vertex* second, Face* face){
	//gathering terms//
	//iterating for second vertex in the face
	double derivativeValue; //variable to store derivative
	unsigned int faceid = face->getID();
	Vertex * currentVertex;
	Edge * currentEdge;
	VertexEdgeIterator vertEdges(first);// iterating through all the edges coming out from vertex First
	unsigned int tempfaceid;
	Face * tempface;
	while((currentEdge = vertEdges.next())!=0){
		tempface = currentEdge->Left();
		tempfaceid = tempface->getID(); 
		if (tempfaceid == faceid){ // if the left face is same as Face* face, then we found the outgoing edge from First vert on the Face* face
			// currentEdge :  Org=firstvertex(j)----->-----Dest=vertex(j+1)
			Vertex* vertj1 = currentEdge->Dest();//the j+1 vertex
			double zj1 = vertj1->getZcoordinate();// z coordinate of j+1
			double zj = first->getZcoordinate(); // z coordinate of j
			derivativeValue = (zj1-(face->getZCentralised()))*deltafunction(first->getID(),second->getID())
							  + 1./6.*(zj - zj1)
							  - (zj - face->getZCentralised())*deltafunction(vertj1->getID(),second->getID());
			return derivativeValue;
		}
	}
	// if something is gone wrong
	return 0.;
}
// -------------------------------------------------------------------------------------------------- //

double CentralisedDerivative::betaYDerivative(Vertex* first, Vertex* second, Face* face){
	return 0.0;
}

// -------------------------------------------------------------------------------------------------- //

double CentralisedDerivative::gammaYDerivative(Vertex* first, Vertex* second, Face* face){
	//gathering terms//
	//iterating for second vertex in the face
	double derivativeValue; //variable to store derivative
	unsigned int faceid = face->getID();
	Vertex * currentVertex;
	Edge * currentEdge;
	VertexEdgeIterator vertEdges(first);// iterating through all the edges coming out from vertex First
	unsigned int tempfaceid;
	Face * tempface;
	while((currentEdge = vertEdges.next())!=0){
		tempface = currentEdge->Left();
		tempfaceid = tempface->getID(); 
		if (tempfaceid == faceid){ // if the left face is same as Face* face, then we found the outgoing edge from First vert on the Face* face
			// currentEdge :  Org=firstvertex(j)----->-----Dest=vertex(j+1)
			Vertex* vertj1 = currentEdge->Dest();//the j+1 vertex
			double xj1 = vertj1->getXcoordinate();// z coordinate of j+1
			double xj = first->getXcoordinate(); // z coordinate of j
			derivativeValue = (xj-(face->getXCentralised()))*deltafunction(vertj1->getID(),second->getID())
							  - 1./6.*(xj - xj1)
							  - (xj1 - face->getXCentralised())*deltafunction(first->getID(),second->getID());
			return derivativeValue;
		}
	}
	// if something is gone wrong
	return 0.0;
}
// -------------------------------------------------------------------------------------------------- //
double CentralisedDerivative::ex1YDerivative(Vertex *first, Vertex *second, Face *face){
	//gathering terms // 
	double *pntpivector = face->getPivector();
	double *pntnormal = face->getNormal();
	double cPiNorm = sqrt(pow(pntpivector[0]-face->getXCentralised(),2)+
							pow(pntpivector[1]-face->getYCentralised(),2) + 
							pow(pntpivector[2]-face->getZCentralised(),2));
	double ncyderivative = ncxYDerivative(first, second, face); //calculating the ncxderivative
	double cPiNormyderivative = cpinormYDerivative(first,second,face); //calculating the cpinorm X derivative
	//calculating the derivative value // 
	double derivativeValue = (cPiNorm*(-2.*pntnormal[0]*ncyderivative)-(1-pow(pntnormal[0],2))*cPiNormyderivative)/pow(cPiNorm,2.);
	return derivativeValue;
}
// -------------------------------------------------------------------------------------------------- //
double CentralisedDerivative::ex2YDerivative(Vertex* first, Vertex* second, Face* face){
	// gathering the terms //
	double *pntpivector = face->getPivector();
	double *pntnormal = face->getNormal();
	double cpivector[3] = {pntpivector[0]-face->getXCentralised(),pntpivector[1]-face->getYCentralised(),pntpivector[2]-face->getZCentralised()};
	double cpinorm = sqrt(pow(cpivector[0],2)+
							pow(cpivector[1],2) + 
							pow(cpivector[2],2));
	double piyderivative = piyYDerivative(first, second, face);
	double cpinormderivative = cpinormYDerivative(first, second, face);
	// calculating derivative 
	double derivativeValue = (1./pow(cpinorm,2))*(cpinorm*(piyderivative-1./6.) - (cpivector[1])*cpinormderivative);
	return derivativeValue;
}
// --------------------------------------------------------------------------------- //
double CentralisedDerivative::ex3YDerivative(Vertex* first, Vertex* second, Face* face){
	// gathering the terms //
	double *pntpivector = face->getPivector();
	double *pntnormal = face->getNormal();
	double cpivector[3] = {pntpivector[0]-face->getXCentralised(),pntpivector[1]-face->getYCentralised(),pntpivector[2]-face->getZCentralised()};
	double cpinorm = sqrt(pow(cpivector[0],2)+
							pow(cpivector[1],2) + 
							pow(cpivector[2],2));
	double pizderivative = pizYDerivative(first, second, face);
	double cpinormderivative = cpinormYDerivative(first, second, face);
	// calculating derivative 
	double derivativeValue = (1./pow(cpinorm,2))*(cpinorm*pizderivative - cpivector[2]*cpinormderivative);
	return derivativeValue;
}
// --------------------------------------------------------------------------------- //

double CentralisedDerivative::ey1YDerivative(Vertex* first, Vertex* second, Face* face){
	//gathering terms // 
	double * unitxvector = face->getUnitx();
	double * unitzvector = face->getUnitz();
	double ncyderivative = ncyYDerivative(first, second, face);
	double nczderivative = nczYDerivative(first, second, face);
	double ex2derivative = ex2YDerivative(first, second, face);
	double ex3derivative = ex3YDerivative(first, second, face);
	// calculating terms
	double derivativevalue = unitxvector[1]*nczderivative+unitzvector[2]*ex2derivative-
							 unitxvector[2]*ncyderivative-unitzvector[1]*ex3derivative;
	//returning 
	return derivativevalue;
}
// --------------------------------------------------------------------------------- //
double CentralisedDerivative::ey2YDerivative(Vertex* first, Vertex* second, Face* face){
	//gathering terms // 
	double * unitxvector = face->getUnitx();
	double * unitzvector = face->getUnitz();
	double ncxderivative = ncxYDerivative(first, second, face);
	double nczderivative = nczYDerivative(first, second, face);
	double ex1derivative = ex1YDerivative(first, second, face);
	double ex3derivative = ex3YDerivative(first, second, face);
	// calculating terms
	double derivativevalue = unitxvector[2]*ncxderivative+unitzvector[0]*ex3derivative-
							 unitxvector[0]*nczderivative-unitzvector[2]*ex1derivative;
	//returning 
	return derivativevalue;
}
// --------------------------------------------------------------------------------- //
double CentralisedDerivative::ey3YDerivative(Vertex* first, Vertex* second, Face* face){
	//gathering terms // 
	double * unitxvector = face->getUnitx();
	double * unitzvector = face->getUnitz();
	double ncxderivative = ncxYDerivative(first, second, face);
	double ncyderivative = ncyYDerivative(first, second, face);
	double ex1derivative = ex1YDerivative(first, second, face);
	double ex2derivative = ex2YDerivative(first, second, face);
	// calculating terms
	double derivativevalue = unitxvector[0]*ncyderivative+unitzvector[1]*ex1derivative-
							 unitxvector[1]*ncxderivative-unitzvector[0]*ex2derivative;
	//returning 
	return derivativevalue;
}
// --------------------------------------------------------------------------------- //

double CentralisedDerivative::areatotalYDerivative(Vertex* first, Vertex* second, Face* face){
	// gathering terms //
	double derivativeValue(0), alphatemp(0), betatemp(0), gammatemp(0), normfactor(0);
	double alphaderivative(0),betaderivative(0), gammaderivative(0);
	unsigned int faceid = face->getID();
	FaceEdgeIterator faceEdges(face); // iterating the edges of this face
	Edge * currentEdge;
	Vertex* currentVertex;
	while ((currentEdge = faceEdges.next())!= 0){
			currentVertex = currentEdge->Dest(); //grabbing the destination vertex
			// gathering all terms //
			alphatemp = currentVertex->getAlpha(faceid);
			betatemp = currentVertex->getBeta(faceid);
			gammatemp = currentVertex->getGamma(faceid);
			normfactor = sqrt(pow(alphatemp,2)+pow(betatemp,2)+pow(gammatemp,2));
			alphaderivative = alphaYDerivative(currentVertex, second, face);
			betaderivative = betaYDerivative(currentVertex, second, face);
			gammaderivative = gammaYDerivative(currentVertex, second, face);
			//calculating derivative sum //
			derivativeValue += (1./normfactor)*(alphatemp*alphaderivative+betatemp*betaderivative+gammatemp*gammaderivative);
	}
	derivativeValue = (1./(2))*derivativeValue;
	return derivativeValue;
}

// -------------------------------------------------------------------------------------------------- //

double CentralisedDerivative::yiYDerivative(Vertex* first, Vertex* second, Face* face){
	return deltafunction(first->getID(),second->getID());
}
// -------------------------------------------------------------------------------------------------- //
double CentralisedDerivative::cpinormYDerivative(Vertex *first, Vertex *second, Face *face){
	// gathering terms
	double *pntpivector = face->getPivector();
	double *pntnormal = face->getNormal();
	double cPivector[3] = {pntpivector[0]-face->getXCentralised(),pntpivector[1]-face->getYCentralised(),pntpivector[2]-face->getZCentralised()};
	double cPiNorm = sqrt(pow(cPivector[0],2)+
							pow(cPivector[1],2) + 
							pow(cPivector[2],2));
	double pixderivative = pixYDerivative(first, second, face);
	double piyderivative = piyYDerivative(first, second, face);
	double pizderivative = pizYDerivative(first, second, face);
	// calculating the derivative value //
	double derivativeValue = (0.5/cPiNorm)*(2*(cPivector[0])*(pixderivative) +
											2*(cPivector[1])*(piyderivative - 1./6.) +
											2*(cPivector[2])*(pizderivative));
	return derivativeValue;
}

// -------------------------------------------------------------------------------------------------- //
double CentralisedDerivative::ncxYDerivative(Vertex* first, Vertex* second, Face* face){
	// gathering terms //
	double *pntnormaltilde = face->getNormalTilde();
	double normalTildeNorm = sqrt(pow(pntnormaltilde[0],2)+pow(pntnormaltilde[1],2)+pow(pntnormaltilde[2],2));
	double ncxtildederivative = ncxtildeYDerivative(first, second, face);
	double nctildenormderivative = nctildenormYDerivative(first, second, face);
	// calculating derivative //
	double derivativeValue = 1./(pow(normalTildeNorm,2))*(normalTildeNorm*ncxtildederivative
														 - pntnormaltilde[0]*nctildenormderivative);
	return derivativeValue;	
}
// -------------------------------------------------------------------------------------------------- //

double CentralisedDerivative::ncyYDerivative(Vertex* first, Vertex* second, Face* face){
	// gathering terms //
	double * normaltilde = face->getNormalTilde(); //getting the normal tilde vector of tha face
	double normaltildenorm = sqrt(pow(normaltilde[0],2)+pow(normaltilde[1],2)+pow(normaltilde[2],2)); // norm of normaltilde
	double ncytildederivative = ncytildeYDerivative(first, second, face);
	double normaltildenormderivative = nctildenormYDerivative(first, second, face);
	// calculating the derivative 
	double derivativeValue = (1./pow(normaltildenorm,2))*(normaltildenorm*ncytildederivative - normaltilde[1]*normaltildenormderivative);
	return derivativeValue;
}

// --------------------------------------------------------------------------------- //

double CentralisedDerivative::nczYDerivative(Vertex* first, Vertex* second, Face* face){
	// gathering terms //
	double * normaltilde = face->getNormalTilde(); //getting the normal tilde vector of tha face
	double normaltildenorm = sqrt(pow(normaltilde[0],2)+pow(normaltilde[1],2)+pow(normaltilde[2],2)); // norm of normaltilde
	double ncztildederivative = ncztildeYDerivative(first, second, face);
	double normaltildenormderivative = nctildenormYDerivative(first, second, face);
	// calculating the derivative 
	double derivativeValue = (1./pow(normaltildenorm,2))*(normaltildenorm*ncztildederivative - normaltilde[2]*normaltildenormderivative);
	return derivativeValue;
}
// -------------------------------------------------------------------------------------------------- //
double CentralisedDerivative::nctildenormYDerivative(Vertex* first, Vertex* second, Face* face){
		// gathering terms //
	double * normaltilde = face->getNormalTilde(); //getting the normal tilde vector of tha face
	double normaltildenorm = sqrt(pow(normaltilde[0],2)+pow(normaltilde[1],2)+pow(normaltilde[2],2));
	double ncxtildederivative = ncxtildeYDerivative(first, second, face);
	double ncytildederivative = ncytildeYDerivative(first, second, face);
	double ncztildederivative = ncztildeYDerivative(first, second, face);
	// calculating the derivative //
	double derivativeValue = (1./normaltildenorm) * (normaltilde[0]*ncxtildederivative +
													 normaltilde[1]*ncytildederivative +
													 normaltilde[2]*ncztildederivative);
	return derivativeValue;
}
// -------------------------------------------------------------------------------------------------- //
double CentralisedDerivative::ncxtildeYDerivative(Vertex* first, Vertex* second, Face* face){
	// gathering terms // 
	double areatotal = face->getAreaOfFace();
	double areatotalderivative = areatotalYDerivative(first, second, face); 
	// iterating the faces to get the all the alpha X derivatives
	double alphaderivative(0);
	double alphasum(0);
	unsigned int faceid = face->getID();//ID of the face
	Vertex* currentVertex; // a keeper for current Vertex
	Edge* currentEdge; // a keeper for current Edge
	FaceEdgeIterator faceEdges(face);//iterator to iterate through the vertex for outgoign edge
    while((currentEdge = faceEdges.next())!= 0){//runnign through the edges in the face again
     		currentVertex = currentEdge->Dest();//vertex to be derivatived by second vertex
     		// summing up alpha X derivatives of all vertices
     		alphaderivative += alphaYDerivative(currentVertex, second, face);
     		//summing up all the alpha values of the vertex
     		alphasum += currentVertex->getAlpha(faceid);
     	}
     // ******** Evaluating the value of derivative *********** //
     double derivativeValue = 0.5*(1./pow(areatotal,2))*(areatotal*alphaderivative - 
     													alphasum*areatotalderivative);
     return derivativeValue;
}

// -------------------------------------------------------------------------------------------------- //
double CentralisedDerivative::ncytildeYDerivative(Vertex* first, Vertex* second, Face* face){
	// gathering terms // 
	double areatotal = face->getAreaOfFace();
	double areatotalderivative = areatotalYDerivative(first, second, face); 
	// iterating the faces to get the all the alpha X derivatives
	double betaderivative(0);
	double betasum(0);
	unsigned int faceid = face->getID();//ID of the face
	Vertex* currentVertex; // a keeper for current Vertex
	Edge* currentEdge; // a keeper for current Edge
	FaceEdgeIterator faceEdges(face);//iterator to iterate through the vertex for outgoign edge
    while((currentEdge = faceEdges.next())!= 0){//runnign through the edges in the face again
     		currentVertex = currentEdge->Dest();//vertex to be derivatived by second vertex
     		// summing up alpha X derivatives of all vertices
     		betaderivative += betaYDerivative(currentVertex, second, face);
     		//summing up all the alpha values of the vertex
     		betasum += currentVertex->getBeta(faceid);
     	}
     // ******** Evaluating the value of derivative *********** //
     double derivativeValue = 0.5*(1./pow(areatotal,2))*(areatotal*betaderivative - 
     													betasum*areatotalderivative);
     return derivativeValue;
}


// --------------------------------------------------------------------------------- //
double CentralisedDerivative::ncztildeYDerivative(Vertex* first, Vertex* second, Face* face){
	// gathering terms // 
	double areatotal = face->getAreaOfFace();
	double areatotalderivative = areatotalYDerivative(first, second, face); 
	// iterating the faces to get the all the alpha X derivatives
	double gammaderivative(0);
	double gammasum(0);
	unsigned int faceid = face->getID();//ID of the face
	Vertex* currentVertex; // a keeper for current Vertex
	Edge* currentEdge; // a keeper for current Edge
	FaceEdgeIterator faceEdges(face);//iterator to iterate through the vertex for outgoign edge
    while((currentEdge = faceEdges.next())!= 0){//runnign through the edges in the face again
     		currentVertex = currentEdge->Dest();//vertex to be derivatived by second vertex
     		// summing up alpha X derivatives of all vertices
     		gammaderivative += gammaYDerivative(currentVertex, second, face);
     		//summing up all the alpha values of the vertex
     		gammasum += currentVertex->getGamma(faceid);
     	}
     // ******** Evaluating the value of derivative *********** //
     double derivativeValue = 0.5*(1./pow(areatotal,2))*(areatotal*gammaderivative - 
     													gammasum*areatotalderivative);
     return derivativeValue;
}



// --------------------------------------------------------------------------------- //

double CentralisedDerivative::pixYDerivative(Vertex* first, Vertex* second, Face* face){
	// gathering terms //
	double *pntnormal = face->getNormal();//grabing the normal of the face
	double ncxderivative = ncxYDerivative(first, second, face);
	//calcuating the derivative //
	double derivativeValue = -2.*pntnormal[0]*ncxderivative;
	return derivativeValue;
}
// --------------------------------------------------------------------------------- //


double CentralisedDerivative::piyYDerivative(Vertex* first, Vertex* second, Face* face){
	// gathering terms //
	double * normal = face->getNormal(); //grabbing the normal of the face
	double ncyderivative = ncyYDerivative(first, second, face);
	double ncxderivative = ncxYDerivative(first, second, face);
	// calculating the derivative // 
	double derivativeValue = 1./6.- normal[0]*ncyderivative - normal[1]*ncxderivative;
	return derivativeValue;
}

// --------------------------------------------------------------------------------- //

double CentralisedDerivative::pizYDerivative(Vertex* first, Vertex* second, Face* face){
	// gathering terms //
	double * normal = face->getNormal(); //grabbing the normal of the face
	double nczderivative = nczYDerivative(first, second, face);
	double ncxderivative = ncxYDerivative(first, second, face);
	// calculating the derivative // 
	double derivativeValue = - normal[0]*nczderivative - normal[2]*ncxderivative;
	return derivativeValue;
}


// ================================================================================================= //
// ========================= Z DERIVATIVE ========================================================== //
// ================================================================================================= //

double CentralisedDerivative::xtildeZDerivative(Vertex *first, Vertex *second, Face *face){
	// gathering terms //
	if (face->getID() == 1){
		return 0.0;
	}
	else{
	double* unitxvector = face->getUnitx();
	double coordinates[3] = {first->getXcoordinate(),
							first->getYcoordinate(),
							first->getZcoordinate()};
	double centroid[3] ={face->getXCentralised(),
						face->getYCentralised(),
						face->getZCentralised()};
	//double * unitx = face->getUnitx();
	// calculating derivatives //
	double ziderivative = ziZDerivative(first, second, face);
	double ex1derivative =  ex1ZDerivative(first, second, face);
	double ex2derivative = ex2ZDerivative(first, second, face);
	double ex3derivative =  ex3ZDerivative(first, second, face);
	// final calcuation //
	double derivativeValue = (coordinates[0]-centroid[0])*ex1derivative
							 +(coordinates[1]-centroid[1])*ex2derivative
							 +(coordinates[2]-centroid[2])*ex3derivative
							 +unitxvector[2]*(ziderivative-1./6.);
	// returning value //
	return derivativeValue;
	}
	
}
// -------------------------------------------------------------------------------------------------- //

double CentralisedDerivative::ytildeZDerivative(Vertex *first, Vertex *second, Face *face){
	// gathering terms //
	if (face->getID() == 1){
		return 0.0;
	}
	else{
	//double* unitxvector = face->getUnitx();
	double coordinates[3] = {first->getXcoordinate(),
							first->getYcoordinate(),
							first->getZcoordinate()};
	double centroid[3] ={face->getXCentralised(),
						face->getYCentralised(),
						face->getZCentralised()};
	double * unityvector = face->getUnity();
	// calculating derivatives //
	double ziderivative = ziZDerivative(first, second, face);
	double ey1derivative = ey1ZDerivative(first, second, face);
	double ey2derivative = ey2ZDerivative(first, second, face);
	double ey3derivative = ey3ZDerivative(first, second, face);
	// final calcuation //
	double derivativeValue = (coordinates[0]-centroid[0])*ey1derivative
							 +(coordinates[1]-centroid[1])*ey2derivative
							 +(coordinates[2]-centroid[2])*ey3derivative
							 +unityvector[2]*(ziderivative-1./6.);
	// returning value //
	//std::cout<<"Unit Y : "<<unityvector[0]<<unityvector[1]<<unityvector[3]<<std::endl;
	return derivativeValue;
	}
	
}

// -------------------------------------------------------------------------------------------------- //

double CentralisedDerivative::alphaZDerivative(Vertex* first, Vertex* second, Face* face){
	//gathering terms//
	//iterating for second vertex in the face
	double derivativeValue; //variable to store derivative
	unsigned int faceid = face->getID();
	Vertex * currentVertex;
	Edge * currentEdge;
	VertexEdgeIterator vertEdges(first);// iterating through all the edges coming out from vertex First
	unsigned int tempfaceid;
	Face * tempface;
	while((currentEdge = vertEdges.next())!=0){
		tempface = currentEdge->Left();
		tempfaceid = tempface->getID(); 
		if (tempfaceid == faceid){ // if the left face is same as Face* face, then we found the outgoing edge from First vert on the Face* face
			// currentEdge :  Org=firstvertex(j)----->-----Dest=vertex(j+1)
			Vertex* vertj1 = currentEdge->Dest();//the j+1 vertex
			double yj1 = vertj1->getYcoordinate();// y coordinate of j+1
			double yj = first->getYcoordinate(); // y coordinate of j
			derivativeValue = (yj-(face->getYCentralised()))*deltafunction(vertj1->getID(),second->getID())
							  + 1./6.*(yj1 - yj)
							  - (yj1 - face->getYCentralised())*deltafunction(first->getID(),second->getID());
			return derivativeValue;
		}
	}
	// if something is gone wrong
	return 0.;
}
// -------------------------------------------------------------------------------------------------- //

double CentralisedDerivative::betaZDerivative(Vertex* first, Vertex* second, Face* face){
	//gathering terms//
	//iterating for second vertex in the face
	double derivativeValue; //variable to store derivative
	unsigned int faceid = face->getID();
	Vertex * currentVertex;
	Edge * currentEdge;
	VertexEdgeIterator vertEdges(first);// iterating through all the edges coming out from vertex First
	unsigned int tempfaceid;
	Face * tempface;
	while((currentEdge = vertEdges.next())!=0){
		tempface = currentEdge->Left();
		tempfaceid = tempface->getID(); 
		if (tempfaceid == faceid){ // if the left face is same as Face* face, then we found the outgoing edge from First vert on the Face* face
			// currentEdge :  Org=firstvertex(j)----->-----Dest=vertex(j+1)
			Vertex* vertj1 = currentEdge->Dest();//the j+1 vertex
			double xj1 = vertj1->getXcoordinate();// y coordinate of j+1
			double xj = first->getXcoordinate(); // y coordinate of j
			derivativeValue = (xj1-(face->getXCentralised()))*deltafunction(first->getID(),second->getID())
							  + 1./6.*(xj - xj1)
							  - (xj - face->getXCentralised())*deltafunction(vertj1->getID(),second->getID());
			return derivativeValue;
		}
	}
	// if something is gone wrong
	return 0.;
}

// -------------------------------------------------------------------------------------------------- //

double CentralisedDerivative::gammaZDerivative(Vertex* first, Vertex* second, Face* face){
	return 0.0;
}
// -------------------------------------------------------------------------------------------------- //
double CentralisedDerivative::ex1ZDerivative(Vertex *first, Vertex *second, Face *face){
	//gathering terms // 
	double *pntpivector = face->getPivector();
	double *pntnormal = face->getNormal();
	double cPiNorm = sqrt(pow(pntpivector[0]-face->getXCentralised(),2)+
							pow(pntpivector[1]-face->getYCentralised(),2) + 
							pow(pntpivector[2]-face->getZCentralised(),2));
	//gathering derivatives
	double pixderivative = pixZDerivative(first, second, face); //calculating the ncxderivative
	double cPiNormderivative = cpinormZDerivative(first,second,face); //calculating the cpinorm X derivative
	//calculating the derivative value // 
	double derivativeValue = (cPiNorm*pixderivative-(pntpivector[0]-face->getXCentralised())*cPiNormderivative)/pow(cPiNorm,2.);
	return derivativeValue;
}
// -------------------------------------------------------------------------------------------------- //
double CentralisedDerivative::ex2ZDerivative(Vertex* first, Vertex* second, Face* face){
	// gathering the terms //
	double *pntpivector = face->getPivector();
	double *pntnormal = face->getNormal();
	double cpivector[3] = {pntpivector[0]-face->getXCentralised(),pntpivector[1]-face->getYCentralised(),pntpivector[2]-face->getZCentralised()};
	double cpinorm = sqrt(pow(cpivector[0],2)+
							pow(cpivector[1],2) + 
							pow(cpivector[2],2));
	double piderivative = piyZDerivative(first, second, face);
	double cpinormderivative = cpinormZDerivative(first, second, face);
	// calculating derivative 
	double derivativeValue = (1./pow(cpinorm,2))*(cpinorm*(piderivative) - (cpivector[1])*cpinormderivative);
	return derivativeValue;
}
// --------------------------------------------------------------------------------- //
double CentralisedDerivative::ex3ZDerivative(Vertex* first, Vertex* second, Face* face){
	// gathering the terms //
	double *pntpivector = face->getPivector();
	double *pntnormal = face->getNormal();
	double cpivector[3] = {pntpivector[0]-face->getXCentralised(),pntpivector[1]-face->getYCentralised(),pntpivector[2]-face->getZCentralised()};
	double cpinorm = sqrt(pow(cpivector[0],2)+
							pow(cpivector[1],2) + 
							pow(cpivector[2],2));
	double pizderivative = pizZDerivative(first, second, face);
	double cpinormderivative = cpinormZDerivative(first, second, face);
	// calculating derivative 
	double derivativeValue = (1./pow(cpinorm,2))*(cpinorm*(pizderivative-1./6.)- cpivector[2]*cpinormderivative);
	return derivativeValue;
}
// --------------------------------------------------------------------------------- //

double CentralisedDerivative::ey1ZDerivative(Vertex* first, Vertex* second, Face* face){
	//gathering terms // 
	double * unitxvector = face->getUnitx();
	double * unitzvector = face->getUnitz();
	double ncyderivative = ncyZDerivative(first, second, face);
	double nczderivative = nczZDerivative(first, second, face);
	double ex2derivative = ex2ZDerivative(first, second, face);
	double ex3derivative = ex3ZDerivative(first, second, face);
	// calculating terms
	double derivativevalue = unitxvector[1]*nczderivative+unitzvector[2]*ex2derivative-
							 unitxvector[2]*ncyderivative-unitzvector[1]*ex3derivative;
	//returning 
	return derivativevalue;
}
// --------------------------------------------------------------------------------- //
double CentralisedDerivative::ey2ZDerivative(Vertex* first, Vertex* second, Face* face){
	//gathering terms // 
	double * unitxvector = face->getUnitx();
	double * unitzvector = face->getUnitz();
	double ncxderivative = ncxZDerivative(first, second, face);
	double nczderivative = nczZDerivative(first, second, face);
	double ex1derivative = ex1ZDerivative(first, second, face);
	double ex3derivative = ex3ZDerivative(first, second, face);
	// calculating terms
	double derivativevalue = unitxvector[2]*ncxderivative+unitzvector[0]*ex3derivative-
							 unitxvector[0]*nczderivative-unitzvector[2]*ex1derivative;
	//returning 
	return derivativevalue;
}
// --------------------------------------------------------------------------------- //
double CentralisedDerivative::ey3ZDerivative(Vertex* first, Vertex* second, Face* face){
	//gathering terms // 
	double * unitxvector = face->getUnitx();
	double * unitzvector = face->getUnitz();
	double ncxderivative = ncxZDerivative(first, second, face);
	double ncyderivative = ncyZDerivative(first, second, face);
	double ex1derivative = ex1ZDerivative(first, second, face);
	double ex2derivative = ex2ZDerivative(first, second, face);
	// calculating terms
	double derivativevalue = unitxvector[0]*ncyderivative+unitzvector[1]*ex1derivative-
							 unitxvector[1]*ncxderivative-unitzvector[0]*ex2derivative;
	//returning 
	return derivativevalue;
}
// --------------------------------------------------------------------------------- //

double CentralisedDerivative::areatotalZDerivative(Vertex* first, Vertex* second, Face* face){
	// gathering terms //
	double derivativeValue(0), alphatemp(0), betatemp(0), gammatemp(0), normfactor(0);
	double alphaderivative(0),betaderivative(0), gammaderivative(0);
	unsigned int faceid = face->getID();
	FaceEdgeIterator faceEdges(face); // iterating the edges of this face
	Edge * currentEdge;
	Vertex* currentVertex;
	while ((currentEdge = faceEdges.next())!= 0){
			currentVertex = currentEdge->Dest(); //grabbing the destination vertex
			// gathering all terms //
			alphatemp = currentVertex->getAlpha(faceid);
			betatemp = currentVertex->getBeta(faceid);
			gammatemp = currentVertex->getGamma(faceid);
			normfactor = sqrt(pow(alphatemp,2)+pow(betatemp,2)+pow(gammatemp,2));
			alphaderivative = alphaZDerivative(currentVertex, second, face);
			betaderivative = betaZDerivative(currentVertex, second, face);
			gammaderivative = gammaZDerivative(currentVertex, second, face);
			//calculating derivative sum //
			derivativeValue += (1./normfactor)*(alphatemp*alphaderivative+betatemp*betaderivative+gammatemp*gammaderivative);
	}
	derivativeValue = (1./2.)*derivativeValue;
	return derivativeValue;
}

// -------------------------------------------------------------------------------------------------- //

double CentralisedDerivative::ziZDerivative(Vertex* first, Vertex* second, Face* face){
	return deltafunction(first->getID(),second->getID());
}
// -------------------------------------------------------------------------------------------------- //
double CentralisedDerivative::cpinormZDerivative(Vertex *first, Vertex *second, Face *face){
	// gathering terms
	double *pntpivector = face->getPivector();
	double *pntnormal = face->getNormal();
	double cPivector[3] = {pntpivector[0]-face->getXCentralised(),pntpivector[1]-face->getYCentralised(),pntpivector[2]-face->getZCentralised()};
	double cPiNorm = sqrt(pow(cPivector[0],2)+
							pow(cPivector[1],2) + 
							pow(cPivector[2],2));
	double pixderivative = pixZDerivative(first, second, face);
	double piyderivative = piyZDerivative(first, second, face);
	double pizderivative = pizZDerivative(first, second, face);
	// calculating the derivative value //
	double derivativeValue = (0.5/cPiNorm)*(2*(cPivector[0])*(pixderivative) +
											2*(cPivector[1])*(piyderivative) +
											2*(cPivector[2])*(pizderivative - 1./6.));
	return derivativeValue;
}

// -------------------------------------------------------------------------------------------------- //
double CentralisedDerivative::ncxZDerivative(Vertex* first, Vertex* second, Face* face){
	// gathering terms //
	double *pntnormaltilde = face->getNormalTilde();
	double normalTildeNorm = sqrt(pow(pntnormaltilde[0],2)+pow(pntnormaltilde[1],2)+pow(pntnormaltilde[2],2));
	double ncxtildederivative = ncxtildeZDerivative(first, second, face);
	double nctildenormderivative = nctildenormZDerivative(first, second, face);
	// calculating derivative //
	double derivativeValue = 1./(pow(normalTildeNorm,2))*(normalTildeNorm*ncxtildederivative
														 - pntnormaltilde[0]*nctildenormderivative);
	return derivativeValue;	
}
// -------------------------------------------------------------------------------------------------- //

double CentralisedDerivative::ncyZDerivative(Vertex* first, Vertex* second, Face* face){
	// gathering terms //
	double * normaltilde = face->getNormalTilde(); //getting the normal tilde vector of tha face
	double normaltildenorm = sqrt(pow(normaltilde[0],2)+pow(normaltilde[1],2)+pow(normaltilde[2],2)); // norm of normaltilde
	double ncytildederivative = ncytildeZDerivative(first, second, face);
	double normaltildenormderivative = nctildenormZDerivative(first, second, face);
	// calculating the derivative 
	double derivativeValue = (1./pow(normaltildenorm,2))*(normaltildenorm*ncytildederivative - normaltilde[1]*normaltildenormderivative);
	return derivativeValue;
}

// --------------------------------------------------------------------------------- //

double CentralisedDerivative::nczZDerivative(Vertex* first, Vertex* second, Face* face){
	// gathering terms //
	double * normaltilde = face->getNormalTilde(); //getting the normal tilde vector of tha face
	double normaltildenorm = sqrt(pow(normaltilde[0],2)+pow(normaltilde[1],2)+pow(normaltilde[2],2)); // norm of normaltilde
	double ncztildederivative = ncztildeZDerivative(first, second, face);
	double normaltildenormderivative = nctildenormZDerivative(first, second, face);
	// calculating the derivative 
	double derivativeValue = (1./pow(normaltildenorm,2))*(normaltildenorm*ncztildederivative - normaltilde[2]*normaltildenormderivative);
	return derivativeValue;
}
// -------------------------------------------------------------------------------------------------- //
double CentralisedDerivative::nctildenormZDerivative(Vertex* first, Vertex* second, Face* face){
		// gathering terms //
	double * normaltilde = face->getNormalTilde(); //getting the normal tilde vector of tha face
	double normaltildenorm = sqrt(pow(normaltilde[0],2)+pow(normaltilde[1],2)+pow(normaltilde[2],2));
	double ncxtildederivative = ncxtildeZDerivative(first, second, face);
	double ncytildederivative = ncytildeZDerivative(first, second, face);
	double ncztildederivative = ncztildeZDerivative(first, second, face);
	// calculating the derivative //
	double derivativeValue = (1./normaltildenorm) * (normaltilde[0]*ncxtildederivative +
													 normaltilde[1]*ncytildederivative +
													 normaltilde[2]*ncztildederivative);
	return derivativeValue;
}
// -------------------------------------------------------------------------------------------------- //
double CentralisedDerivative::ncxtildeZDerivative(Vertex* first, Vertex* second, Face* face){
	// gathering terms // 
	double areatotal = face->getAreaOfFace();
	double areatotalderivative = areatotalZDerivative(first, second, face); 
	// iterating the faces to get the all the alpha X derivatives
	double alphaderivative(0);
	double alphasum(0);
	unsigned int faceid = face->getID();//ID of the face
	Vertex* currentVertex; // a keeper for current Vertex
	Edge* currentEdge; // a keeper for current Edge
	FaceEdgeIterator faceEdges(face);//iterator to iterate through the vertex for outgoign edge
    while((currentEdge = faceEdges.next())!= 0){//runnign through the edges in the face again
     		currentVertex = currentEdge->Dest();//vertex to be derivatived by second vertex
     		// summing up alpha X derivatives of all vertices
     		alphaderivative += alphaZDerivative(currentVertex, second, face);
     		//summing up all the alpha values of the vertex
     		alphasum += currentVertex->getAlpha(faceid);
     	}
     // ******** Evaluating the value of derivative *********** //
     double derivativeValue = 0.5*(1./pow(areatotal,2))*(areatotal*alphaderivative - 
     													alphasum*areatotalderivative);
     return derivativeValue;
}

// -------------------------------------------------------------------------------------------------- //
double CentralisedDerivative::ncytildeZDerivative(Vertex* first, Vertex* second, Face* face){
	// gathering terms // 
	double areatotal = face->getAreaOfFace();
	double areatotalderivative = areatotalZDerivative(first, second, face); 
	// iterating the faces to get the all the alpha X derivatives
	double betaderivative(0);
	double betasum(0);
	unsigned int faceid = face->getID();//ID of the face
	Vertex* currentVertex; // a keeper for current Vertex
	Edge* currentEdge; // a keeper for current Edge
	FaceEdgeIterator faceEdges(face);//iterator to iterate through the vertex for outgoign edge
    while((currentEdge = faceEdges.next())!= 0){//runnign through the edges in the face again
     		currentVertex = currentEdge->Dest();//vertex to be derivatived by second vertex
     		// summing up alpha X derivatives of all vertices
     		betaderivative += betaZDerivative(currentVertex, second, face);
     		//summing up all the alpha values of the vertex
     		betasum += currentVertex->getBeta(faceid);
     	}
     // ******** Evaluating the value of derivative *********** //
     double derivativeValue = 0.5*(1./pow(areatotal,2))*(areatotal*betaderivative - 
     													betasum*areatotalderivative);
     return derivativeValue;
}


// --------------------------------------------------------------------------------- //
double CentralisedDerivative::ncztildeZDerivative(Vertex* first, Vertex* second, Face* face){
	// gathering terms // 
	double areatotal = face->getAreaOfFace();
	double areatotalderivative = areatotalZDerivative(first, second, face); 
	// iterating the faces to get the all the alpha X derivatives
	double gammaderivative(0);
	double gammasum(0);
	unsigned int faceid = face->getID();//ID of the face
	Vertex* currentVertex; // a keeper for current Vertex
	Edge* currentEdge; // a keeper for current Edge
	FaceEdgeIterator faceEdges(face);//iterator to iterate through the vertex for outgoign edge
    while((currentEdge = faceEdges.next())!= 0){//runnign through the edges in the face again
     		currentVertex = currentEdge->Dest();//vertex to be derivatived by second vertex
     		// summing up alpha X derivatives of all vertices
     		gammaderivative += gammaZDerivative(currentVertex, second, face);
     		//summing up all the alpha values of the vertex
     		gammasum += currentVertex->getGamma(faceid);
     	}
     // ******** Evaluating the value of derivative *********** //
     double derivativeValue = 0.5*(1./pow(areatotal,2))*(areatotal*gammaderivative - 
     													gammasum*areatotalderivative);
     return derivativeValue;
}



// --------------------------------------------------------------------------------- //

double CentralisedDerivative::pixZDerivative(Vertex* first, Vertex* second, Face* face){
	// gathering terms //
	double *pntnormal = face->getNormal();//grabing the normal of the face
	double ncxderivative = ncxZDerivative(first, second, face);
	//calcuating the derivative //
	double derivativeValue = -2.*pntnormal[0]*ncxderivative;
	return derivativeValue;
}
// --------------------------------------------------------------------------------- //


double CentralisedDerivative::piyZDerivative(Vertex* first, Vertex* second, Face* face){
	// gathering terms //
	double * normal = face->getNormal(); //grabbing the normal of the face
	double ncyderivative = ncyZDerivative(first, second, face);
	double ncxderivative = ncxZDerivative(first, second, face);
	// calculating the derivative // 
	double derivativeValue = -normal[0]*ncyderivative - normal[1]*ncxderivative;
	return derivativeValue;
}

// --------------------------------------------------------------------------------- //

double CentralisedDerivative::pizZDerivative(Vertex* first, Vertex* second, Face* face){
	// gathering terms //
	double * normal = face->getNormal(); //grabbing the normal of the face
	double nczderivative = nczZDerivative(first, second, face);
	double ncxderivative = ncxZDerivative(first, second, face);
	// calculating the derivative // 
	double derivativeValue = 1./6.- normal[0]*nczderivative - normal[2]*ncxderivative;
	return derivativeValue;
}
//  ------------------------------------------------------------------------------------//
//  ---------------------NUMERICAL DERIVATIVE CALCULATORS-------------------------//
// ------------------------------------------------------------------------------------//
// --------------------------------------------------------------------------------- //


// -----------------------X Derivatives--------------------------------------------- //
double CentralisedDerivative::numericalXtildeXDerivative(Vertex* first, Vertex* second, Face* face, double stepsize = pow(10.,-8.)){
	// ********************************************************************** //
	// calculating the numerical derivative with Newton method
	// f'(x) = (f(x+h)-f(x))/h
	// ********************************************************************** //
	// getting face ID 
	unsigned int faceid = face->getID();
	// getting the initial xtilde 
	double xtilde_initial = first->getProjectedXcoordinate(faceid);
	// now changing the second vertex by given stepsize and then setting new xtilde (x projected coordinate)
	double xsecond = second->getXcoordinate();
 	second->setXcoordinate(xsecond+stepsize);
 	face->setProjectedCoordinate();
 	// now new projected coordinate has been calculated with change in second vertex by (x+h)
 	double xtilde_final = first->getProjectedXcoordinate(faceid);
 	//now calculating the numerical derivative by Newtown method
 	double derivative = (xtilde_final-xtilde_initial)/stepsize;
 	// also setting the xcoordinate of second vertex back to previous value and projected coordinate also back to original
 	second->setXcoordinate(xsecond);
 	face->setProjectedCoordinate();
 	//return the derivative
 	return derivative;
}
// --------------------------------------------------------------------------------- //

double CentralisedDerivative::numericalYtildeXDerivative(Vertex* first, Vertex* second, Face* face, double stepsize = pow(10.,-8.)){
	// ********************************************************************** //
	// calculating the numerical derivative with Newton method
	// f'(x) = (f(x+h)-f(x))/h
	// ********************************************************************** //
	// getting face ID 
	unsigned int faceid = face->getID();
	// getting the initial xtilde 
	double ytilde_initial = first->getProjectedYcoordinate(faceid);
	// now changing the second vertex by given stepsize and then setting new xtilde (x projected coordinate)
	double xsecond = second->getXcoordinate();
 	second->setXcoordinate(xsecond+stepsize);
 	face->setProjectedCoordinate();
 	// now new projected coordinate has been calculated with change in second vertex by (x+h)
 	double ytilde_final = first->getProjectedYcoordinate(faceid);
 	//now calculating the numerical derivative by Newtown method
 	double derivative = (ytilde_final-ytilde_initial)/stepsize;
 	// also setting the xcoordinate of second vertex back to previous value and projected coordinate also back to original
 	second->setXcoordinate(xsecond);
 	face->setProjectedCoordinate();
 	//return the derivative
 	return derivative;
}

// -------------------  Y Derivatives ----------------------------------------- //
// --------------------------------------------------------------------------------- //
double CentralisedDerivative::numericalXtildeYDerivative(Vertex* first, Vertex* second, Face* face, double stepsize = pow(10.,-8.)){
	// ********************************************************************** //
	// calculating the numerical derivative with Newton method
	// f'(x) = (f(x+h)-f(x))/h
	// ********************************************************************** //
	// getting face ID 
	unsigned int faceid = face->getID();
	// getting the initial xtilde 
	double xtilde_initial = first->getProjectedXcoordinate(faceid);
	// now changing the second vertex by given stepsize and then setting new xtilde (x projected coordinate)
	double ysecond = second->getYcoordinate();
 	second->setYcoordinate(ysecond+stepsize);
 	face->setProjectedCoordinate();
 	// now new projected coordinate has been calculated with change in second vertex by (x+h)
 	double xtilde_final = first->getProjectedXcoordinate(faceid);
 	//now calculating the numerical derivative by Newtown method
 	double derivative = (xtilde_final-xtilde_initial)/stepsize;
 	// also setting the xcoordinate of second vertex back to previous value and projected coordinate also back to original
 	second->setYcoordinate(ysecond);
 	face->setProjectedCoordinate();
 	//return the derivative
 	return derivative;
}
// --------------------------------------------------------------------------------- //
double CentralisedDerivative::numericalYtildeYDerivative(Vertex* first, Vertex* second, Face* face, double stepsize = pow(10.,-8.)){
	// ********************************************************************** //
	// calculating the numerical derivative with Newton method
	// f'(x) = (f(x+h)-f(x))/h
	// ********************************************************************** //
	// getting face ID 
	unsigned int faceid = face->getID();
	// getting the initial xtilde 
	double ytilde_initial = first->getProjectedYcoordinate(faceid);
	// now changing the second vertex by given stepsize and then setting new xtilde (x projected coordinate)
	double ysecond = second->getYcoordinate();
 	second->setYcoordinate(ysecond+stepsize);
 	face->setProjectedCoordinate();
 	// now new projected coordinate has been calculated with change in second vertex by (x+h)
 	double ytilde_final = first->getProjectedYcoordinate(faceid);
 	//now calculating the numerical derivative by Newtown method
 	double derivative = (ytilde_final-ytilde_initial)/stepsize;
 	// also setting the xcoordinate of second vertex back to previous value and projected coordinate also back to original
 	second->setYcoordinate(ysecond);
 	face->setProjectedCoordinate();
 	//return the derivative
 	return derivative;
}
// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------- //

// -------------------------- Z Derivatives  -------------------------------------------- //

double CentralisedDerivative::numericalXtildeZDerivative(Vertex* first, Vertex* second, Face* face, double stepsize = pow(10.,-8.)){
	// ********************************************************************** //
	// calculating the numerical derivative with Newton method
	// f'(x) = (f(x+h)-f(x))/h
	// ********************************************************************** //
	// getting face ID 
	unsigned int faceid = face->getID();
	// getting the initial xtilde 
	double xtilde_initial = first->getProjectedXcoordinate(faceid);
	// now changing the second vertex by given stepsize and then setting new xtilde (x projected coordinate)
	double zsecond = second->getZcoordinate();
 	second->setZcoordinate(zsecond+stepsize);
 	face->setProjectedCoordinate();
 	// now new projected coordinate has been calculated with change in second vertex by (x+h)
 	double xtilde_final = first->getProjectedXcoordinate(faceid);
 	//now calculating the numerical derivative by Newtown method
 	double derivative = (xtilde_final-xtilde_initial)/stepsize;
 	// also setting the xcoordinate of second vertex back to previous value and projected coordinate also back to original
 	second->setZcoordinate(zsecond);
 	face->setProjectedCoordinate();
 	//return the derivative
 	return derivative;
}
// --------------------------------------------------------------------------------- //
double CentralisedDerivative::numericalYtildeZDerivative(Vertex* first, Vertex* second, Face* face, double stepsize = pow(10.,-8.)){
	// ********************************************************************** //
	// calculating the numerical derivative with Newton method
	// f'(x) = (f(x+h)-f(x))/h
	// ********************************************************************** //
	// getting face ID 
	unsigned int faceid = face->getID();
	// getting the initial xtilde 
	double ytilde_initial = first->getProjectedYcoordinate(faceid);
	// now changing the second vertex by given stepsize and then setting new xtilde (x projected coordinate)
	double zsecond = second->getZcoordinate();
 	second->setZcoordinate(zsecond+stepsize);
 	face->setProjectedCoordinate();
 	// now new projected coordinate has been calculated with change in second vertex by (x+h)
 	double ytilde_final = first->getProjectedYcoordinate(faceid);
 	//now calculating the numerical derivative by Newtown method
 	double derivative = (ytilde_final-ytilde_initial)/stepsize;
 	// also setting the xcoordinate of second vertex back to previous value and projected coordinate also back to original
 	second->setZcoordinate(zsecond);
 	face->setProjectedCoordinate();
 	//return the derivative
 	return derivative;
}
// ------------------------------------------------------------------------------------ //
double CentralisedDerivative::numericalFirstTermXDerivative(Vertex* second){
	//getting the initial x coordinate value
	double initialxk = second->getXcoordinate();
	//stepsize 
	double stepsize = pow(10.,-7.);
	std::cout<<"step size : "<<stepsize<<std::endl; //printing stepsize
	//iterating the faces of this vertex to collect initial value of first term
	VertexEdgeIterator edges(second); //iterating the edges around this vertex
	Edge * currentEdge;//pointer to current edge being worked on
	Face * currentFace;//pointer to current face
	double initialfirstvalues[3], finalfirstvalues[3];//array that stores the values of first term 
	unsigned int faceids[3];//storage for face ids
	int counter;
	//getting initial values
	counter = 0;
	while ((currentEdge = edges.next())!= 0){
		currentFace = currentEdge->Left();
		initialfirstvalues[counter] = currentFace->getFirstTerm();
		counter++;
	}
	// ******** Changin Coordinate ********** //
	//now changing the coordiante of the vertex
	second->setXcoordinate(initialxk + stepsize);
	//now setting the projected coordinates 
	VertexEdgeIterator edges1(second); //iterating the edges around this vertex
	while ((currentEdge = edges1.next())!= 0){
		currentFace = currentEdge->Left();
		currentFace->setProjectedCoordinate();
	}
	// ********** Getting new energy values ********* //
	//now setting the energy values
	VertexEdgeIterator edges2(second); //iterating the edges around this vertex
	counter = 0;
	while ((currentEdge = edges2.next())!= 0){
		currentFace = currentEdge->Left();
		currentFace->setEnergyTerms();
		faceids[counter] = currentFace->getID();
		finalfirstvalues[counter] = currentFace->getFirstTerm();
		counter++;
	}
	// *********** changing back the Coordinate ********** //
	//now changing the coordiante of the vertex
	second->setXcoordinate(initialxk);
	//now setting the projected coordinates 
	VertexEdgeIterator edges3(second); //iterating the edges around this vertex
	while ((currentEdge = edges3.next())!= 0){
		currentFace = currentEdge->Left();
		currentFace->setProjectedCoordinate();
		currentFace->setEnergyTerms();
	}
	//-------------------------------------------------------------------------- //
	//  Now calculating the derivatives // 
	double totalvalue = 0;
	for (int i = 0; i<counter ; i++){
		std::cout<<"face id "<<faceids[i]<<" first term initial : "<< finalfirstvalues[i]<<" final value : "<<initialfirstvalues[i]<<std::endl;
		totalvalue += (finalfirstvalues[i]-initialfirstvalues[i])/stepsize;
 	}
 	// final value // 
 	return totalvalue;
}



// ------------------------------------------------------------------------------------ //
double CentralisedDerivative::numericalSecondTermXDerivative(Vertex* second){
	//getting the initial x coordinate value
	double initialxk = second->getXcoordinate();
	//stepsize 
	double stepsize = pow(10.,-7.);
	std::cout<<"step size : "<<stepsize<<std::endl; //printing stepsize
	//iterating the faces of this vertex to collect initial value of first term
	VertexEdgeIterator edges(second); //iterating the edges around this vertex
	Edge * currentEdge;//pointer to current edge being worked on
	Face * currentFace;//pointer to current face
	double initialfirstvalues[3], finalfirstvalues[3];//array that stores the values of first term 
	unsigned int faceids[3];//storage for face ids
	int counter;
	//getting initial values
	counter = 0;
	while ((currentEdge = edges.next())!= 0){
		currentFace = currentEdge->Left();
		initialfirstvalues[counter] = currentFace->getSecondTerm();
		counter++;
	}
	// ******** Changin Coordinate ********** //
	//now changing the coordiante of the vertex
	second->setXcoordinate(initialxk + stepsize);
	//now setting the projected coordinates 
	VertexEdgeIterator edges1(second); //iterating the edges around this vertex
	while ((currentEdge = edges1.next())!= 0){
		currentFace = currentEdge->Left();
		currentFace->setProjectedCoordinate();
	}
	// ********** Getting new energy values ********* //
	//now setting the energy values
	VertexEdgeIterator edges2(second); //iterating the edges around this vertex
	counter = 0;
	while ((currentEdge = edges2.next())!= 0){
		currentFace = currentEdge->Left();
		currentFace->setEnergyTerms();
		faceids[counter] = currentFace->getID();
		finalfirstvalues[counter] = currentFace->getSecondTerm();
		counter++;
	}
	// *********** changing back the Coordinate ********** //
	//now changing the coordiante of the vertex
	second->setXcoordinate(initialxk);
	//now setting the projected coordinates 
	VertexEdgeIterator edges3(second); //iterating the edges around this vertex
	while ((currentEdge = edges3.next())!= 0){
		currentFace = currentEdge->Left();
		currentFace->setProjectedCoordinate();
		currentFace->setEnergyTerms();
	}
	//-------------------------------------------------------------------------- //
	//  Now calculating the derivatives // 
	double totalvalue = 0;
	for (int i = 0; i<counter ; i++){
		std::cout<<"first term initial : "<< finalfirstvalues[i]<<" final value : "<<initialfirstvalues[i]<<std::endl;
		totalvalue += (finalfirstvalues[i]-initialfirstvalues[i])/stepsize;
 	}
 	// final value // 
 	return totalvalue;
}

// ------------------------------------------------------------------------------------ //
double CentralisedDerivative::numericalThirdTermXDerivative(Vertex* second){
	//getting the initial x coordinate value
	double initialxk = second->getXcoordinate();
	//stepsize 
	double stepsize = pow(10.,-7.);
	std::cout<<"step size : "<<stepsize<<std::endl; //printing stepsize
	//iterating the faces of this vertex to collect initial value of first term
	VertexEdgeIterator edges(second); //iterating the edges around this vertex
	Edge * currentEdge;//pointer to current edge being worked on
	Face * currentFace;//pointer to current face
	double initialfirstvalues[3], finalfirstvalues[3];//array that stores the values of first term 
	unsigned int faceids[3];//storage for face ids
	int counter;
	//getting initial values
	counter = 0;
	while ((currentEdge = edges.next())!= 0){
		currentFace = currentEdge->Left();
		initialfirstvalues[counter] = currentFace->getThirdTerm();// *** CHANGE HERE *** //
		counter++;
	}
	// ******** Changin Coordinate ********** //
	//now changing the coordiante of the vertex
	second->setXcoordinate(initialxk + stepsize);
	//now setting the projected coordinates 
	VertexEdgeIterator edges1(second); //iterating the edges around this vertex
	while ((currentEdge = edges1.next())!= 0){
		currentFace = currentEdge->Left();
		currentFace->setProjectedCoordinate();
	}
	// ********** Getting new energy values ********* //
	//now setting the energy values
	VertexEdgeIterator edges2(second); //iterating the edges around this vertex
	counter = 0;
	while ((currentEdge = edges2.next())!= 0){
		currentFace = currentEdge->Left();
		currentFace->setEnergyTerms();
		faceids[counter] = currentFace->getID();
		finalfirstvalues[counter] = currentFace->getThirdTerm();// *** CHANGE HERE *** //
		counter++;
	}
	// *********** changing back the Coordinate ********** //
	//now changing the coordiante of the vertex
	second->setXcoordinate(initialxk);
	//now setting the projected coordinates 
	VertexEdgeIterator edges3(second); //iterating the edges around this vertex
	while ((currentEdge = edges3.next())!= 0){
		currentFace = currentEdge->Left();
		currentFace->setProjectedCoordinate();
		currentFace->setEnergyTerms();
	}
	//-------------------------------------------------------------------------- //
	//  Now calculating the derivatives // 
	double totalvalue = 0;
	for (int i = 0; i<counter ; i++){
		std::cout<<"first term initial : "<< finalfirstvalues[i]<<" final value : "<<initialfirstvalues[i]<<std::endl;
		totalvalue += (finalfirstvalues[i]-initialfirstvalues[i])/stepsize;
 	}
 	// final value // 
 	return totalvalue;
}



// ------------------------------------------------------------------------------------ //
double CentralisedDerivative::numericalFirstTermYDerivative(Vertex* second){
	//getting the initial x coordinate value
	double initialxk = second->getYcoordinate();
	//stepsize 
	double stepsize = pow(10.,-7.);
	std::cout<<"step size : "<<stepsize<<std::endl; //printing stepsize
	//iterating the faces of this vertex to collect initial value of first term
	VertexEdgeIterator edges(second); //iterating the edges around this vertex
	Edge * currentEdge;//pointer to current edge being worked on
	Face * currentFace;//pointer to current face
	double initialfirstvalues[3], finalfirstvalues[3];//array that stores the values of first term 
	unsigned int faceids[3];//storage for face ids
	int counter;
	//getting initial values
	counter = 0;
	while ((currentEdge = edges.next())!= 0){
		currentFace = currentEdge->Left();
		initialfirstvalues[counter] = currentFace->getFirstTerm();
		counter++;
	}
	// ******** Changin Coordinate ********** //
	//now changing the coordiante of the vertex
	second->setYcoordinate(initialxk + stepsize);
	//now setting the projected coordinates 
	VertexEdgeIterator edges1(second); //iterating the edges around this vertex
	while ((currentEdge = edges1.next())!= 0){
		currentFace = currentEdge->Left();
		currentFace->setProjectedCoordinate();
	}
	// ********** Getting new energy values ********* //
	//now setting the energy values
	VertexEdgeIterator edges2(second); //iterating the edges around this vertex
	counter = 0;
	while ((currentEdge = edges2.next())!= 0){
		currentFace = currentEdge->Left();
		currentFace->setEnergyTerms();
		faceids[counter] = currentFace->getID();
		finalfirstvalues[counter] = currentFace->getFirstTerm();
		counter++;
	}
	// *********** changing back the Coordinate ********** //
	//now changing the coordiante of the vertex
	second->setYcoordinate(initialxk);
	//now setting the projected coordinates 
	VertexEdgeIterator edges3(second); //iterating the edges around this vertex
	while ((currentEdge = edges3.next())!= 0){
		currentFace = currentEdge->Left();
		currentFace->setProjectedCoordinate();
		currentFace->setEnergyTerms();
	}
	//-------------------------------------------------------------------------- //
	//  Now calculating the derivatives // 
	double totalvalue = 0;
	for (int i = 0; i<counter ; i++){
		std::cout<<"face id "<<faceids[i]<<" first term initial : "<< finalfirstvalues[i]<<" final value : "<<initialfirstvalues[i]<<std::endl;
		totalvalue += (finalfirstvalues[i]-initialfirstvalues[i])/stepsize;
 	}
 	// final value // 
 	return totalvalue;
}



// ------------------------------------------------------------------------------------ //
double CentralisedDerivative::numericalSecondTermYDerivative(Vertex* second){
	//getting the initial x coordinate value
	double initialxk = second->getYcoordinate();
	//stepsize 
	double stepsize = pow(10.,-7.);
	std::cout<<"step size : "<<stepsize<<std::endl; //printing stepsize
	//iterating the faces of this vertex to collect initial value of first term
	VertexEdgeIterator edges(second); //iterating the edges around this vertex
	Edge * currentEdge;//pointer to current edge being worked on
	Face * currentFace;//pointer to current face
	double initialfirstvalues[3], finalfirstvalues[3];//array that stores the values of first term 
	unsigned int faceids[3];//storage for face ids
	int counter;
	//getting initial values
	counter = 0;
	while ((currentEdge = edges.next())!= 0){
		currentFace = currentEdge->Left();
		initialfirstvalues[counter] = currentFace->getSecondTerm();
		counter++;
	}
	// ******** Changin Coordinate ********** //
	//now changing the coordiante of the vertex
	second->setYcoordinate(initialxk + stepsize);
	//now setting the projected coordinates 
	VertexEdgeIterator edges1(second); //iterating the edges around this vertex
	while ((currentEdge = edges1.next())!= 0){
		currentFace = currentEdge->Left();
		currentFace->setProjectedCoordinate();
	}
	// ********** Getting new energy values ********* //
	//now setting the energy values
	VertexEdgeIterator edges2(second); //iterating the edges around this vertex
	counter = 0;
	while ((currentEdge = edges2.next())!= 0){
		currentFace = currentEdge->Left();
		currentFace->setEnergyTerms();
		faceids[counter] = currentFace->getID();
		finalfirstvalues[counter] = currentFace->getSecondTerm();
		counter++;
	}
	// *********** changing back the Coordinate ********** //
	//now changing the coordiante of the vertex
	second->setYcoordinate(initialxk);
	//now setting the projected coordinates 
	VertexEdgeIterator edges3(second); //iterating the edges around this vertex
	while ((currentEdge = edges3.next())!= 0){
		currentFace = currentEdge->Left();
		currentFace->setProjectedCoordinate();
		currentFace->setEnergyTerms();
	}
	//-------------------------------------------------------------------------- //
	//  Now calculating the derivatives // 
	double totalvalue = 0;
	for (int i = 0; i<counter ; i++){
		std::cout<<"first term initial : "<< finalfirstvalues[i]<<" final value : "<<initialfirstvalues[i]<<std::endl;
		totalvalue += (finalfirstvalues[i]-initialfirstvalues[i])/stepsize;
 	}
 	// final value // 
 	return totalvalue;
}

// ------------------------------------------------------------------------------------ //
double CentralisedDerivative::numericalThirdTermYDerivative(Vertex* second){
	//getting the initial x coordinate value
	double initialxk = second->getYcoordinate();
	//stepsize 
	double stepsize = pow(10.,-7.);
	std::cout<<"step size : "<<stepsize<<std::endl; //printing stepsize
	//iterating the faces of this vertex to collect initial value of first term
	VertexEdgeIterator edges(second); //iterating the edges around this vertex
	Edge * currentEdge;//pointer to current edge being worked on
	Face * currentFace;//pointer to current face
	double initialfirstvalues[3], finalfirstvalues[3];//array that stores the values of first term 
	unsigned int faceids[3];//storage for face ids
	int counter;
	//getting initial values
	counter = 0;
	while ((currentEdge = edges.next())!= 0){
		currentFace = currentEdge->Left();
		initialfirstvalues[counter] = currentFace->getThirdTerm();// *** CHANGE HERE *** //
		counter++;
	}
	// ******** Changin Coordinate ********** //
	//now changing the coordiante of the vertex
	second->setYcoordinate(initialxk + stepsize);
	//now setting the projected coordinates 
	VertexEdgeIterator edges1(second); //iterating the edges around this vertex
	while ((currentEdge = edges1.next())!= 0){
		currentFace = currentEdge->Left();
		currentFace->setProjectedCoordinate();
	}
	// ********** Getting new energy values ********* //
	//now setting the energy values
	VertexEdgeIterator edges2(second); //iterating the edges around this vertex
	counter = 0;
	while ((currentEdge = edges2.next())!= 0){
		currentFace = currentEdge->Left();
		currentFace->setEnergyTerms();
		faceids[counter] = currentFace->getID();
		finalfirstvalues[counter] = currentFace->getThirdTerm();// *** CHANGE HERE *** //
		counter++;
	}
	// *********** changing back the Coordinate ********** //
	//now changing the coordiante of the vertex
	second->setYcoordinate(initialxk);
	//now setting the projected coordinates 
	VertexEdgeIterator edges3(second); //iterating the edges around this vertex
	while ((currentEdge = edges3.next())!= 0){
		currentFace = currentEdge->Left();
		currentFace->setProjectedCoordinate();
		currentFace->setEnergyTerms();
	}
	//-------------------------------------------------------------------------- //
	//  Now calculating the derivatives // 
	double totalvalue = 0;
	for (int i = 0; i<counter ; i++){
		std::cout<<"first term initial : "<< finalfirstvalues[i]<<" final value : "<<initialfirstvalues[i]<<std::endl;
		totalvalue += (finalfirstvalues[i]-initialfirstvalues[i])/stepsize;
 	}
 	// final value // 
 	return totalvalue;
}





// ------------------------------------------------------------------------------------ //
double CentralisedDerivative::numericalFirstTermZDerivative(Vertex* second){
	//getting the initial x coordinate value
	double initialxk = second->getZcoordinate();
	//stepsize 
	double stepsize = pow(10.,-7.);
	std::cout<<"step size : "<<stepsize<<std::endl; //printing stepsize
	//iterating the faces of this vertex to collect initial value of first term
	VertexEdgeIterator edges(second); //iterating the edges around this vertex
	Edge * currentEdge;//pointer to current edge being worked on
	Face * currentFace;//pointer to current face
	double initialfirstvalues[3], finalfirstvalues[3];//array that stores the values of first term 
	unsigned int faceids[3];//storage for face ids
	int counter;
	//getting initial values
	counter = 0;
	while ((currentEdge = edges.next())!= 0){
		currentFace = currentEdge->Left();
		initialfirstvalues[counter] = currentFace->getFirstTerm();
		counter++;
	}
	// ******** Changin Coordinate ********** //
	//now changing the coordiante of the vertex
	second->setZcoordinate(initialxk + stepsize);
	//now setting the projected coordinates 
	VertexEdgeIterator edges1(second); //iterating the edges around this vertex
	while ((currentEdge = edges1.next())!= 0){
		currentFace = currentEdge->Left();
		currentFace->setProjectedCoordinate();
	}
	// ********** Getting new energy values ********* //
	//now setting the energy values
	VertexEdgeIterator edges2(second); //iterating the edges around this vertex
	counter = 0;
	while ((currentEdge = edges2.next())!= 0){
		currentFace = currentEdge->Left();
		currentFace->setEnergyTerms();
		faceids[counter] = currentFace->getID();
		finalfirstvalues[counter] = currentFace->getFirstTerm();
		counter++;
	}
	// *********** changing back the Coordinate ********** //
	//now changing the coordiante of the vertex
	second->setZcoordinate(initialxk);
	//now setting the projected coordinates 
	VertexEdgeIterator edges3(second); //iterating the edges around this vertex
	while ((currentEdge = edges3.next())!= 0){
		currentFace = currentEdge->Left();
		currentFace->setProjectedCoordinate();
		currentFace->setEnergyTerms();
	}
	//-------------------------------------------------------------------------- //
	//  Now calculating the derivatives // 
	double totalvalue = 0;
	for (int i = 0; i<counter ; i++){
		std::cout<<"face id "<<faceids[i]<<" first term initial : "<< finalfirstvalues[i]<<" final value : "<<initialfirstvalues[i]<<std::endl;
		totalvalue += (finalfirstvalues[i]-initialfirstvalues[i])/stepsize;
 	}
 	// final value // 
 	return totalvalue;
}



// ------------------------------------------------------------------------------------ //
double CentralisedDerivative::numericalSecondTermZDerivative(Vertex* second){
	//getting the initial x coordinate value
	double initialxk = second->getZcoordinate();
	//stepsize 
	double stepsize = pow(10.,-7.);
	std::cout<<"step size : "<<stepsize<<std::endl; //printing stepsize
	//iterating the faces of this vertex to collect initial value of first term
	VertexEdgeIterator edges(second); //iterating the edges around this vertex
	Edge * currentEdge;//pointer to current edge being worked on
	Face * currentFace;//pointer to current face
	double initialfirstvalues[3], finalfirstvalues[3];//array that stores the values of first term 
	unsigned int faceids[3];//storage for face ids
	int counter;
	//getting initial values
	counter = 0;
	while ((currentEdge = edges.next())!= 0){
		currentFace = currentEdge->Left();
		initialfirstvalues[counter] = currentFace->getSecondTerm();
		counter++;
	}
	// ******** Changin Coordinate ********** //
	//now changing the coordiante of the vertex
	second->setZcoordinate(initialxk + stepsize);
	//now setting the projected coordinates 
	VertexEdgeIterator edges1(second); //iterating the edges around this vertex
	while ((currentEdge = edges1.next())!= 0){
		currentFace = currentEdge->Left();
		currentFace->setProjectedCoordinate();
	}
	// ********** Getting new energy values ********* //
	//now setting the energy values
	VertexEdgeIterator edges2(second); //iterating the edges around this vertex
	counter = 0;
	while ((currentEdge = edges2.next())!= 0){
		currentFace = currentEdge->Left();
		currentFace->setEnergyTerms();
		faceids[counter] = currentFace->getID();
		finalfirstvalues[counter] = currentFace->getSecondTerm();
		counter++;
	}
	// *********** changing back the Coordinate ********** //
	//now changing the coordiante of the vertex
	second->setZcoordinate(initialxk);
	//now setting the projected coordinates 
	VertexEdgeIterator edges3(second); //iterating the edges around this vertex
	while ((currentEdge = edges3.next())!= 0){
		currentFace = currentEdge->Left();
		currentFace->setProjectedCoordinate();
		currentFace->setEnergyTerms();
	}
	//-------------------------------------------------------------------------- //
	//  Now calculating the derivatives // 
	double totalvalue = 0;
	for (int i = 0; i<counter ; i++){
		std::cout<<"first term initial : "<< finalfirstvalues[i]<<" final value : "<<initialfirstvalues[i]<<std::endl;
		totalvalue += (finalfirstvalues[i]-initialfirstvalues[i])/stepsize;
 	}
 	// final value // 
 	return totalvalue;
}

// ------------------------------------------------------------------------------------ //
double CentralisedDerivative::numericalThirdTermZDerivative(Vertex* second){
	//getting the initial x coordinate value
	double initialxk = second->getZcoordinate();
	//stepsize 
	double stepsize = pow(10.,-7.);
	std::cout<<"step size : "<<stepsize<<std::endl; //printing stepsize
	//iterating the faces of this vertex to collect initial value of first term
	VertexEdgeIterator edges(second); //iterating the edges around this vertex
	Edge * currentEdge;//pointer to current edge being worked on
	Face * currentFace;//pointer to current face
	double initialfirstvalues[3], finalfirstvalues[3];//array that stores the values of first term 
	unsigned int faceids[3];//storage for face ids
	int counter;
	//getting initial values
	counter = 0;
	while ((currentEdge = edges.next())!= 0){
		currentFace = currentEdge->Left();
		initialfirstvalues[counter] = currentFace->getThirdTerm();// *** CHANGE HERE *** //
		counter++;
	}
	// ******** Changin Coordinate ********** //
	//now changing the coordiante of the vertex
	second->setZcoordinate(initialxk + stepsize);
	//now setting the projected coordinates 
	VertexEdgeIterator edges1(second); //iterating the edges around this vertex
	while ((currentEdge = edges1.next())!= 0){
		currentFace = currentEdge->Left();
		currentFace->setProjectedCoordinate();
	}
	// ********** Getting new energy values ********* //
	//now setting the energy values
	VertexEdgeIterator edges2(second); //iterating the edges around this vertex
	counter = 0;
	while ((currentEdge = edges2.next())!= 0){
		currentFace = currentEdge->Left();
		currentFace->setEnergyTerms();
		faceids[counter] = currentFace->getID();
		finalfirstvalues[counter] = currentFace->getThirdTerm();// *** CHANGE HERE *** //
		counter++;
	}
	// *********** changing back the Coordinate ********** //
	//now changing the coordiante of the vertex
	second->setZcoordinate(initialxk);
	//now setting the projected coordinates 
	VertexEdgeIterator edges3(second); //iterating the edges around this vertex
	while ((currentEdge = edges3.next())!= 0){
		currentFace = currentEdge->Left();
		currentFace->setProjectedCoordinate();
		currentFace->setEnergyTerms();
	}
	//-------------------------------------------------------------------------- //
	//  Now calculating the derivatives // 
	double totalvalue = 0;
	for (int i = 0; i<counter ; i++){
		std::cout<<"first term initial : "<< finalfirstvalues[i]<<" final value : "<<initialfirstvalues[i]<<std::endl;
		totalvalue += (finalfirstvalues[i]-initialfirstvalues[i])/stepsize;
 	}
 	// final value // 
 	return totalvalue;
}





















