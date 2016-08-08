
#include "cell.hh"
#include "edge.hh"
#include "vertex.hh"
#include "face.hh"
#include "derivative.hh"
#include <math.h>

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


// -----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%//
// -----%-----%-----%----Private Members Only-----%-----%-----%-----%-----%-----%-----%//
// -----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%//

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
double CentralisedDerivative::pixXDerivative(Vertex* first, Vertex* second, Face* face){
	// gathering terms //
	double *pntnormal = face->getNormal();//grabing the normal of the face
	double ncxderivative = ncxXDerivative(first, second, face);
	//calcuating the derivative //
	double derivativeValue = 1./6. - 2*pntnormal[0]*ncxderivative;
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
	derivativeValue = (1./(12.))*derivativeValue;
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
























