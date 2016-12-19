
/* ============================================================================
 * p2/cell/face.cc
 * ========================================================================= */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>//for math operations sqrt and others
#include <algorithm>//for std::max
#include "cell.hh"
#include "edge.hh"
#include "face.hh"
#include "vertex.hh"

/* ----------------------------------------------------------------------------
 * Face
 * ------------------------------------------------------------------------- */

/* -- public class methods ------------------------------------------------- */

Face *Face::make(Cell *cell)
{
  assert(cell!=0);
  return new Face(cell);
}

void Face::kill(Face *face)
{
  assert(face!=0);

  delete face;
}

/* -- public instance methods ---------------------------------------------- */
void Face::setID(unsigned int id)
{
  assert(id>0);

  this->id = id;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
void Face::addEdge(Edge *edge)
{
  assert(edge!=0);

  // only keep track of one edge in the orbit--this one is as good as any

  this->edge = edge;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
void Face::removeEdge(Edge *edge)
{
  assert(edge!=0);

  // replace the arbitrary edge with another edge in the orbit
  // use null if this is the only edge
  // assumes that the edge hasn't been actually removed yet

  Edge *next = edge->Onext();

  this->edge = next!=edge ? next : 0;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
void Face::setCentralisedCoordinate(double xcent, double ycent, double zcent){
  
  this->xCentralised = xcent;
  this->yCentralised = ycent; 
  this->zCentralised = zcent;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
void Face::setProjectedCoordinate(){
  //setting the projected coordinates of all vertices in this face
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
    // getting the mean center position of this face //
    unsigned int faceid = this->getID();// first grabbing id of  current face
    if (faceid == 1){//if faceid == 1 then this is external face, that doesnot need projection or calculation
        // setting the central coordinate of this face in terms of cartisian coordinate
          this->xCentralised = 0.;
          this->yCentralised = 0.;
          this->zCentralised = 0.;
        return;
    } 
    FaceEdgeIterator faceEdges(this);//iterator to iterate through the vertex for outgoign edge
    double xCentroid(0), yCentroid(0), zCentroid(0); // coordinate of the cnetroid
    // array of vertices
    double xcood[6], ycood[6], zcood[6];// coordinate of the vertices of this face
    double xTriCen[6], yTriCen[6], zTriCen[6], areaTri;// coordinate of center of triangluated triangles of face
    //areaTri : area of triangles
    double xmean(0), ymean(0), zmean(0);//mean center of the face
    int counter;//counter to keep track of vertices
    Edge *currentEdge;//pointer to keep track of a edge that is iterated 
    Vertex *currentVertex; // vertex to store the dest() of currentEdge
    counter = 0;
    while ((currentEdge = faceEdges.next())!=0){//iterating the edges in the face
          currentVertex = currentEdge->Dest();
          // getting the vertices of the current face and storring it in array
          xcood[counter] = currentVertex->getXcoordinate();
          ycood[counter] = currentVertex->getYcoordinate();
          zcood[counter] = currentVertex->getZcoordinate();
          xmean += xcood[counter];
          ymean += ycood[counter];
          zmean += zcood[counter];
          counter += 1;//increasing the counter value
    }
    //printf(" face id : %u \n", faceid);
    //printf(" calculated means : X = %F ; Y = %F ; Z = %F \n", xmean, ymean,zmean );
    //printf("counter %h \n",counter);
    // divinding by the number of vertices to get the mean, which is also the centroid
    xCentroid = (1./counter)*xmean;
    yCentroid = (1./counter)*ymean;
    zCentroid = (1./counter)*zmean;
    // setting the central coordinate of this face in terms of cartisian coordinate
    this->xCentralised = xCentroid;
    this->yCentralised = yCentroid;
    this->zCentralised = zCentroid;
    //printf("%s %u %F \n ", "calculating total area of Face :",faceid, totalarea);
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
    //Calculating normal of the triangles of the shape and then calculating weighted normal
    //double normalTempX[6],normalTempY[6],normalTempZ[6],totalarea;//vertices of normal of triangles
    double vector1[3], vector2[3], crossProductVector[3],crossProductMagnitude;
    double totalarea = 0;
    double normalX(0), normalY(0), normalZ(0); // vertices of weighted normal X, Y, Z 
    // iterating the faces also
    FaceEdgeIterator faceEdges1(this);//iterator to iterate through the vertex for outgoign edge
    // curentVertex and currentEdge have been already declared above
    currentEdge = faceEdges1.next();// grabbing the one edge of face = this
    for (int counter = 0; counter<6; counter++){//triangulating the vertices with centroid of this face
          currentVertex = currentEdge->Dest();//getting the destination of the edge
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
          // now saving the above calculation of alpha, beta, gamma for the vertex
          currentVertex->setAlpha(faceid, crossProductVector[0]);
          currentVertex->setBeta(faceid, crossProductVector[1]);
          currentVertex->setGamma(faceid, crossProductVector[2]);
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
    //setting area of the face
    // this totalarea is not the exact area of face, 
    //so instead calling setareaofface fucntion
    // however, this area is the total area for derivative calculation
    this->areaOfFace = totalarea;
    //this->setAreaOfFace(); // to set the area
    {//right now normalx,y,z is un-normalised so it is normalTilde
      //setting normalTilde
        double normaltildeOfFace[3] = {normalX, normalY, normalZ};
        double * pntnormaltilde = normaltildeOfFace;
        this->setNormalTilde(pntnormaltilde);
    }
    //normalising the normal vector 
    double normalMagnitude = sqrt(pow(normalX,2.0)+pow(normalY,2.0)+pow(normalZ,2.0));
    normalX = normalX/normalMagnitude;
    normalY = normalY/normalMagnitude;
    normalZ = normalZ/normalMagnitude;
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
    double normalOfFace[3] = {normalX, normalY, normalZ};
    double *pntNormalOfFace = normalOfFace;
    this->setNormal(pntNormalOfFace);//saving the normal of this face
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
    //now need to get this projected 3d plane coordinates onto rotated 2d coordinate system
    // getting the unit vector of 2d Plane
    double unitx[3], unity[3], vectorVertex[3];//new unit vector of x, y on the plane
    double dotproduct, xvertex, yvertex, zvertex;
    //to get the unit vector in X direction, lets project X axis (1,0,0) on to the plane
    //vector from (1,0,0)+(xCentroid,Ycentroid,Zcentroid) until xCentroid
    vectorVertex[0] = 1.;//1.+xCentroid-xCentroid; just directly writing the resulting answer
    vectorVertex[1] = 0.;//0.+yCentroid-yCentroid;
    vectorVertex[2] = 0.;//0.+zCentroid-zCentroid;
    //dot product of vectorVertex and normal, which is equal to normalX; this is n^c_x = \phi_P
    dotproduct = normalX;//vectorVertex[0]*normalX+vectorVertex[1]*normalY+vectorVertex[2]*normalZ;
    //now calculating the projected vertices --at this stage unitx is storing pi_vector
    unitx[0]  = 1.0 +xCentroid- dotproduct*normalX;
    unitx[1] = 0.+yCentroid - dotproduct*normalY;
    unitx[2] = 0.+zCentroid - dotproduct*normalZ;
    {//saving the pi_vector of this face
    double * pntPivector = unitx;
    this->setPivector(pntPivector);
    }
    //getting the unitx = normalised[Projectedvertex1 - Origin]
    unitx[0] = unitx[0]- xCentroid;
    unitx[1] = unitx[1]- yCentroid;
    unitx[2] = unitx[2]- zCentroid;
    // normalising unitx
    double normUnitx   = sqrt(pow(unitx[0],2)+pow(unitx[1],2)+pow(unitx[2],2));
    unitx[0] = unitx[0]/normUnitx;
    unitx[1] = unitx[1]/normUnitx;
    unitx[2] = unitx[2]/normUnitx;
    // now getting unity : it is the cross product of  normal and unitx to the plane. 
    unity[0] = unitx[1]*normalZ-unitx[2]*normalY;
    unity[1] = unitx[2]*normalX-unitx[0]*normalZ;
    unity[2] = unitx[0]*normalY-unitx[1]*normalX;
    //normalising unity
    double normUnity = sqrt(pow(unity[0],2)+pow(unity[1],2)+pow(unity[2],2));
    unity[0] = unity[0]/normUnity;
    unity[1] = unity[1]/normUnity;
    unity[2] = unity[2]/normUnity;
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
    // now projecting the 3d coordinates to the 2D coordinates and assinging the vertices
    double xprojection, yprojection, zprojection;
    FaceEdgeIterator faceEdges2(this);//iterator to iterate through the vertex for outgoign edge
    while((currentEdge = faceEdges2.next())!= 0){//runnign through the edges in the face again
          currentVertex = currentEdge->Dest();// the current vertex, iterated
          xvertex = currentVertex->getXcoordinate();
          yvertex = currentVertex->getYcoordinate();
          zvertex = currentVertex->getZcoordinate();
          //getting the vector form the centroid or new origin
          // CENTRALISED PorjectedCoordinate
          vectorVertex[0] = xvertex-xCentroid;
          vectorVertex[1] = yvertex-yCentroid;
          vectorVertex[2] = zvertex-zCentroid;
          //now getting the new x, y coordinates, dot product of unit vector with the vectorVertex
          xprojection = unitx[0]*vectorVertex[0]+unitx[1]*vectorVertex[1]+unitx[2]*vectorVertex[2];
          yprojection = unity[0]*vectorVertex[0]+unity[1]*vectorVertex[1]+unity[2]*vectorVertex[2];
          zprojection = normalX*vectorVertex[0]+normalY*vectorVertex[1]+normalZ*vectorVertex[2];
          //now setting the projected coordinates in the vertex properties
          currentVertex->insertProjectedXcoordinate(faceid,xprojection);
          currentVertex->insertProjectedYcoordinate(faceid,yprojection); 
          currentVertex->insertProjectedZcoordinate(faceid,zprojection); 
          //getting the vector form the real origin on new vertex
          // NON CENTRALISED Projected Coordainte
          vectorVertex[0] = xvertex;
          vectorVertex[1] = yvertex;
          vectorVertex[2] = zvertex;
          //now getting the new x, y coordinates, dot product of unit vector with the vectorVertex
          xprojection = unitx[0]*vectorVertex[0]+unitx[1]*vectorVertex[1]+unitx[2]*vectorVertex[2];
          yprojection = unity[0]*vectorVertex[0]+unity[1]*vectorVertex[1]+unity[2]*vectorVertex[2];
          zprojection = normalX*vectorVertex[0]+normalY*vectorVertex[1]+normalZ*vectorVertex[2];
          //now setting the projected coordinates in the vertex properties
          currentVertex->insertNonCentralisedProjectedXcoordinate(faceid,xprojection);
          currentVertex->insertNonCentralisedProjectedYcoordinate(faceid,yprojection); 
          currentVertex->insertNonCentralisedProjectedZcoordinate(faceid,zprojection); 
        }
    //saving unit vectors
   {
    double * pntunitx = unitx;
    double * pntunity = unity;
    double * pntunitz = normalOfFace;
    this->setUnitx(pntunitx);
    this->setUnity(pntunity);
    this->setUnitz(pntunitz);
    }
    //Calculating the Angle of tilt to the cartesian x axis
    // dor product between the unitx vector and cartesian x unit vector [1,0,0]
    dotproduct = unitx[0]; // as second and third terms are multiplied by 0
    double theta = acos(dotproduct);//since both unitx and x unit vector have 1 magnitude, just dividing by 1
    this->setAngleOfTilt(theta);
    // setting Mu values 
    //setting Mu values 
    this->setMu();
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
void Face::setNormal(double * tempnormal){
   for (int i = 0; i<3; i++){
      normal[i] = tempnormal[i];
   }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
double *Face::getNormal(){
   double * pntnormal = normal;
   return pntnormal;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
void Face::setNormalTilde(double * tempnormal){
   for (int i = 0; i<3; i++){
      normaltilde[i] = tempnormal[i];
   }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
double *Face::getNormalTilde(){
   double * pntnormal = normaltilde;
   return pntnormal;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
void Face::setPivector(double * tempnormal){
   for (int i = 0; i<3; i++){
      pivector[i] = tempnormal[i];
   }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
double *Face::getPivector(){
   double * pntnormal = pivector;
   return pntnormal;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
void Face::setUnitx(double * tempunit){
  for (int i = 0; i<3; i++){
      unitx[i] = tempunit[i];
   }
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
double * Face::getUnitx(){
  double * pntunit = unitx;
  return pntunit;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
void Face::setUnity(double * tempunit){
  for (int i = 0; i<3; i++){
      unity[i] = tempunit[i];
   }
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
double * Face::getUnity(){
  double * pntunit = unity;
  return pntunit;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
void Face::setUnitz(double * tempunit){
  for (int i = 0; i<3; i++){
      unitz[i] = tempunit[i];
   }
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
double * Face::getUnitz(){
  double * pntunit = unitz;
  return pntunit;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
double Face::getAreaOfFace(){
	return this->areaOfFace;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
void Face::setAreaOfFace(){
	FaceEdgeIterator edges(this);//iterator to iterate through the face
	Edge *newedge;//a pointer to keep track of current edge
	double areasum = 0;//variable to store the area of the face
	unsigned int faceid = this->id;//getting the faceid
  Vertex *vertexOrg;
	//looping through all the edges in the face until exhausted
	while ((newedge = edges.next())!=0){
    vertexOrg = newedge->Org();
		areasum += vertexOrg->getAk(faceid);//summing up all the Ak values for this face
	}
	this->areaOfFace = 0.5*areasum;//storing the area of the face in areaOfFace variable
}
//******************added features********************************//
 void Face::addVertex(Vertex *vertex)
   {
    assert(vertex!=0);

  // expand the vertex array, if necessary

  if (vertexCount>=vertexSize)
  {
    unsigned int vertexSizeNew = vertexSize*2;
    Vertex     **verticesNew   = new Vertex*[vertexSizeNew];

    for (unsigned int i = 0; i<vertexCount; i++)
      verticesNew[i] = vertices[i];

    delete[] vertices;

    vertices   = verticesNew;
    vertexSize = vertexSizeNew;
  }
  
  // add the vertex to the end of the array

  vertices[vertexCount++] = vertex;
  }

  void Face::removeVertex(Vertex *vertex)
  {
  assert(vertex!=0);

  // locate the vertex in the array and replace it with the current last vertex
  // if already the last vertex, just overwrite it
  // slow: should make this a doubly-linked list ???

  for (unsigned int i = vertexCount; i>0; i--)
    if (vertices[i-1]==vertex)
    {
      vertices[i-1] = vertices[--vertexCount];
      return;
    }

  assert(0);
  }
  void Face::setMu(){
      //setting the values of TargetFormMatrix first, to calculate the value of Mu
      //this->setTargetFormMatrix();
      unsigned int faceid = this->getID();//id of this face
      //iterating this face for the vertices
      FaceEdgeIterator faceIterator(this);
      Edge *currentEdge;//to keep the current edge
      Vertex *vertOrg;//to keep the current vertex
      double ak,f1,f2,f3;//to store the values of functions
      double term11(0), term12(0), term22(0); // to store the values of summing terms to calculate mu
      //iterating through all the edges
      while ((currentEdge = faceIterator.next())!=0){
            vertOrg = currentEdge->Org();
            ak = vertOrg->getAk(faceid);
            f1 = vertOrg->getFunction1(faceid);
            f2 = vertOrg->getFunction2(faceid);
            f3 = vertOrg->getFunction3(faceid);
            //calculating the sumterms
            term11 += (1./12.)*ak*f1;
            term12 += (1./24.)*ak*f2;
            term22 += (1./12.)*ak*f3;
      }
      //setting the current form matrix
      this->currentFormMatrix[0][0] = term11;
      this->currentFormMatrix[0][1] = term12;
      this->currentFormMatrix[1][0] = term12;
      this->currentFormMatrix[1][1] = term22;
      //setting mu
      this->mu1 = term11 - targetFormMatrix[0][0];
      this->mu2 = term12 - targetFormMatrix[1][0];
      this->mu3 = this->mu2;//as mu2 and mu3 are equal
      this->mu4 = term22 - targetFormMatrix[1][1];
  }
void Face::setTargetFormMatrix(){
    //*******************************************//
    // this should be edited to include time update
    // at each growth step
    //*******************************************//
     //iterating this face for the vertices
    /*
    unsigned int faceid = this->getID();//id of this face
    FaceEdgeIterator faceIterator(this);
    Edge *currentEdge;//to keep the current edge
    Vertex *vertOrg;//to keep the current vertex
    double ak,f1,f2,f3;//to store the values of functions
    double term11(0), term12(0), term22(0); // to store the values of matrix element
    //iterating through all the edges
    while ((currentEdge = faceIterator.next())!=0){
          vertOrg = currentEdge->Org();
          ak = vertOrg->getAk(faceid);
          f1 = vertOrg->getFunction1(faceid);
          f2 = vertOrg->getFunction2(faceid);
          f3 = vertOrg->getFunction3(faceid);
          //calculating the terms of the form matrix
          term11 += (1./12.)*ak*f1;
          term12 += (1./24.)*ak*f2;
          term22 += (1./12.)*ak*f3;
    }
    //putting the values in the formatrix
    this->targetFormMatrix[0][0] = term11;
    this->targetFormMatrix[1][0] = term12;
    this->targetFormMatrix[0][1] = term12;
    this->targetFormMatrix[1][1] = term22;
    this->setTraceSquaredTargetFormMatrix();
    */
    //just setting a circle target form matrix
  /*identity matrix
    this->targetFormMatrix[0][0] = 1.;
    this->targetFormMatrix[1][0] = 0.;
    this->targetFormMatrix[0][1] = 0.;
    this->targetFormMatrix[1][1] = 1.;
    */
    //setting this equal to the Mc
    this->targetFormMatrix[0][0] = (this->getMu1() + 0. );
    this->targetFormMatrix[1][0] = (this->getMu2() + 0. );
    this->targetFormMatrix[0][1] = (this->getMu3() + 0. );
    this->targetFormMatrix[1][1] = (this->getMu4() + 0. );
    this->setTraceSquaredTargetFormMatrix();
  }
//***************************************************************************** //
void Face::printTargetFormMatrix(){
  std::cout<< "face id : "<< this->getID() << std::endl;
  std::cout<< "TargetFormMatrix"<<std::endl;
  std::cout<< this->targetFormMatrix[0][0] << "   "<<this->targetFormMatrix[0][1]<<std::endl;
  std::cout<< this->targetFormMatrix[0][1] << "   "<<this->targetFormMatrix[1][1]<<std::endl;
  std::cout<< "Current Form Matrix" <<std::endl;
  std::cout<< this->currentFormMatrix[0][0] << "   "<<this->currentFormMatrix[0][1]<<std::endl;
  std::cout<< this->currentFormMatrix[0][1] << "   "<<this->currentFormMatrix[1][1]<<std::endl;
  std::cout<< "difference of targetFormMatrix and current form matrix" <<std::endl;
  std::cout<< this->targetFormMatrix[0][0]-currentFormMatrix[0][0] << "   "<<this->targetFormMatrix[0][1]-currentFormMatrix[0][1]<<std::endl;
  std::cout<< this->targetFormMatrix[0][1]-currentFormMatrix[0][1] << "   "<<this->targetFormMatrix[1][1]-currentFormMatrix[1][1]<<std::endl;
}

//***************************************************************************** //

void Face::setTempTargetFormMatrix(){
    this->targetFormMatrix[0][0] = (this->getMu1() + this->targetFormMatrix[0][0] );
    this->targetFormMatrix[1][0] = (this->getMu2() + this->targetFormMatrix[1][0] );
    this->targetFormMatrix[0][1] = (this->getMu3() + this->targetFormMatrix[0][1] );
    this->targetFormMatrix[1][1] = (this->getMu4() + this->targetFormMatrix[1][1] );
    this->setTraceSquaredTargetFormMatrix();
}
//***************************************************************************** //
 void Face::setTraceSquaredTargetFormMatrix(){
    double lambda1 = this->targetFormMatrix[0][0];
    double lambda2 = this->targetFormMatrix[1][1];
    double trace =  (lambda1+lambda2)*(lambda1+lambda2);
    this->traceSquaredTargetFormMatrix = trace;
 }
 //************ ENERGY VALUE SETTER ******************** //
 // First, Second, Third terms are set here ! //
 void Face::setEnergyTerms(){
  //getting need properties from this face
  if (this->getID() != 1){
      double area = this->getAreaOfFace();//area of face
      Cell *cell = this->getCell();//getting the cell on which this face lies
      double alpha = cell->getAlpha();
      double beta = cell->getBeta();
      double pressure = cell->getPressure();
      unsigned int faceid = this->getID();//getting the face id
      // ****************************************************************************************** //
      //setting the target form matrix of this face
      //this->setTargetFormMatrix();
      double (*targetFormMatrix)[2] = this->targetFormMatrix;
      //std::copy(&face->targetFormMatrix[0][0],&face->targetFormMatrix+2*2,&targetFormMatrix[0][0]);
      double targetFormMatrixTraceSquared = this->getTraceSquaredTargetFormMatrix();
      // ****************************************************************************************** //
      //temporary variables to calculate centralised Coordinates
      double firstterm(0), secondterm(0), thirdterm(0);// to store the values of each term of energy
      // ****************************************************************************************** //
      // calculating the Form matrix of this face
        FaceEdgeIterator faceIterator(this);
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
        firstterm = (1./targetFormMatrixTraceSquared)*area*(pow(currentFormMatrix[0][0]-targetFormMatrix[0][0],2)+
                  pow(currentFormMatrix[0][1]-targetFormMatrix[0][1],2)+
                  pow(currentFormMatrix[1][0]-targetFormMatrix[1][0],2)+
                  pow(currentFormMatrix[1][1]-targetFormMatrix[1][1],2));
        //Second Term
        secondterm = (1./targetFormMatrixTraceSquared)*area*(pow(currentFormMatrix[0][0]-targetFormMatrix[0][0]+
                    currentFormMatrix[1][1]-targetFormMatrix[1][1],2));
        //Third term
        thirdterm = area;
        // ****************************************************************************************** //
        //calculating energy
        double energytemp = alpha*firstterm + beta*secondterm - pressure*thirdterm;
        //setting the energy values of face 
        this->firstTerm = firstterm;
        this->secondTerm = secondterm;
        this->thirdTerm = thirdterm;
        this->energy = energytemp;
  }
 }
 // *************************************************************** //
 void Face::grow(){
  //currently using lockhardt model 
  //thresholdmatrix : property of cell
  Cell *thiscell = this->getCell();//getting the cell's threshold matrix
  double kappa = thiscell->getKappa();
  // performing the growth of targetFormMatrix
  this->targetFormMatrix[0][0] = this->targetFormMatrix[0][0] + (thiscell->hstepsize)*kappa*
                                    std::max(0.,(this->getMu1()-thiscell->thresholdMatrix[0][0]));
  this->targetFormMatrix[1][0] = this->targetFormMatrix[1][0] + (thiscell->hstepsize)*kappa*
                                    std::max(0.,(this->getMu3()-thiscell->thresholdMatrix[1][0]));
  this->targetFormMatrix[0][1] = this->targetFormMatrix[0][1] + (thiscell->hstepsize)*kappa*
                                    std::max(0.,(this->getMu2()-thiscell->thresholdMatrix[0][1]));
  this->targetFormMatrix[1][1] = this->targetFormMatrix[1][1] + (thiscell->hstepsize)*kappa*
                                    std::max(0.,(this->getMu4()-thiscell->thresholdMatrix[1][1]));

 }

  //****************** end added features********************************//
/* -- protected instance methods ------------------------------------------- */

Face::Face(Cell *cell)
{
  assert(cell!=0);

  this->cell = cell;
  this->id   = cell->makeFaceID();
  this->data = 0;
  this->edge = 0;
  //***************added features*******************************************//
  this->vertices = new Vertex*[8];
  this->vertexCount = 0;
  this->vertexSize = 8;
  //***************end added features*******************************************//
  cell->addFace(this);
}

Face::~Face()
{
  //***************added features***************************************//
  {
    for (unsigned int i = vertexCount; i>0; i--)
      Vertex::kill(vertices[i-1]);
  }
  delete[] vertices;
  //***************end added features*************************************//
  cell->removeFace(this);
}

