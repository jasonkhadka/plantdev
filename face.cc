
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


//random number generating
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>

//for eigenvalue conputation
#include "./eigen/Eigen/Dense"
#include "./eigen/Eigen/Eigenvalues"

//struct to store verticies
struct vertex_coordinate {
    unsigned int id;
    double x,y,z;
};
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
void Face::setVertexCount(){
  Edge * currentEdge;
  FaceEdgeIterator edges(this); //iterating through this face
  int counter = 0;
  while ((currentEdge = edges.next())!= 0){
    //iterating this face until exhausted
      counter++;//every time edge is switched counter is updated
  }
  this->vertexCount = counter;//putting the number of vertices in the vertexCount
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
void Face::setKappa(double newkappa){
  kappa = newkappa;
  this->lastGrowthRate = this->getGrowthRandomNumber();
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
  //First : getting the number of vertices in this face
  this->setVertexCount();
  int numOfVertex = this->getVertexCount();
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
    // getting the mean center position of this face //
    unsigned int faceid = this->getID();// first grabbing id of  current face
    //std::cout<<"face id : "<<faceid<< "numberOfvertex :"<<numOfVertex<<std::endl;
    double xCentroid(0), yCentroid(0), zCentroid(0); // coordinate of the cnetroid
    // array of vertices
    double xcood[numOfVertex], ycood[numOfVertex], zcood[numOfVertex];// coordinate of the vertices of this face
    double xTriCen[numOfVertex], yTriCen[numOfVertex], zTriCen[numOfVertex], triangulatedArea[numOfVertex];// coordinate of center of triangluated triangles of face
    //areaTri : area of triangles
    double xmean(0), ymean(0), zmean(0);//mean center of the face
    int counter(0);//counter to keep track of vertices
    Edge *currentEdge;//pointer to keep track of a edge that is iterated 
    Vertex *currentVertex; // vertex to store the dest() of currentEdge
    //getting the mean center of the face
    FaceEdgeIterator faceEdges(this);//iterator to iterate through the vertex for outgoign edge
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
    //std::cout<<"Numberofvertex : "<<numOfVertex << "|  counter : "<< counter <<"|  ";
    xmean = xmean/numOfVertex;
    ymean = ymean/numOfVertex;
    zmean = zmean/numOfVertex;
    //std::cout<<"xmean : "<< xmean<<" ymean : "<< ymean<<" zmean : "<< zmean<<std::endl;
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
    // triangulating the face and getting the center of each triangle
    //vectors of the triangles 
    double vector1[3],vector2[3];
    double totalarea = 0;
    double crossProductVector[3];//cross product of the vector
    double crossProductMagnitude;//cross product magintude
    for (int counter = 0; counter<numOfVertex; counter++){
          // calculating center of the triangles
          xTriCen[counter] = 1./3.*(xmean+xcood[counter]+xcood[(counter+1)%numOfVertex]);
          yTriCen[counter] = 1./3.*(ymean+ycood[counter]+ycood[(counter+1)%numOfVertex]);
          zTriCen[counter] = 1./3.*(zmean+zcood[counter]+zcood[(counter+1)%numOfVertex]);
          //printf("counter = %d; center of triangles : %F ; %F ; %F \n", counter, xTriCen[counter], yTriCen[counter], zTriCen[counter]);
          //getting two vectors of this triangle
          vector1[0] = xcood[counter]-xmean;
          vector1[1] = ycood[counter]-ymean;
          vector1[2] = zcood[counter]-zmean;
          vector2[0] = xcood[(counter+1)%numOfVertex]-xmean;
          vector2[1] = ycood[(counter+1)%numOfVertex]-ymean;
          vector2[2] = zcood[(counter+1)%numOfVertex]-zmean;
          //cross product of the two vectors
          crossProductVector[0] = vector1[1]*vector2[2] - vector1[2]*vector2[1];
          crossProductVector[1] = vector1[2]*vector2[0] - vector1[0]*vector2[2];
          crossProductVector[2] = vector1[0]*vector2[1] - vector1[1]*vector2[0];
          //maginitude of cross product
          crossProductMagnitude = sqrt(abs(crossProductVector[0]*crossProductVector[0] + crossProductVector[1]*crossProductVector[1]+
                                  crossProductVector[2]*crossProductVector[2]));
          triangulatedArea[counter] = 0.5*crossProductMagnitude;//area of this triangle
          xCentroid += xTriCen[counter]*triangulatedArea[counter];//adding the weigthed centroid
          yCentroid += yTriCen[counter]*triangulatedArea[counter];
          zCentroid += zTriCen[counter]*triangulatedArea[counter]; 
          //adding total area
          totalarea += triangulatedArea[counter];// calculating total area of triange
      }
    //std::cout<<"xCentroid: "<< xCentroid<<" yCentroid : "<< yCentroid<<" zCentroid : "<< zCentroid<<std::endl;
    //std::cout<<"faceID: "<<this->getID()<<"| totalarea : "<< totalarea<<std::endl;
    //now calculating weighted center
    xCentroid = xCentroid/totalarea;
    yCentroid = yCentroid/totalarea;
    zCentroid = zCentroid/totalarea;
    //std::cout<< "mean : "<<xmean<<"  " << ymean <<"  "<< zmean <<std::endl;
    //std::cout<< "centroid : "<<xCentroid<<"  " << yCentroid <<"  "<< zCentroid<<std::endl;
    // setting the central coordinate of this face in terms of cartisian coordinate
    this->xCentralised = xCentroid;
    this->yCentralised = yCentroid;
    this->zCentralised = zCentroid;
    if (this->getID() == 1){
      this->xCentralised = 0.;
      this->yCentralised = 0.;
      this->zCentralised = 0.;
    }
    //printf("%s %u %F \n ", "calculating total area of Face :",faceid, totalarea);
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
    //Calculating normal of the triangles of the shape and then calculating weighted normal
    //double normalTempX[numOfVertex],normalTempY[numOfVertex],normalTempZ[numOfVertex],totalarea;//vertices of normal of triangles
    double normalX(0), normalY(0), normalZ(0), areaTri(0); // vertices of weighted normal X, Y, Z 
    totalarea = 0;//resetting total area for calculation
    // iterating the faces also
    FaceEdgeIterator faceEdges1(this);//iterator to iterate through the vertex for outgoign edge
    // curentVertex and currentEdge have been already declared above
    currentEdge = faceEdges1.next();// grabbing the one edge of face = this
    for (int counter = 0; counter<numOfVertex; counter++){//triangulating the vertices with centroid of this face
          currentVertex = currentEdge->Dest();//getting the destination of the edge
          //getting two vectors of this triangle
          vector1[0] = xcood[counter]-xCentroid;
          vector1[1] = ycood[counter]-yCentroid;
          vector1[2] = zcood[counter]-zCentroid;
          vector2[0] = xcood[(counter+1)%numOfVertex]-xCentroid;
          vector2[1] = ycood[(counter+1)%numOfVertex]-yCentroid;             
          vector2[2] = zcood[(counter+1)%numOfVertex]-zCentroid;
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
    //std::cout<<"xCentroid: "<< xCentroid<<" yCentroid : "<< yCentroid<<" zCentroid : "<< zCentroid<<std::endl;
    //std::cout<<"totalarea : "<< totalarea<<std::endl;
    
    //printf("%s %u %F \n ", "calculating total area of Face :",faceid, totalarea);
    //weighted normal to this face
    normalX = normalX/totalarea;
    normalY = normalY/totalarea;
    normalZ = normalZ/totalarea;
    //setting area of the face
    // this totalarea is not the exact area of face, 
    //so instead calling setareaofface fucntion
    // however, this area is the total area for derivative calculation
    //this->areaOfFace = totalarea;
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
    // Now choosing Unit vector in intrinsic X - Direction
    // For that we need a perpendicular vector to Normal
    // this can be achieved by switching any two coordinate of normal and 
    // add a minus sign to one of them, and complete the vector with zeroes.
    if (pow(normalX*normalX + normalY*normalY,0.5) > pow(10.,-6)){//normal is not parrelel to z axis
        vectorVertex[0] = -1*normalY;
        vectorVertex[1] = normalX;
        vectorVertex[2] = 0;
    }else{//normal is parrelel to z axis
        vectorVertex[0] = normalZ;
        vectorVertex[1] = 0;
        vectorVertex[2] = -1*normalX;
    }
    //dot product of vectorVertex and normal
    dotproduct = vectorVertex[0]*normalX+vectorVertex[1]*normalY+vectorVertex[2]*normalZ;
    //now calculating the projected vertices --at this stage unitx is storing pi_vector
    unitx[0]  = vectorVertex[0] +xCentroid- dotproduct*normalX;
    unitx[1] = vectorVertex[1]+yCentroid - dotproduct*normalY;
    unitx[2] = vectorVertex[2]+zCentroid - dotproduct*normalZ;
    {//saving the pi_vector of this face
    double * pntPivector = unitx;
    this->setPivector(pntPivector);
    }
    ////getting the unitx = normalised[Projectedvertex1 - Origin]
   //    Failed one : -> //now getting the vertex on the plane to take it as the x-direction
    //unitx = Normalised[Origin + wvector]
    unitx[0] = unitx[0]- xCentroid;
    unitx[1] = unitx[1]- yCentroid;
    unitx[2] = unitx[2]- zCentroid;
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
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
    // now projecting the 3d coordinates to the 2D coordinates and assinging the vertices
    double xprojection, yprojection, zprojection;
    double gamma1, gamma2, gamma3;
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
          // NON CENTRALISED Projected Coordinate or coordinates in terms of CARTESIAN SYSTEM
          vectorVertex[0] = xprojection;
          vectorVertex[1] = yprojection;
          //vectorVertex[2] = zprojection; //taking zproj means fully changing 3D shape to cartesian system
          vectorVertex[2] = 0.; //taking zproj = 0. means changing only x-y projection to cartesian system
          //now getting the new x, y coordinates, dot product of unit vector with the vectorVertex
          //This converts fully the projected coordinates to cartesian
          gamma1 = unitx[0]*vectorVertex[0]+unity[0]*vectorVertex[1]+unitz[0]*vectorVertex[2] + xCentroid;
          gamma2 = unitx[1]*vectorVertex[0]+unity[1]*vectorVertex[1]+unitz[1]*vectorVertex[2] + yCentroid;
          gamma3 = unitx[2]*vectorVertex[0]+unity[2]*vectorVertex[1]+unitz[2]*vectorVertex[2] + zCentroid;
          //now setting the projected coordinates in the vertex properties
          currentVertex->insertNonCentralisedProjectedXcoordinate(faceid,gamma1);
          currentVertex->insertNonCentralisedProjectedYcoordinate(faceid,gamma2); 
          currentVertex->insertNonCentralisedProjectedZcoordinate(faceid,gamma3); 
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
    //double theta = acos(dotproduct);//since both unitx and x unit vector have 1 magnitude, just dividing by 1
    //this->setAngleOfTilt(theta);
    //setting areas for the face
    //this->setAreaOfFace();
    // setting Mu values 
    //setting Mu values 
    //this->setMu();
    // ***************************************************************************** //
    //Setting Trace of Target Form Matrix Squared
    //this->setTraceSquaredTargetFormMatrix();
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
bool Face::isConvex(){
    Edge * edge;
    Vertex * first, * second, * third;
    unsigned int initialVertID;
    double angleThreshold = (this->getCell())->getConvexAngleThreshold();
    unsigned int faceid = this->getID();
    /*
    FaceEdgeIterator edges(this);
    edge = edges.next();
    first = edge->Org();
    second = edge->Dest();
    double x1,y1,x2,y2,zCross; // variables declaration for z-component of cross product calcuation
    //double vec1[3],vec2[3];//variables declaration for angle calculation
    while ((edge = edges.next())!= 0){//NEED TO IMPLEMENT THE MECHANISM TO CALCULATE angle for all vertex
            // NEED IMPROVED LOOPING ! 
          third = edge->Dest();
          // Algorithm to calculate the z-component of cross-product between two vecs
          //Gathering the two vectors from first & second Edge
          x1 = second->getProjectedXcoordinate(this->getID())-first->getProjectedXcoordinate(this->getID());
          y1 = second->getProjectedYcoordinate(this->getID())-first->getProjectedYcoordinate(this->getID());
          x2 = third->getProjectedXcoordinate(this->getID())-second->getProjectedXcoordinate(this->getID());
          y2 = third->getProjectedYcoordinate(this->getID())-second->getProjectedYcoordinate(this->getID());
          //now calculating & checking the z-component of cross-product of this
          zCross = x1*y2-y1*x2;
          if (zCross<0){ return false;};
          
          //if above condition is not satisfied, this polygon can still be Convex, so continuing with search
          first = edge->Org();
          second = edge->Dest();
    }
    */
    //std::cout<<" ============================== \n" << faceid << "\n ============================== \n";

    edge = this->getEdge();
    first = edge->Org();
    second = edge->Dest();
    initialVertID = first->getID();
    //arrays to use for angle calculation
    double firstvec[2], secondvec[2];
    double dot, det, angle;
    while (true){
        edge = edge->Lnext();
        third = edge->Dest();
        // *****  Calculation of Angle  ***** //
        // getting the two vectors for calculations
        firstvec[0] = third->getProjectedXcoordinate(faceid) - second->getProjectedXcoordinate(faceid);
        firstvec[1] = third->getProjectedYcoordinate(faceid) - second->getProjectedYcoordinate(faceid);
        secondvec[0] = first->getProjectedXcoordinate(faceid) - second->getProjectedXcoordinate(faceid);
        secondvec[1] = first->getProjectedYcoordinate(faceid) - second->getProjectedYcoordinate(faceid);
        // dot product of two vectors //
        dot = firstvec[0]*secondvec[0] + firstvec[1]*secondvec[1];
        det = firstvec[0]*secondvec[1]-firstvec[1]*secondvec[0];
        // now calculating the angle [0,360]
        angle = fmod(((atan2(det,dot))/(M_PI)*180.+360.),360.);
        //std::cout << " Vert ID : "<< second->getID() << " Angle  : "<< angle << std::endl;
        if (angle > angleThreshold){//check if angle is greater than the threshold
          return false;
        }
        //std::cout<< " initialVertID : "<< initialVertID << "edge->Org()"<< (edge->Org())->getID()<<std::endl;
        // checking the exit condition
        if ((edge->Org())->getID() == initialVertID){
            break;
        }else
        {
          first = edge->Org();
          second = edge->Dest();
        }
    }
  //if all vertex has been checked for having angle >180 (here condition that z-cross product is positive)
  //then this is a Convex polygon
  return true;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //

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
double * Face::getStressEigenVector1(){
  double * pntunit = this->stressEigenVector1;
  return pntunit;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
double * Face::getStressEigenVector2(){
  double * pntunit = this->stressEigenVector2;
  return pntunit;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
double * Face::getStrainEigenVector1(){
  double * pntunit = this->strainEigenVector1;
  return pntunit;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
double * Face::getStrainEigenVector2(){
  double * pntunit = this->strainEigenVector2;
  return pntunit;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
double Face::getAreaOfFace(){
	return this->areaOfFace;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
double * Face::getNoncentralisedUnitx(){
  double * pntunit = this->nonCentralisedUnitx;
  return pntunit;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
double * Face::getNoncentralisedUnity(){
  double * pntunit = this->nonCentralisedUnity;
  return pntunit;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
double * Face::getNoncentralisedUnitz(){
  double * pntunit = this->nonCentralisedUnitz;
  return pntunit;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
void Face::setNoncentralisedunitx(double * tempunit){
  for (int i = 0; i<3; i++){
      nonCentralisedUnitx[i] = tempunit[i];
   }
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
void Face::setNoncentralisedunity(double * tempunit){
  for (int i = 0; i<3; i++){
      nonCentralisedUnity[i] = tempunit[i];
   }
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
void Face::setNoncentralisedunitz(double * tempunit){
  for (int i = 0; i<3; i++){
      nonCentralisedUnitz[i] = tempunit[i];
   }
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
void Face::setAreaOfFace(){
  /*
      * Polygon calculation of Area of Face
  */
  //double vertexnumber =(double) this->countVertices();
  double total[3] ={0.,0.,0.};
  double cross[3] = {0.,0.,0.};
  // Iterating through the edges of this face
  FaceEdgeIterator edges(this);//iterator to iterate through the face
  Edge *newedge;//a pointer to keep track of current edge
  Vertex *first, *second;
  //looping through all the edges in the face until exhausted
  while ((newedge = edges.next())!=0){
      first = newedge->Org();
      second = newedge->Dest();
      //cross product of two vertex vectors
      cross[0] = first->getYcoordinate()*second->getZcoordinate()-first->getZcoordinate()*second->getYcoordinate();
      cross[1] = first->getZcoordinate()*second->getXcoordinate()-first->getXcoordinate()*second->getZcoordinate();
      cross[2] = first->getXcoordinate()*second->getYcoordinate()-first->getYcoordinate()*second->getXcoordinate();
      // adding to total vector
      total[0]+= cross[0];
      total[1]+= cross[1];
      total[2]+= cross[2];   
    }
  // now getting the unit normal to the face
  double * unitnormal = this->getNormal();
  double result = unitnormal[0]*total[0]+unitnormal[1]*total[1]+unitnormal[2]*total[2];
this->areaOfFace = abs(result)/2;


  /*
  if (this->getID() == 1){
    return;
  }
	FaceEdgeIterator edges(this);//iterator to iterate through the face
	Edge *newedge;//a pointer to keep track of current edge
	double areasum = 0;//variable to store the area of the face
	unsigned int faceid = this->id;//getting the faceid
  Vertex *vertexOrg;
	//looping through all the edges in the face until exhausted
	while ((newedge = edges.next())!=0){
    vertexOrg = newedge->Org();
    vertexOrg->setparameters();//setting parameters required to calculate area
		areasum += vertexOrg->getAk(faceid);//summing up all the Ak values for this face
	}
	this->areaOfFace = 0.5*areasum;//storing the area of the face in areaOfFace variable
  //in case the area is calculated negative (for example the crossing of vertices)
  //then assing max value to area
  if (this->areaOfFace <0.){
      this->areaOfFace = HUGE_VAL;
  }
  */
}
//******************added features********************************//
/*
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
*/

  void Face::setMu(){
    if (this->getID() == 1){
    return;
  }
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
            //calculating the sumterms // switched term11 and term22's definition as now 
            // term11 takes in to account area distribution in x - direction 
            // term22 takes in to account area distribution in y - direction
            term22 += (1./12.)*ak*f1;
            term12 += (1./24.)*ak*f2;
            term11 += (1./12.)*ak*f3;
      }
      //setting the current form matrix
      this->currentFormMatrix[0][0] = term11;
      this->currentFormMatrix[0][1] = term12;
      this->currentFormMatrix[1][0] = term12;
      this->currentFormMatrix[1][1] = term22;
      //setting mu
      this->mu1 = term11 - targetFormMatrix[0][0];
      this->mu2 = term12 - targetFormMatrix[0][1];
      this->mu3 = term12 - targetFormMatrix[1][0];//this->mu2;//as mu2 and mu3 are equal
      this->mu4 = term22 - targetFormMatrix[1][1];
      //Also calculating the area of this face
      this->setAreaOfFace();
  }
void Face::setTargetFormMatrix(){
  if (this->getID() == 1){
    return;
  }
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
  // This is master // 
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

void Face::setInitialTargetFormMatrixCurrent(){
    double facestrain = (this->getCell())->getInitialStrain();//getting the Initial strain to be used
    this->targetFormMatrix[0][0] = ((1.-facestrain)*this->getMu1());
    this->targetFormMatrix[1][0] = ((1.-facestrain)*this->getMu2());
    this->targetFormMatrix[0][1] = ((1.-facestrain)*this->getMu3());
    this->targetFormMatrix[1][1] = ((1.-facestrain)*this->getMu4());
    this->setTraceSquaredTargetFormMatrix();
    // Setting the intrinsic constant growth rate for FeedbackConstantGrow
    // currently, constantgrowthmatrix = 1% of TFM
    this->constantGrowthMatrix[0][0] = 0.01*this->targetFormMatrix[0][0];
    this->constantGrowthMatrix[0][1] = 0.01*this->targetFormMatrix[0][1];
    this->constantGrowthMatrix[1][0] = 0.01*this->targetFormMatrix[1][0];
    this->constantGrowthMatrix[1][1] = 0.01*this->targetFormMatrix[1][1];
    this->setTargetArea(this->getAreaOfFace());//area of current form matrix as target area
    this->lastGrowthRate = this->getGrowthRandomNumber();
}
//***************************************************************************** //

void Face::setSkewedTargetFormMatrixCurrent(){
  if (this->domePosition == true ){// IF this is dome face
    this->setInitialTargetFormMatrixCurrent();
  }
  else if (this->getMu1() > this->getMu4()){//Ixx > Iyy
      this->targetFormMatrix[0][0] = (this->getMu1());
      this->targetFormMatrix[1][0] = (this->getMu2());
      this->targetFormMatrix[0][1] = (this->getMu3());
      this->targetFormMatrix[1][1] = (0.5*this->getMu4());
  }else{//Iyy>Ixx
    this->targetFormMatrix[0][0] = (0.5*this->getMu1());
    this->targetFormMatrix[1][0] = (this->getMu2());
      this->targetFormMatrix[0][1] = (this->getMu3());
      this->targetFormMatrix[1][1] = (this->getMu4());
    }
    this->setTraceSquaredTargetFormMatrix();
    this->setTargetArea(this->getAreaOfFace());//area of current form matrix as target area
}


//***************************************************************************** //

void Face::setTempTargetFormMatrixIdentity(){
    this->targetFormMatrix[0][0] = 1.;
    this->targetFormMatrix[1][0] = 0.;
    this->targetFormMatrix[0][1] = 0.;
    this->targetFormMatrix[1][1] = 1.;
    this->setTraceSquaredTargetFormMatrix();
    this->setTargetArea(1.);//area of unit circle
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
      double cellalpha;
      if (this->alpha == 0){//means not update directly to face
          cellalpha = cell->getAlpha();
      }
      else{
          cellalpha = this->alpha;
      }
      double beta = cell->getBeta();
      double pressure = cell->getPressure();
      //unsigned int faceid = this->getID();//getting the face id
      double firstterm, secondterm, thirdterm;
      // ****************************************************************************************** //
      //std::copy(&face->targetFormMatrix[0][0],&face->targetFormMatrix+2*2,&targetFormMatrix[0][0]);
      double targetFormMatrixTraceSquared = this->getTraceSquaredTargetFormMatrix();
      // ****************************************************************************************** //
      // ****************************************************************************************** //
      //First Term
      firstterm = (1./targetFormMatrixTraceSquared)*area*(
                pow(this->getMu1(),2)+
                pow(this->getMu2(),2)+
                pow(this->getMu3(),2)+
                pow(this->getMu4(),2)
                );
      //Second Term
      secondterm = (1./targetFormMatrixTraceSquared)*area*(
                pow(this->getMu1()+this->getMu4(),2));
      //Third term
      thirdterm = area;
      // ****************************************************************************************** //
      //calculating energy
      double energytemp = cellalpha*firstterm + beta*secondterm - pressure*thirdterm;
      //removing the area term
      //double energytemp = cellalpha*firstterm + beta*secondterm;
      //setting the energy values of face 
      this->firstTerm = firstterm;
      this->secondTerm = secondterm;
      this->thirdTerm = thirdterm;
      this->energy = energytemp;
  }
 }
 // *************************************************************** //
 // Calculate and set the sum of Edge Length of this face
 void Face::setSumEdgeLength(){
    double sumlength = 0.;
    Vertex *first, *second;
    double x1,y1,z1,x2,y2,z2;
    FaceEdgeIterator edges(this);
    // Summing over all the edges
    while ((edge=edges.next())!=0){
          first = edge->Org();
          second = edge->Dest();
          // getting the coordinates of the vertices
          // first vertex
          x1 = first->getXcoordinate();
          y1 = first->getYcoordinate();
          z1 = first->getZcoordinate();
          // second vertex
          x2 = second->getXcoordinate();
          y2 = second->getYcoordinate();
          z2 = second->getZcoordinate();
          // calculate and sum the length of edge
          sumlength += pow(pow(x2-x1,2)+pow(y2-y1,2)+pow(z2-z1,2),2);
        }
    // setting the sum Edge Length 
    this->sumEdgeLength = sumlength;
 }
 // *************************************************************** //
 void Face::calculateStress(){
    if (this->getID() == 1){
    this->stress << 0.,0.,
                    0.,0.;
    return;
    }
    /*
    //Declarations
    Vertex *first, *second;
    //double to store the values
    double sxx(0.), syy(0.), sxy(0.);
    // calculating normal force on the face 
    double vertexNum = 0;
    //coordinates
    double x1,y1,z1, x2,y2,z2, length,normalforce(0.);
    //unsigned int initialVertID;
    Edge * edge;
    //Vectors of cartesianForce calculated with respect to cartesian coordinates
    Eigen::Vector3d forces1;
    Eigen::Vector3d forces2;
    Eigen::Vector3d intrinsicForces1;
    Eigen::Vector3d intrinsicForces2;
    //Transformation Matrix
    Eigen::Matrix3d transformationMatrix;
    transformationMatrix << this->unitx[0] , this->unitx[1], this->unitx[2],
                            this->unity[0] , this->unity[1], this->unity[2],
                            this->unitz[0] , this->unitz[1], this->unitz[2];
    //******Start of Calculations ******** //
    FaceEdgeIterator edges(this);
    while ((edge=edges.next())!=0){
          first = edge->Org();
          second = edge->Dest();
          //grabbing the coordinates
          x1 = first->getProjectedXcoordinate(this->getID());
          y1 = first->getProjectedYcoordinate(this->getID());
          z1 = first->getProjectedZcoordinate(this->getID());
          x2 = second->getProjectedXcoordinate(this->getID());
          y2 = second->getProjectedYcoordinate(this->getID());
          z2 = second->getProjectedZcoordinate(this->getID());
          //calculating the length of the side
          length = pow(pow((x2-x1),2)+pow((y2-y1),2),0.5);
          //getting the cartesianForce
          forces1 << first->cartesianForce[0],first->cartesianForce[1],first->cartesianForce[2];
          forces2 << second->cartesianForce[0],second->cartesianForce[1],second->cartesianForce[2];
          //now converting this too intrinsic coordinates
          intrinsicForces1 = transformationMatrix*forces1;
          intrinsicForces2 = transformationMatrix*forces2;
          // now calculating the terms 
          sxx += length*((intrinsicForces1[0]*x1 + intrinsicForces2[0]*x2)/3. + 
                 (intrinsicForces1[0]*x2 + intrinsicForces2[0]*x1)/6.);
          syy += length*((intrinsicForces1[1]*y1 + intrinsicForces2[1]*y2)/3. + 
                 (intrinsicForces1[1]*y2 + intrinsicForces2[1]*y1)/6.);
          sxy += length*((intrinsicForces1[0]*y1+intrinsicForces1[1]*x1+intrinsicForces2[0]*y2+intrinsicForces2[1]*x2)/6. +
                 (intrinsicForces1[0]*y2+intrinsicForces1[1]*x2+intrinsicForces2[0]*y1+intrinsicForces2[1]*x1)/12.);
          normalforce += intrinsicForces1[2];//sum normal forces to calculate average force
          vertexNum += 1.;
      }
      sxx /= this->getAreaOfFace();
      syy /= this->getAreaOfFace();
      sxy /= this->getAreaOfFace();
    //Now setting the stress Matrix on the 
      this->stress << sxx, sxy,
                      sxy, syy;
    */
this->calculateStrain();
Cell * cell = this->getCell();
double eta = cell->getEta();
double cellalpha;
if (this->alpha == 0){//means not update directly to face
          cellalpha = cell->getAlpha();
      }
      else{
          cellalpha = this->alpha;
      }
//Vectors of cartesianForce calculated with respect to cartesian coordinates
    Eigen::Vector3d forces1;
    Eigen::Vector3d forces2;
//Transformation Matrix
    Eigen::Matrix3d transformationMatrix;
    transformationMatrix << this->unitx[0] , this->unitx[1], this->unitx[2],
                            this->unity[0] , this->unity[1], this->unity[2],
                            this->unitz[0] , this->unitz[1], this->unitz[2];
Eigen::Matrix2d stressMatrix = cellalpha*(this->getAreaOfFace())*(this->strain);
// Saving Stress Matrix
this->stress<<stressMatrix(0,0),stressMatrix(0,1),
              stressMatrix(1,0),stressMatrix(1,1);
    //std::cout<< "====================== Calculation of Stress : facid :"<<this->getID()<<"====================== \n"<<
    //sxx <<"   "<< sxy <<"\n"
    //sxy <<"   "<< syy<<"\n"<<std::endl;

    // Now computing the Eigen values of the stress Matrix
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eigensolver;
    eigensolver.compute(this->stress);//computing the eigenvalues of stress
    // Converting EigenVector in intrinsic coordinate to the back to cartesian
    forces1 << eigensolver.eigenvectors().col(0)[0], eigensolver.eigenvectors().col(0)[1], 0.;
    forces1 = (transformationMatrix.transpose())*forces1;
    forces2 << eigensolver.eigenvectors().col(1)[0], eigensolver.eigenvectors().col(1)[1], 0.;
    forces2 = (transformationMatrix.transpose())*forces2;
    //saving the eigen values 
    this->stressEigenValue1 = eigensolver.eigenvalues()[0];
    this->stressEigenValue2 = eigensolver.eigenvalues()[1];
    this->stressEigenVector1[0] = forces1[0];
    this->stressEigenVector1[1] = forces1[1];
    this->stressEigenVector1[2] = forces1[2];
    this->stressEigenVector2[0] = forces2[0];
    this->stressEigenVector2[1] = forces2[1];
    this->stressEigenVector2[2] = forces2[2];
    //returning//
// Before calculation of Growth : Calculate the stress-strain of this cell
/* --------------------------------------------------------------------------------------
// CHANGE : This has been edited to calculate the stress matrix from Strain
// Stess = alpha*Area*Strain
   as, Stiffness = alpha*Area
   --------------------------------------------------------------------------------------
*/

    return;
 }
 // *************************************************************** //
 void Face::calculateStrain(){
  if (this->getID() == 1){
    this->strain << 0.,0.,
                    0.,0.;
    return;
    }
     //Vectors of cartesianForce calculated with respect to cartesian coordinates
    Eigen::Vector3d forces1;
    Eigen::Vector3d forces2;
    //Transformation Matrix
    Eigen::Matrix3d transformationMatrix;
    transformationMatrix << this->unitx[0] , this->unitx[1], this->unitx[2],
                            this->unity[0] , this->unity[1], this->unity[2],
                            this->unitz[0] , this->unitz[1], this->unitz[2];
  double traceofTargetForm = (this->targetFormMatrix[0][0]+ this->targetFormMatrix[1][1]);
  // Strain matrix = (M-M0)/trace(M0)
  this->strain<< 1./traceofTargetForm*(this->getMu1()),  1./traceofTargetForm*(this->getMu2()),
            1./traceofTargetForm*(this->getMu3()), 1./traceofTargetForm*(this->getMu4()); 
  // saving the determinant 
  this->strainDeterminant = strain.determinant();
  this->strainTrace = strain.trace();
  // Now computing the Eigen values of the Strain Matrix
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eigensolver;
    eigensolver.compute(this->strain);//computing the eigenvalues of strain
    // Converting EigenVector in intrinsic coordinate to the back to cartesian
    forces1 << eigensolver.eigenvectors().col(0)[0], eigensolver.eigenvectors().col(0)[1], 0.;
    forces1 = (transformationMatrix.transpose())*forces1;
    forces2 << eigensolver.eigenvectors().col(1)[0],eigensolver.eigenvectors().col(1)[1], 0.;
    forces2 = (transformationMatrix.transpose())*forces2;
    //saving the eigen values 
    this->strainEigenValue1 = eigensolver.eigenvalues()[0];
    this->strainEigenValue2 = eigensolver.eigenvalues()[1];
    this->strainEigenVector1[0] = forces1[0];
    this->strainEigenVector1[1] = forces1[1];
    this->strainEigenVector1[2] = forces1[2];
    this->strainEigenVector2[0] = forces2[0];
    this->strainEigenVector2[1] = forces2[1];
    this->strainEigenVector2[2] = forces2[2];

  return;
 }
 // *************************************************************** //
void Face::feedbackStrainGrow(){
  if(this->getID()== 1){
      return;
    }
// Before calculation of Growth : Calculate the stress-strain of this cell
/* --------------------------------------------------------------------------------------
// CHANGE : This has been edited to calculate the stress matrix from Strain
// Stess = alpha*Area*Strain
   as, Stiffness = alpha*Area
   --------------------------------------------------------------------------------------
*/

this->calculateStrain();
this->calculateStress();
Cell * cell = this->getCell();
double eta = cell->getEta();
double cellalpha;
if (this->alpha == 0){//means not update directly to face
          cellalpha = cell->getAlpha();
      }
      else{
          cellalpha = this->alpha;
      }
Eigen::Matrix2d stressMatrix = cellalpha*(this->getAreaOfFace())*(this->strain);
//getting traceless deviatoric matrix
Eigen::Matrix2d deviatoric = (stressMatrix) - 0.5*(stressMatrix.trace())*Eigen::Matrix2d::Identity();
//current Form Matrix
Eigen::Matrix2d M0;
M0 << this->targetFormMatrix[0][0],this->targetFormMatrix[0][1],
      this->targetFormMatrix[1][0],this->targetFormMatrix[1][1];
//get the feedback matrix
Eigen::Matrix2d feedback = deviatoric*M0 + M0*deviatoric;
//printing Feed back matrix
/*
std::cout<<"-------------------------face id "<<this->getID()<< "---------------------------"<<
"\n identity*trace"<<
"\n"<< 0.5*((this->stress).trace())*Eigen::Matrix2d::Identity()(0,0)<<"  "<< 0.5*((this->stress).trace())*Eigen::Matrix2d::Identity()(0,1)<<
"\n"<< 0.5*((this->stress).trace())*Eigen::Matrix2d::Identity()(1,0)<<"  "<< 0.5*((this->stress).trace())*Eigen::Matrix2d::Identity()(1,1)<<
"\n =========================================="<<
"\n Stress"<<
"\n"<<this->stress(0,0)<<"  "<<this->stress(0,1)<<
"\n"<<this->stress(1,0)<<"  "<<this->stress(1,1)<<
"\n =========================================="<<
"\n Strain"<<
"\n"<<this->strain(0,0)<<"  "<<this->strain(0,1)<<
"\n"<<this->strain(1,0)<<"  "<<this->strain(1,1)<<
"\n =========================================="<<
"\n TargetFormMatrix"<<
"\n"<<M0(0,0)<<"  "<<M0(0,1)<<
"\n"<<M0(1,0)<<"  "<<M0(1,1)<<
"\n =========================================="<<
"\n Deviatoric"<<
"\n"<<deviatoric(0,0)<<"  "<<deviatoric(0,1)<<
"\n"<<deviatoric(1,0)<<"  "<<deviatoric(1,1)<<
"\n =========================================="<<
"\n feedback matrix "<<
"\n"<<feedback(0,0)<<"  "<<feedback(0,1)<<
"\n"<<feedback(1,0)<<"  "<<feedback(1,1)<<std::endl;
*/
//growth rate of faces : randomized number between (kappa-0.5 to kappa + 0.5)
double growthfactor = this->getGrowthRandomNumber();
lastGrowthRate = growthfactor; //saving growth rate for plotting
//std::cout<<"Kappa : "<< kappa <<"/n Actual Growth Var  : "<<growthvar <<std::endl;
Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eigensolver;
// Growth Matrix
Eigen::Matrix2d growthMatrix;
// dM0/dt = kappa*STRAIN - n/2*feedback
growthMatrix<< growthfactor*(this->strain(0,0)) - eta/2.*feedback(0,0),
               growthfactor*(this->strain(0,1)) - eta/2.*feedback(0,1),
               growthfactor*(this->strain(1,0)) - eta/2.*feedback(1,0),
               growthfactor*(this->strain(1,1)) - eta/2.*feedback(1,1); 
eigensolver.compute(growthMatrix);//computing the eigenvalues of growthMatrix, to make sure it is always growing
//to calculate the individual growth eigen direction
  Eigen::Matrix2d eigen1;
  Eigen::Matrix2d eigen2;
  eigen1 = std::max(eigensolver.eigenvalues()[0]-cell->thresholdMatrix[0][0],0.0)*
                      ((eigensolver.eigenvectors().col(0))*(eigensolver.eigenvectors().col(0)).transpose());
  eigen2 = std::max(eigensolver.eigenvalues()[1]-cell->thresholdMatrix[1][1],0.0)*
                      ((eigensolver.eigenvectors().col(1))*(eigensolver.eigenvectors().col(1)).transpose());
// Now combining growth in both direction
  growthMatrix = eigen1 + eigen2;
//now setting the new targetFormMatrix with feedback matrix
this->targetFormMatrix[0][0] += growthMatrix(0,0);
this->targetFormMatrix[1][0] += growthMatrix(1,0);
this->targetFormMatrix[0][1] += growthMatrix(0,1);
this->targetFormMatrix[1][1] += growthMatrix(1,1);
/*
std::cout<<
"\n =========================================="<<
"\n Growth matrix "<<
"\n"<<growthMatrix(0,0)<<"  "<<growthMatrix(0,1)<<
"\n"<<growthMatrix(1,0)<<"  "<<growthMatrix(1,1)<<
"\n =============================================================================================="<<
std::endl;
*/
//now setting tracesq of Target Form Matrix
this->setTraceSquaredTargetFormMatrix();
}

 // *************************************************************** //
 void Face::grow(){
  if (this->getID() == 1){
    return;
  }
  //cell for this Face
  Cell * cell = this->getCell();
  //getting the strain matrix for this face
  Eigen::Matrix2d strain;
  //assigning values to strain values
  double traceofTargetForm = (this->targetFormMatrix[0][0]+ this->targetFormMatrix[1][1]);
  // Strain matrix = Mu-Matrix
  strain<< 1./traceofTargetForm*(this->getMu1()),  1./traceofTargetForm*(this->getMu2()),
            1./traceofTargetForm*(this->getMu3()), 1./traceofTargetForm*(this->getMu4()); 
  //Eigensolver for strain
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eigensolver;
  eigensolver.compute(strain);//computing the eigenvalues of strain
//growth rate of faces : randomized number between (kappa-0.5 to kappa + 0.5)
double growthfactor = this->getGrowthRandomNumber();
lastGrowthRate = growthfactor; //saving growth rate for plotting
//calculating the time derivative now

  /*
  std::cout<< "***************************************************************************"<<std::endl;
  std::cout<<"face id : "<<this->getID()<<std::endl;
  std::cout<< "strain matrix : "<<std::endl<<strain<<std::endl;
  std::cout<<"trace of target Form Matrix : "<<traceofTargetForm << std::endl;
  std::cout<< "Determinant of strain Matrix : "<<strain(0,0)*strain(1,1)-strain(0,1)*strain(1,0)<<std::endl;
  std::cout<< "***************************************************************************"<<std::endl;
  std::cout<<"Growth fluctiuation"<<fluctuation<<std::endl;
  std::cout<<"growthrate (kappa) : "<<kappa << " | after fluctuation (k*growthfluc+other) : "<< kappa*(1+(2*growthvar*fluctuation-growthvar))<<std::endl;
  */
  //for the growthrate (timederivative) calcualtion -> New_Mo = Old_mo + growthRate
  Eigen::Matrix2d growthRate;
  //to calculate the individual eigen direction of strain 
  Eigen::Matrix2d eigen1;
  Eigen::Matrix2d eigen2;
  //printing eigen values & vector//
  /*
  std::cout<<"The eigenvalues are "<<std::endl<<eigensolver.eigenvalues()<<std::endl;
  std::cout<<"The vectors are "<<std::endl<<eigensolver.eigenvectors()<<std::endl;
  std::cout<<"The thresholdMatrix "<<std::endl;
  std::cout<<cell->thresholdMatrix[0][0]<<"    "<<cell->thresholdMatrix[0][1]<<std::endl;
  std::cout<<cell->thresholdMatrix[1][0]<<"    "<<cell->thresholdMatrix[1][1]<<std::endl;
  */
  //calcuating the time derivative
  eigen1 = std::max(eigensolver.eigenvalues()[0]-cell->thresholdMatrix[0][0],0.0)*
                      ((eigensolver.eigenvectors().col(0))*(eigensolver.eigenvectors().col(0)).transpose());
  eigen2 = std::max(eigensolver.eigenvalues()[1]-cell->thresholdMatrix[1][1],0.0)*
                      ((eigensolver.eigenvectors().col(1))*(eigensolver.eigenvectors().col(1)).transpose());
  //calculating the time derivative now
  growthRate = growthfactor*(eigen1+eigen2);
  /*
  std::cout<<"Growth Rate addition to TargetFormMatrix : (from Face::grow() : Faceid :"<<this->getID()<<std::endl;
  std::cout<<growthRate<<std::endl;
  */
  //now setting the new targetFormMatrix
  this->targetFormMatrix[0][0] += growthRate(0,0);
  this->targetFormMatrix[1][0] += growthRate(1,0);
  this->targetFormMatrix[0][1] += growthRate(0,1);
  this->targetFormMatrix[1][1] += growthRate(1,1);
  //now setting tracesq of Target Form Matrix
  this->setTraceSquaredTargetFormMatrix();
 }

// *************************************************************** //
void Face::feedbackInflatedGrow(){
  if(this->getID()== 1){
      return;
    }
// Before calculation of Growth : Calculate the stress-strain of this cell
  this->calculateStress();
  this->calculateStrain();
Cell * cell = this->getCell();
double eta = cell->getEta();
//getting traceless deviatoric matrix
Eigen::Matrix2d deviatoric = (this->stress) - 0.5*((this->stress).trace())*Eigen::Matrix2d::Identity();
//current Form Matrix
Eigen::Matrix2d M0;
M0 << this->targetFormMatrix[0][0],this->targetFormMatrix[0][1],
      this->targetFormMatrix[1][0],this->targetFormMatrix[1][1];
//get the feedback matrix
Eigen::Matrix2d feedback = deviatoric*M0 + M0*deviatoric;
//printing Feed back matrix
/*
std::cout<<"-------------------------face id "<<this->getID()<< "---------------------------"<<
"\n identity*trace"<<
"\n"<< 0.5*((this->stress).trace())*Eigen::Matrix2d::Identity()(0,0)<<"  "<< 0.5*((this->stress).trace())*Eigen::Matrix2d::Identity()(0,1)<<
"\n"<< 0.5*((this->stress).trace())*Eigen::Matrix2d::Identity()(1,0)<<"  "<< 0.5*((this->stress).trace())*Eigen::Matrix2d::Identity()(1,1)<<
"\n =========================================="<<
"\n Stress"<<
"\n"<<this->stress(0,0)<<"  "<<this->stress(0,1)<<
"\n"<<this->stress(1,0)<<"  "<<this->stress(1,1)<<
"\n =========================================="<<
"\n TargetFormMatrix"<<
"\n"<<M0(0,0)<<"  "<<M0(0,1)<<
"\n"<<M0(1,0)<<"  "<<M0(1,1)<<
"\n =========================================="<<
"\n Deviatoric"<<
"\n"<<deviatoric(0,0)<<"  "<<deviatoric(0,1)<<
"\n"<<deviatoric(1,0)<<"  "<<deviatoric(1,1)<<
"\n =========================================="<<
"\n feedback matrix "<<
"\n"<<feedback(0,0)<<"  "<<feedback(0,1)<<
"\n"<<feedback(1,0)<<"  "<<feedback(1,1)<<std::endl;
*/
//growth rate of faces : randomized number between (kappa-0.5 to kappa + 0.5)
double growthfactor = this->getGrowthRandomNumber();
lastGrowthRate = growthfactor; //saving growth rate for plotting
/*
std::cout<<"Kappa :: "<< kappa <<"\n Actual Growth Var  : "<<growthvar <<std::endl;
std::cout<<"Feedback parameter : Eta :: " << eta <<std::endl;
*/
Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eigensolver;
// Growth Matrix
Eigen::Matrix2d growthMatrix;
// dM0/dt = kappa*M0 - n/2*feedback
growthMatrix<< growthfactor*(this->targetFormMatrix[0][0]) - eta/2.*feedback(0,0),
               growthfactor*(this->targetFormMatrix[0][1]) - eta/2.*feedback(0,1),
               growthfactor*(this->targetFormMatrix[1][0]) - eta/2.*feedback(1,0),
               growthfactor*(this->targetFormMatrix[1][1]) - eta/2.*feedback(1,1); 
eigensolver.compute(growthMatrix);//computing the eigenvalues of growthMatrix, to make sure it is always growing
//to calculate the individual growth eigen direction
  Eigen::Matrix2d eigen1;
  Eigen::Matrix2d eigen2;
  eigen1 = std::max(eigensolver.eigenvalues()[0]-cell->thresholdMatrix[0][0],0.0)*
                      ((eigensolver.eigenvectors().col(0))*(eigensolver.eigenvectors().col(0)).transpose());
  eigen2 = std::max(eigensolver.eigenvalues()[1]-cell->thresholdMatrix[1][1],0.0)*
                      ((eigensolver.eigenvectors().col(1))*(eigensolver.eigenvectors().col(1)).transpose());
// Now combining growth in both direction
  growthMatrix = eigen1 + eigen2;
//now setting the new targetFormMatrix with feedback matrix
this->targetFormMatrix[0][0] += growthMatrix(0,0);
this->targetFormMatrix[1][0] += growthMatrix(1,0);
this->targetFormMatrix[0][1] += growthMatrix(0,1);
this->targetFormMatrix[1][1] += growthMatrix(1,1);
/*
std::cout<<
"\n =========================================="<<
"\n Face ID : : "<< this->getID() <<
"\n eta : : "<< eta <<
"\n =========================================="<<
"\n Growth matrix "<<
"\n"<<growthMatrix(0,0)<<"  "<<growthMatrix(0,1)<<
"\n"<<growthMatrix(1,0)<<"  "<<growthMatrix(1,1)<<
"\n =========================================="<<
"\n New Target Form matrix "<<
"\n"<<this->targetFormMatrix[0][0]<<"  "<<this->targetFormMatrix[0][1]<<
"\n"<< this->targetFormMatrix[1][0]<<"  "<< this->targetFormMatrix[1][1]<<std::endl;
*/
//now setting tracesq of Target Form Matrix
this->setTraceSquaredTargetFormMatrix();
}

// *************************************************************** //

void Face::feedbackConstantGrow(){
  if(this->getID()== 1){
      return;
    }
// Before calculation of Growth : Calculate the stress-strain of this cell
  this->calculateStrain();
  this->calculateStress();
Cell * cell = this->getCell();
double eta = cell->getEta();
//getting traceless deviatoric matrix
Eigen::Matrix2d deviatoric = (this->stress) - 0.5*((this->stress).trace())*Eigen::Matrix2d::Identity();
// Constant Growth Matrix
Eigen::Matrix2d MG;
MG << this->constantGrowthMatrix[0][0],this->constantGrowthMatrix[0][1],
      this->constantGrowthMatrix[1][0],this->constantGrowthMatrix[1][1];
//get the feedback matrix
Eigen::Matrix2d feedback = deviatoric*MG + MG*deviatoric;
//printing Feed back matrix
/*
std::cout<<"-------------------------face id "<<this->getID()<< "---------------------------"<<
"\n identity*trace"<<
"\n"<< 0.5*((this->stress).trace())*Eigen::Matrix2d::Identity()(0,0)<<"  "<< 0.5*((this->stress).trace())*Eigen::Matrix2d::Identity()(0,1)<<
"\n"<< 0.5*((this->stress).trace())*Eigen::Matrix2d::Identity()(1,0)<<"  "<< 0.5*((this->stress).trace())*Eigen::Matrix2d::Identity()(1,1)<<
"\n =========================================="<<
"\n Stress"<<
"\n"<<this->stress(0,0)<<"  "<<this->stress(0,1)<<
"\n"<<this->stress(1,0)<<"  "<<this->stress(1,1)<<
"\n =========================================="<<
"\n TargetFormMatrix"<<
"\n"<<M0(0,0)<<"  "<<M0(0,1)<<
"\n"<<M0(1,0)<<"  "<<M0(1,1)<<
"\n =========================================="<<
"\n Deviatoric"<<
"\n"<<deviatoric(0,0)<<"  "<<deviatoric(0,1)<<
"\n"<<deviatoric(1,0)<<"  "<<deviatoric(1,1)<<
"\n =========================================="<<
"\n feedback matrix "<<
"\n"<<feedback(0,0)<<"  "<<feedback(0,1)<<
"\n"<<feedback(1,0)<<"  "<<feedback(1,1)<<std::endl;
*/
//growth rate of faces : randomized number between (kappa-0.5 to kappa + 0.5)
double growthfactor = this->getGrowthRandomNumber();
this->lastGrowthRate = growthfactor; //saving growth rate for plotting
/*
std::cout<<"Kappa :: "<< kappa <<"\n Actual Growth Var  : "<<growthvar <<std::endl;
std::cout<<"Feedback parameter : Eta :: " << eta <<std::endl;
*/
Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eigensolver;
// Growth Matrix
Eigen::Matrix2d growthMatrix;
// dM0/dt = kappa*M0 - n/2*feedback
growthMatrix<< growthfactor*(MG(0,0)) - eta/2.*feedback(0,0),
               growthfactor*(MG(0,1)) - eta/2.*feedback(0,1),
               growthfactor*(MG(1,0)) - eta/2.*feedback(1,0),
               growthfactor*(MG(1,1)) - eta/2.*feedback(1,1); 
eigensolver.compute(growthMatrix);//computing the eigenvalues of growthMatrix, to make sure it is always growing
//to calculate the individual growth eigen direction
  Eigen::Matrix2d eigen1;
  Eigen::Matrix2d eigen2;
  eigen1 = std::max(eigensolver.eigenvalues()[0]-cell->thresholdMatrix[0][0],0.0)*
                      ((eigensolver.eigenvectors().col(0))*(eigensolver.eigenvectors().col(0)).transpose());
  eigen2 = std::max(eigensolver.eigenvalues()[1]-cell->thresholdMatrix[1][1],0.0)*
                      ((eigensolver.eigenvectors().col(1))*(eigensolver.eigenvectors().col(1)).transpose());
// Now combining growth in both direction
  growthMatrix = eigen1 + eigen2;
//now setting the new targetFormMatrix with feedback matrix
this->targetFormMatrix[0][0] += growthMatrix(0,0);
this->targetFormMatrix[1][0] += growthMatrix(1,0);
this->targetFormMatrix[0][1] += growthMatrix(0,1);
this->targetFormMatrix[1][1] += growthMatrix(1,1);
/*
std::cout<<
"\n =========================================="<<
"\n Face ID : : "<< this->getID() <<
"\n eta : : "<< eta <<
"\n =========================================="<<
"\n Growth matrix "<<
"\n"<<growthMatrix(0,0)<<"  "<<growthMatrix(0,1)<<
"\n"<<growthMatrix(1,0)<<"  "<<growthMatrix(1,1)<<
"\n =========================================="<<
"\n New Target Form matrix "<<
"\n"<<this->targetFormMatrix[0][0]<<"  "<<this->targetFormMatrix[0][1]<<
"\n"<< this->targetFormMatrix[1][0]<<"  "<< this->targetFormMatrix[1][1]<<std::endl;
*/
//now setting tracesq of Target Form Matrix
this->setTraceSquaredTargetFormMatrix();
}

// *************************************************************** //
 
 void Face::inflatedGrow(){
  if (this->getID() == 1){
    return;
  }
  
//growth rate of faces : randomized number between (kappa-0.5 to kappa + 0.5)
double growthfactor = this->getGrowthRandomNumber();
lastGrowthRate = growthfactor; //saving growth rate for plotting
  //std::cout<<"Kappa : "<< kappa <<"/n Actual Growth Var  : "<<growthvar <<std::endl;

  /*
  std::cout<<"Growth Rate addition to TargetFormMatrix : (from Face::grow() : Faceid :"<<this->getID()<<std::endl;
  std::cout<<growthRate<<std::endl;
  */
  //now setting the new targetFormMatrix
  this->targetFormMatrix[0][0] += growthfactor*(this->targetFormMatrix[0][0]);
  this->targetFormMatrix[1][0] += growthfactor*(this->targetFormMatrix[1][0]);
  this->targetFormMatrix[0][1] += growthfactor*(this->targetFormMatrix[0][1]);
  this->targetFormMatrix[1][1] += growthfactor*(this->targetFormMatrix[1][1]);
  //now setting tracesq of Target Form Matrix
  this->setTraceSquaredTargetFormMatrix();
 }




// *************************************************************** //
 void Face::setDivisionThreshold(){
     this->divisionThreshold = (this->divisionFactor)*(this->getAreaOfFace());
 }
 // ************************************************************** //
 void Face::divideRandom(){
    if (this->id == 1){
      return;
    };
    // Checking if the Area of Face has crossed the division Threshold or not
    // If the Area is greater than threshold then the cell can divide
    if (this->getAreaOfFace() < this->divisionThreshold){
        return;
    }
    //getting cell of this face
    Cell * cell = this->getCell();
    double phi = M_PI*(cell->getCellDivisionRandomNumber());//get a random angle for the cells to divide
    unsigned int faceid = this->getID();
    double vdirect[2] = {cos(phi), sin(phi)}; //taking the direction of the new wall
    Edge * edge;
    Edge * intersectedEdge[2];
    Vertex * vertexA;
    Vertex * vertexB;
    double Cx(0.), Cy(0.), Ax,Ay,Bx,By,s;
    double xintersect, yintersect, zintersect;
    int counter;
    vertex_coordinate intersection_points[2];
    {//iterating through the edge to see the division
      FaceEdgeIterator edges(this);
      counter = 0;
      while ((edge = edges.next())!= 0){
            vertexA = edge->Org();
            vertexB = edge->Dest();
            //getting coordinates of A and B
            Ax = vertexA->getProjectedXcoordinate(faceid);
            Ay = vertexA->getProjectedYcoordinate(faceid);
            Bx = vertexB->getProjectedXcoordinate(faceid);
            By = vertexB->getProjectedYcoordinate(faceid);
            //now calculating the intersection of the old walls to new walls
            //calculating the denominator
            double denominator = (By-Ay)*vdirect[0]-(Bx-Ax)*vdirect[1];
            s = ((Cy-Ay)*vdirect[0]-(Cx-Ax)*vdirect[1])/denominator;
            // now checking if the new cellwall intersects the current edge
            if ((s>0) && (s<1)){
                assert(counter<2);
                intersectedEdge[counter] = edge;
                xintersect = vertexA->getXcoordinate() + s*(vertexB->getXcoordinate()-vertexA->getXcoordinate());
                yintersect = vertexA->getYcoordinate() + s*(vertexB->getYcoordinate()-vertexA->getYcoordinate());
                zintersect = vertexA->getZcoordinate() + s*(vertexB->getZcoordinate()-vertexA->getZcoordinate());
                //assinging the intersection points to the intertion_points array
                intersection_points[counter].x = xintersect;
                intersection_points[counter].y = yintersect;
                intersection_points[counter].z = zintersect;
                counter += 1;
            }
        }
      }
    //now we have two edges that are taking part in division
    Face * left, * right;
    Vertex * origin;
    Vertex * addedVertex[2];
    for (int i = 0; i <2 ; i++){
        left = intersectedEdge[i]->Left();
        right = intersectedEdge[i]->Right();
        origin = intersectedEdge[i]->Org();
        //making new edge on the origin of this edge
        edge =cell->makeVertexEdge(origin, left,right);
        origin = edge->Dest();
        origin->setXcoordinate(intersection_points[i].x);
        origin->setYcoordinate(intersection_points[i].y);
        origin->setZcoordinate(intersection_points[i].z);
        addedVertex[i] = origin;
        }
    // Now making the connecting vertex for new wall
    //                   v2
    //    left = this     |    right = new face
    //                    v1
    edge = cell->makeFaceEdge(this, addedVertex[0], addedVertex[1]); 
    right = edge->Right();//new face created
    vertexA = edge->Org();
    vertexB = edge->Dest();
    //making the new vertex on the wall now.
    edge =  cell->makeVertexEdge(vertexA, edge->Left(), edge->Right());
    origin = edge->Dest();
    //now setting hte vertex of this new vertex
    origin->setXcoordinate((vertexA->getXcoordinate()+vertexB->getXcoordinate())/2.);
    origin->setYcoordinate((vertexA->getYcoordinate()+vertexB->getYcoordinate())/2.);
    origin->setZcoordinate((vertexA->getZcoordinate()+vertexB->getZcoordinate())/2.);
    // saving the needed quantitites for the calculation of new Form Matrix of daugther cells
    // Daughter cells share same quantities, inherited from Mother cell
    // Newest Daughter cell (new face created above) inherits the same properties from Mother cell 
    // Daughter cell 2
    right->mu1 = this->mu1;
    right->mu2 = this->mu2;
    right->mu3 = this->mu3;
    right->mu4 = this->mu4;
    right->currentFormMatrix[0][0] = this->currentFormMatrix[0][0];
    right->currentFormMatrix[0][1] = this->currentFormMatrix[0][1];
    right->currentFormMatrix[1][0] = this->currentFormMatrix[1][0];
    right->currentFormMatrix[1][1] = this->currentFormMatrix[1][1];
    //Also updating the division threshold and target Area for the new daughter cell to be same as mother cells
    right->setDivisionThreshold(this->divisionThreshold);
    right->setTargetArea(this->getAreaOfFace());
    right->setKappa(this->getKappa()); // Updating the Daughter cell with same growth rate as  Mother cell
    right->setGrowthVar(this->getGrowthVar());
    // Updating this cell's terms for Form Matrix calcualtion
    this->oldMuMatrix[0][0] = this->mu1;
    this->oldMuMatrix[0][1] = this->mu2;
    this->oldMuMatrix[1][0] = this->mu3;
    this->oldMuMatrix[1][1] = this->mu4;
    this->oldcurrentFormMatrixTrace = this->currentFormMatrix[0][0]+this->currentFormMatrix[1][1];
    // Now updating this cells neighbours to calculate the Form Matrix 
    {
      FaceEdgeIterator edges(this);
      while ((edge = edges.next())!= 0 ){
            right = edge->Right();//getting the right face
            if (right->getID() == 1) continue;// continue if the right face is the external face
            right->oldMuMatrix[0][0] = right->mu1;
            right->oldMuMatrix[0][1] = right->mu2;
            right->oldMuMatrix[1][0] = right->mu3;
            right->oldMuMatrix[1][1] = right->mu4;
            right->oldcurrentFormMatrixTrace = right->currentFormMatrix[0][0]+right->currentFormMatrix[1][1];
      }
    }
    // now as all the needed old need measures are stored 
    //we update the parameters for this face and the neighbouring cells
    this->setDivideFormMatrix();
    cell->countCellDivision();

 }
 // *************************************************************** // 
 void Face::setDivideFormMatrix(){
  // update the parameters of all vertex and faces of this cell
  Cell * cell = this->getCell();
  cell->setParameters();
  
  //updating Form Matrix of this face
  double traceRatio(0.);
  //right = this;
  traceRatio = (this->currentFormMatrix[0][0]+this->currentFormMatrix[1][1])/(this->oldcurrentFormMatrixTrace);
  this->targetFormMatrix[0][0] = this->currentFormMatrix[0][0] - traceRatio*(this->oldMuMatrix[0][0]);
  this->targetFormMatrix[0][1] = this->currentFormMatrix[0][1] - traceRatio*(this->oldMuMatrix[0][1]);
  this->targetFormMatrix[1][0] = this->currentFormMatrix[1][0] - traceRatio*(this->oldMuMatrix[1][0]);
  this->targetFormMatrix[1][1] = this->currentFormMatrix[1][1] - traceRatio*(this->oldMuMatrix[1][1]);
  // updating the trace squared 
  this->setTraceSquaredTargetFormMatrix();

  //updating Form Matrix of neighbouring faces
  {
      Face * newright;
      FaceEdgeIterator newedges(this);
      Edge * newedge;
      while ((newedge = newedges.next())!= 0 ){
            newright = newedge->Right();//getting the right face
            if (newright->getID() == 1) continue;// continue if the right face is the external face
            traceRatio = (newright->currentFormMatrix[0][0]+newright->currentFormMatrix[1][1])/(newright->oldcurrentFormMatrixTrace);
            newright->targetFormMatrix[0][0] = newright->currentFormMatrix[0][0] - traceRatio*(newright->oldMuMatrix[0][0]);
            newright->targetFormMatrix[0][1] = newright->currentFormMatrix[0][1] - traceRatio*(newright->oldMuMatrix[0][1]);
            newright->targetFormMatrix[1][0] = newright->currentFormMatrix[1][0] - traceRatio*(newright->oldMuMatrix[1][0]);
            newright->targetFormMatrix[1][1] = newright->currentFormMatrix[1][1] - traceRatio*(newright->oldMuMatrix[1][1]);
            // updating trace squared
            newright->setTraceSquaredTargetFormMatrix();
      }
    }
    // Now setting final parameters for this total tissue
    cell->setParameters();
 }
  //****************** end added features********************************//
/* -- protected instance methods ------------------------------------------- */

Face::Face(Cell *cell):gaussianWidth(0.125), randomNumberGeneratorType(gsl_rng_default) 
{
  assert(cell!=0);

  this->cell = cell;
  this->id   = cell->makeFaceID();
  this->data = 0;
  this->edge = 0;
  this->sumEdgeLength = 0;
  //***************added features*******************************************//
  this->vertices = new Vertex*[8];
  this->vertexCount = 0;
  this->vertexSize = 8;
  this->divisionFactor = cell->getDivisionFactor();
  this->domePosition = true;//seting position to dome as True in default
  this->lastGrowthRate = 0.;
  this->targetArea = 1.;
  this->growthVar = cell->getGrowthVar();
  this->constantGrowthMatrix[0][0] = 0.;
  this->constantGrowthMatrix[0][1] = 0.;
  this->constantGrowthMatrix[1][0] = 0.;
  this->constantGrowthMatrix[1][1] = 0.;
  //setting the random number generator
  // intialised in Initialising list :-> randomNumberGeneratorType = gsl_rng_default;//this is Mersenne Twister algorithm
  randomNumberGenerator = gsl_rng_alloc(randomNumberGeneratorType);
  gsl_rng_set(randomNumberGenerator,cell->getSeedRandomInteger());//some number as seed-> this can be set with another random number/time or /dev/random /urandom

  // RANDOM NUMBER GENERATOR FOR CELL DIVISION
  cellDivisionRandomNumberGenerator = gsl_rng_alloc(randomNumberGeneratorType);
  //gsl_rng_set(cellDivisionRandomNumberGenerator, random_seed());//using /dev/urandom to seed this generator
  gsl_rng_set(cellDivisionRandomNumberGenerator, cell->getCellDivisionRandomNumber());
  //calculating the square root of epsilon
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
  //release the randomNumberGenerator
  gsl_rng_free(randomNumberGenerator);
  gsl_rng_free(cellDivisionRandomNumberGenerator);
  //***************end added features*************************************//
  cell->removeFace(this);
}

