
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
    if (faceid == 1){//if faceid == 1 then this is external face, that doesnot need projection or calculation
        // setting the central coordinate of this face in terms of cartisian coordinate
          this->xCentralised = 0.;
          this->yCentralised = 0.;
          this->zCentralised = 0.;
        return;
    } 
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
    //new way of getting intrinsic x axis
    //double wvector[3];
    //If normal is parralel to X-axis (or in this normalised case, eqaul to x-axis) then the projection is done using Y-Axis
    if (((normalX - 1)<0.001) and((normalY - 0)<0.001) and ((normalZ - 0)<0.001) ){
        vectorVertex[0] = 0.;
        vectorVertex[1] = 1.;
        vectorVertex[2] = 0.;
        //std::cout<<"normal X"<<std::endl;   
    } else {
        vectorVertex[0] = 1.;
        vectorVertex[1] = 0.;
        vectorVertex[2] = 0.;  
        //std::cout<<"normal Y"<<std::endl;
      };

    //to get the unit vector in X direction, lets project X axis (1,0,0) on to the plane
    /*
    double cross_u[3];
    //cross product of the two vectors
    cross_u[0] = wvector[1]*normalZ - wvector[2]*normalY;
    cross_u[1]=  wvector[2]*normalX - wvector[0]*normalZ;
    cross_u[2] = wvector[0]*normalY - wvector[1]*normalX;
    //final cross - product to given vector perpendicular to normal vector,
    //cross product of the two vectors
    wvector[0] = cross_u[1]*normalZ - cross_u[2]*normalY;
    wvector[1]=  cross_u[2]*normalX - cross_u[0]*normalZ;
    wvector[2] = cross_u[0]*normalY - cross_u[1]*normalX;
    */
    //vector from (1,0,0)+(xCentroid,Ycentroid,Zcentroid) until xCentroid
    //vectorVertex[0] = 1.;//1.+xCentroid-xCentroid; just directly writing the resulting answer
    //vectorVertex[1] = 0.;//0.+yCentroid-yCentroid;
    //vectorVertex[2] = 0.;//0.+zCentroid-zCentroid;
    //dot product of vectorVertex and normal, which is equal to normalX; this is n^c_x = \phi_P
    dotproduct = normalX;//vectorVertex[0]*normalX+vectorVertex[1]*normalY+vectorVertex[2]*normalZ;
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
          // NON CENTRALISED Projected Coordainte or coordinates in terms of CARTESIAN SYSTEM
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
    this->setTraceSquaredTargetFormMatrix();
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

void Face::setTempTargetFormMatrixCurrent(){
    this->targetFormMatrix[0][0] = (this->getMu1() + this->targetFormMatrix[0][0] );
    this->targetFormMatrix[1][0] = (this->getMu2() + this->targetFormMatrix[1][0] );
    this->targetFormMatrix[0][1] = (this->getMu3() + this->targetFormMatrix[0][1] );
    this->targetFormMatrix[1][1] = (this->getMu4() + this->targetFormMatrix[1][1] );
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
      double alpha = cell->getAlpha();
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
      double energytemp = alpha*firstterm + beta*secondterm - pressure*thirdterm;
      //removing the area term
      //double energytemp = alpha*firstterm + beta*secondterm;
      //setting the energy values of face 
      this->firstTerm = firstterm;
      this->secondTerm = secondterm;
      this->thirdTerm = thirdterm;
      this->energy = energytemp;
  }
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
  //growth fluctuation : calculated by same randomnumber generator set as property of cell
  double fluctuation = cell->getRandomNumber();
  //growth variation of face : Amplitude of fluctuation
  double growthvar = cell->getGrowthVar();
  //growth rate of faces : kappa
  double kappa = cell->getKappa();
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
  growthRate = kappa*(1+(2*growthvar*fluctuation-growthvar))*(eigen1+eigen2);
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
    //Also updating the division threshold for the new daughter cell to be same as mother cells
    right->setDivisionThreshold(this->divisionThreshold);
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

      }
    }
    // Now setting final parameters for this total tissue
    cell->setParameters();
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
  this->divisionFactor = cell->getDivisionFactor();
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

