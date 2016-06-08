
/* ============================================================================
 * p2/cell/face.cc
 * ========================================================================= */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>//for math operations sqrt and others

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
void Face::setCentralisedCoordinate(){
  FaceEdgeIterator edges(this);//iterator to iterate through the face
  Edge *newedge;//a pointer to keep current edge
  double area = this->areaOfFace;//area of face
  unsigned int faceid = this->id;//getting the face id
  //temporary variables to calculate centralised Coordinates
  double xTemp = 0;
  double yTemp = 0;
  //temporary vertex pointer
  Vertex *vertexOrg;
  Vertex *vertexDest;
  //now looping through all the edges of this face and calculating the cnetralised coordinate
  while ((newedge = edges.next())!=0){//until exhausted
    vertexOrg = newedge->Org();
    vertexDest = newedge->Dest();
    xTemp += (vertexOrg->getProjectedXcoordinate(id)+vertexDest->getProjectedXcoordinate(id))*
              vertexOrg->getAk(faceid);//summing (x_i + x_i+1)*Ak
    yTemp += (vertexOrg->getProjectedYcoordinate(id)+vertexDest->getProjectedYcoordinate(id))*
              vertexOrg->getAk(faceid);//summing (y_i + y_i+1)*Ak
  }
  this->xCentralised = xTemp/(6*area);
  this->yCentralised = yTemp/(6*area); 
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
void Face::setProjectedCoordinate(){
    //setting the projected coordinates of all vertices in this face
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
    // getting the mean center position of this face //
    unsigned int faceid = this->getID();// first grabbing id of  current face
    if (faceid == 1){//if faceid == 1 then this is external face, that doesnot need projection or calculation
        return;
    } 
    FaceEdgeIterator faceEdges(this);//iterator to iterate through the vertex for outgoign edge
    double nx,ny;//normal coordinate
    double xCentroid(0), yCentroid(0), zCentroid(0), totalarea(0); // coordinate of the cnetroid
    // array of vertices
    double xcood[6], ycood[6], zcood[6];// coordinate of the vertices of this face
    double xTriCen[6], yTriCen[6], zTriCen[6], areaTri[6];// coordinate of center of triangluated triangles of face
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
    // divinding by the number of vertices to get the mean
    xmean = (1./counter)*xmean;
    ymean = (1./counter)*ymean;
    zmean = (1./counter)*zmean;
    //printf(">>>>new output start<<<< face id = %u \n", faceid);
    //printf("calculated means : X = %F ; Y = %F ; Z = %F \n", xmean, ymean,zmean );
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
    // triangulating the face and getting the center of each triangle
    //vectors of the triangles 
    double vector1[3],vector2[3];
    double crossProductVector[3];//cross product of the vector
    double crossProductMagnitude;//cross product magintude
    for (int counter = 0; counter<6; counter++){
          // calculating center of the triangles
          xTriCen[counter] = 1./3.*(xmean+xcood[counter]+xcood[(counter+1)%6]);
          yTriCen[counter] = 1./3.*(ymean+ycood[counter]+ycood[(counter+1)%6]);
          zTriCen[counter] = 1./3.*(zmean+zcood[counter]+zcood[(counter+1)%6]);
          //printf("counter = %d; center of triangles : %F ; %F ; %F \n", counter, xTriCen[counter], yTriCen[counter], zTriCen[counter]);
          //getting two vectors of this triangle
          vector1[0] = xcood[counter]-xTriCen[counter];
          vector1[1] = ycood[counter]-yTriCen[counter];
          vector1[2] = zcood[counter]-zTriCen[counter];
          vector2[0] = xcood[(counter+1)%6]-xTriCen[counter];
          vector2[1] = ycood[(counter+1)%6]-yTriCen[counter];
          vector2[2] = zcood[(counter+1)%6]-zTriCen[counter];
          //cross product of the two vectors
          crossProductVector[0] = vector1[1]*vector2[2] - vector1[2]*vector2[1];
          crossProductVector[1] = vector1[2]*vector2[0] - vector1[0]*vector2[2];
          crossProductVector[2] = vector1[0]*vector2[1] - vector1[1]*vector2[0];
          //maginitude of cross product
          crossProductMagnitude = sqrt(crossProductVector[0]*crossProductVector[0] + crossProductVector[1]+crossProductVector[1]+
                                  crossProductVector[2]*crossProductVector[2]);
          areaTri[counter] = 0.5*crossProductMagnitude;//area of this triangle
          xCentroid += xTriCen[counter]*areaTri[counter];//adding the weigthed centroid
          yCentroid += yTriCen[counter]*areaTri[counter];
          zCentroid += zTriCen[counter]*areaTri[counter]; 
          //adding total area
          totalarea += areaTri[counter];// calculating total area of triange
      }
    //now calculating weighted center
    xCentroid = xCentroid/totalarea;
    yCentroid = yCentroid/totalarea;
    zCentroid = zCentroid/totalarea;
    //printf("%s : ", "calculating total area");
    //printf("%F", totalarea );
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
    //Calculating normal of the triangles of the shape and then calculating weighted normal
    double normalTempX[6],normalTempY[6],normalTempZ[6];//vertices of normal of triangles
    double normalX(0), normalY(0), normalZ(0); // vertices of weighted normal X, Y, Z 
    for (int counter = 0; counter<6; counter++){//triangulating the vertices with centroid of this face
          //getting two vectors of this triangle
          vector1[0] = xcood[counter]-xCentroid;
          vector1[1] = ycood[counter]-yCentroid;
          vector1[2] = zcood[counter]-zCentroid;
          vector2[0] = xcood[(counter+1)%6]-xTriCen[counter];
          vector2[1] = ycood[(counter+1)%6]-yTriCen[counter];
          vector2[2] = zcood[(counter+1)%6]-zTriCen[counter];
          //cross product of the two vectors
          crossProductVector[0] = vector1[1]*vector2[2] - vector1[2]*vector2[1];
          crossProductVector[1] = vector1[2]*vector2[0] - vector1[0]*vector2[2];
          crossProductVector[2] = vector1[0]*vector2[1] - vector1[1]*vector2[0];
          //maginitude of cross product
          crossProductMagnitude = sqrt(crossProductVector[0]*crossProductVector[0] + crossProductVector[1]+crossProductVector[1]+
                                  crossProductVector[2]*crossProductVector[2]);
          areaTri[counter] = 0.5*crossProductMagnitude;//area of this triangle
          //adding the weigthed centroid 
          normalX += crossProductVector[0]*areaTri[counter];
          normalY += crossProductVector[1]*areaTri[counter];
          normalZ += crossProductVector[2]*areaTri[counter]; 
          //total area is same as before calculated for this face
      }
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
    //now projecting the vertices on to the plane
    FaceEdgeIterator faceEdges1(this);//iterator to iterate through the vertex for outgoign edge
    double xvertex, yvertex, zvertex,vectorVertex[3],dotproduct;
    //projected vertices: xproj is 6 element array, with element 1 giving x coordinate of vertex 1, same with yproj, zproj
    double xproj[6], yproj[6], zproj[6];
    int counterproj(0); //counter 
    while((currentEdge = faceEdges1.next())!= 0){//runnign through the edges in the face again
          currentVertex = currentEdge->Dest();// the current vertex, iterated
          xvertex = currentVertex->getXcoordinate();
          yvertex = currentVertex->getYcoordinate();
          zvertex = currentVertex->getZcoordinate();
          //getting the vector form the centroid
          vectorVertex[0] = xvertex-xCentroid;
          vectorVertex[1] = yvertex-yCentroid;
          vectorVertex[2] = zvertex-zCentroid;
          //dot product of vectorVertex and normal
          dotproduct = vectorVertex[0]*normalX+vectorVertex[1]*normalY+vectorVertex[2]*normalZ;
          //now calculating the projected vertices 
          xproj[counterproj]  = xvertex - dotproduct*normalX;
          yproj[counterproj]  = yvertex - dotproduct*normalY;
          zproj[counterproj]  = zvertex - dotproduct*normalZ;
          counterproj += 1;
        }
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
    //now need to get this projected 3d plane coordinates onto rotated 2d coordinate system
    // getting the unit vector of 2d Plane
    double unitx[3], unity[3];//new unit vector of x, y on the plane
    //getting the unitx = normalised[Projectedvertex1 - Origin]
    unitx[0] = xproj[0]- xCentroid;
    unitx[1] = yproj[0]- yCentroid;
    unitx[2] = zproj[0]- zCentroid;
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
    double xprojection, yprojection;
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
          //now setting the projected coordinates in the vertex properties
          currentVertex->insertProjectedXcoordinate(faceid,xprojection);
          currentVertex->insertProjectedYcoordinate(faceid,yprojection); 
          //getting the vector form the real origin on new vertex
          // NON CENTRALISED Projected Coordainte
          vectorVertex[0] = xvertex;
          vectorVertex[1] = yvertex;
          vectorVertex[2] = zvertex;
          //now getting the new x, y coordinates, dot product of unit vector with the vectorVertex
          xprojection = unitx[0]*vectorVertex[0]+unitx[1]*vectorVertex[1]+unitx[2]*vectorVertex[2];
          yprojection = unity[0]*vectorVertex[0]+unity[1]*vectorVertex[1]+unity[2]*vectorVertex[2];
          //now setting the projected coordinates in the vertex properties
          currentVertex->insertNonCentralisedProjectedXcoordinate(faceid,xprojection);
          currentVertex->insertNonCentralisedProjectedYcoordinate(faceid,yprojection); 
        }
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
      this->setTargetFormMatrix();
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
      this->mu1 = term11 - targetFormMatrix[0][0];
      this->mu2 = term11 - targetFormMatrix[1][0];
      this->mu3 = this->mu2;//as mu2 and mu3 are equal
      this->mu4 = term11 - targetFormMatrix[1][1];
  }
  void Face::setTargetFormMatrix(){
    //*******************************************//
    // this should be edited to include time update
    // at each growth step
    //*******************************************//
     //iterating this face for the vertices
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
  }
 void Face::setTraceSquaredTargetFormMatrix(){
    double lambda1 = this->targetFormMatrix[0][0];
    double lambda2 = this->targetFormMatrix[1][1];
    double trace =  (lambda1+lambda2)*(lambda1+lambda2);
    this->traceSquaredTargetFormMatrix = trace;
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

