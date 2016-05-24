
/* ============================================================================
 * p2/cell/face.cc
 * ========================================================================= */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

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
      double term11, term12, term22; // to store the values of summing terms to calculate mu
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
    double term11, term12, term22; // to store the values of matrix element
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

