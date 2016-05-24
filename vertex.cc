
/* ============================================================================
 * p2/cell/vertex.cc
 * ========================================================================= */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

//#include <svl/Vec3.h>
#include "svl-1.5/include/svl/SVL.h"
#include "cell.hh"
#include "edge.hh"
#include "vertex.hh"

/* ----------------------------------------------------------------------------
 * Vertex
 * ------------------------------------------------------------------------- */

/* -- public class methods ------------------------------------------------- */

Vertex *Vertex::make(Cell *cell)
{
  assert(cell!=0);

  return new Vertex(cell);
}

void Vertex::kill(Vertex *vertex)
{
  assert(vertex!=0);
  delete vertex;
}

/* -- public instance methods ---------------------------------------------- */
void Vertex::setXcoordinate(double xcoord){
	Xcoordinate = xcoord;
}
void Vertex::setYcoordinate(double ycoord){
	Ycoordinate = ycoord;
}
void Vertex::setZcoordinate(double zcoord){
	Zcoordinate = zcoord;
}
double Vertex::getXcoordinate(){
	return Xcoordinate;
}
double Vertex::getYcoordinate(){
	return Ycoordinate;
}
double Vertex::getZcoordinate(){
	return Zcoordinate;
}
void Vertex::insertProjectedXcoordinate(unsigned int faceid, double xcood){
	assert(faceid>0);
	this->ProjectedXcoordinate[faceid] = xcood;
}
///////////////////////////////////////////////////////////////////////
void Vertex::insertProjectedYcoordinate(unsigned int faceid, double ycood){
	assert(faceid>0);
	this->ProjectedYcoordinate[faceid] = ycood;
}
///////////////////////////////////////////////////////////////////////
double Vertex::getProjectedXcoordinate(unsigned int faceid){
	return this->ProjectedXcoordinate[faceid];
}
///////////////////////////////////////////////////////////////////////
double Vertex::getProjectedYcoordinate(unsigned int faceid){
	return this->ProjectedYcoordinate[faceid];
}

///////////////////////////////////////////////////////////////////////

void Vertex::setID(unsigned int id)
{
  assert(id>0);

  this->id = id;
}

void Vertex::addEdge(Edge *edge)
{
  assert(edge!=0);

  // only keep track of one edge in the orbit--this one is as good as any

  this->edge = edge;
}

void Vertex::removeEdge(Edge *edge)
{
  assert(edge!=0);

  // replace the arbitrary edge with another edge in the orbit
  // use null if this is the only edge
  // assumes that the edge hasn't been actually removed yet

  Edge *next = edge->Onext();

  this->edge = next!=edge ? next : 0;
}

/* -- protected instance methods ------------------------------------------- */

Vertex::Vertex(Cell *cell)
{
  assert(cell!=0);

  this->pos[0] = 0.0;
  this->pos[1] = 0.0;
  this->pos[2] = 0.0;
  this->cell   = cell;
  this->id     = cell->makeVertexID();
  this->data   = 0;
  this->edge   = 0;

  cell->addVertex(this);
}

Vertex::~Vertex()
{
  cell->removeVertex(this);
}

//========================Derivatives===============================//
void Vertex::setAreaDerivative(){
    VertexEdgeIterator vertEdges(this);//iterator to iterate through the vertex for outgoign edge
    Edge *currentEdge;//pointer to keep track of a edge that is iterated 
    while ((currentEdge = vertEdges.next())!=0){//iterating till exhausted
        Face *left = currentEdge->Left();//grabbing the inner face of the edge
        Edge *prevEdge = currentEdge->Lprev();//edge that is previous on the same left face
        //grabbing j+1 and j-1 vertex
        Vertex *vertex1 = prevEdge->Org();
        Vertex *vertex2 = currentEdge->Dest();
        //setting the dervatives now with key of left face id (inner face id)
        unsigned int innerid = left->getID();
        this->areaXDerivative[innerid] = 0.5*(vertex2->getProjectedYcoordinate(innerid)-vertex1->getProjectedYcoordinate(innerid));//y_{j+1}-y_{j-1}
        this->areaYDerivative[innerid] = 0.5*(vertex1->getProjectedXcoordinate(innerid)-vertex2->getProjectedXcoordinate(innerid));//x_{j-1}-x_{j+1}
    }
}

void Vertex::setAkDerivative(){
    VertexEdgeIterator vertEdges(this);//iterator to iterate through the vertex for outgoign edge
    Edge *currentEdge;//pointer to keep track of a edge that is iterated 
    while ((currentEdge = vertEdges.next())!=0){//iterating till exhausted
        Face *left = currentEdge->Left();//grabbing the inner face of the edge
        Edge *prevEdge = currentEdge->Lprev();//edge that is previous on the same left face
        //grabbing j+1 and j-1 vertex
        Vertex *vertex1 = prevEdge->Org();
        Vertex *vertex2 = currentEdge->Dest();
        //setting the dervatives now with key of left face id (inner face id)
        unsigned int innerid = left->getID();
        this->AkXDerivative[innerid] = vertex2->getProjectedYcoordinate(innerid)-vertex1->getProjectedYcoordinate(innerid);//y_{j+1}-y_{j-1}
        this->AkYDerivative[innerid] = vertex1->getProjectedXcoordinate(innerid)-vertex2->getProjectedXcoordinate(innerid);//x_{j-1}-x_{j+1}
    }
}

double Vertex::getAreaXDerivative(unsigned int faceid){
    return this->areaXDerivative[faceid];
}
double Vertex::getAreaYDerivative(unsigned int faceid){
    return this->areaYDerivative[faceid];
}
double Vertex::getAkXDerivative(unsigned int faceid){
    return this->AkXDerivative[faceid];
}

double Vertex::getAkYDerivative(unsigned int faceid){
    return this->AkYDerivative[faceid];
}

// ******************************************************* //
void Vertex::setFunction1(){
  VertexEdgeIterator vertEdges(this);//iterator to iterate through the vertex for outgoign edge
  double functionvalue, yk, yk1; //variable to store function value
    Edge *currentEdge;//pointer to keep track of a edge that is iterated 
    while ((currentEdge = vertEdges.next())!=0){//iterating till exhausted
        Face *left = currentEdge->Left();//grabbing the inner face of the edge
        unsigned int innerid = left->getID();//id of inner face
        Edge *prevEdge = currentEdge->Lprev();//edge that is previous on the same left face
        //grabbing j+1 and j-1 vertex
        //Vertex *vertex1 = prevEdge->Org();//j-1
        Vertex *vertex2 = currentEdge->Org();//j
        Vertex *vertex3 = currentEdge->Dest();//j+1
        yk = vertex2->getProjectedYcoordinate(innerid);
        yk1 = vertex3->getProjectedYcoordinate(innerid);
        //calculating the fucntion
        functionvalue = yk*yk + yk*yk1+yk1*yk1;
        //putting function1 value in the array
        this->Function1[innerid] = functionvalue;
  }
}
// ******************************************************* //
void Vertex::setFunction2(){
  VertexEdgeIterator vertEdges(this);//iterator to iterate through the vertex for outgoign edge
  double functionvalue, xk, xk1, yk, yk1; //variable to store function value
    Edge *currentEdge;//pointer to keep track of a edge that is iterated 
    while ((currentEdge = vertEdges.next())!=0){//iterating till exhausted
        Face *left = currentEdge->Left();//grabbing the inner face of the edge
        unsigned int innerid = left->getID();//id of inner face
        Edge *prevEdge = currentEdge->Lprev();//edge that is previous on the same left face
        //grabbing j+1 and j-1 vertex
        //Vertex *vertex1 = prevEdge->Org();//j-1
        Vertex *vertex2 = currentEdge->Org();//j
        Vertex *vertex3 = currentEdge->Dest();//j+1
        xk = vertex2->getProjectedXcoordinate(innerid);
        xk1 = vertex3->getProjectedXcoordinate(innerid);
        yk = vertex2->getProjectedYcoordinate(innerid);
        yk1 = vertex3->getProjectedYcoordinate(innerid);
        //calculating the fucntion
        functionvalue = xk*yk1 + 2*xk*yk + 2*xk1*yk1 + xk1*yk;
        //putting function1 value in the array
        this->Function2[innerid] = functionvalue;
  }
}
// ******************************************************* //
void Vertex::setFunction3(){
  VertexEdgeIterator vertEdges(this);//iterator to iterate through the vertex for outgoign edge
  double functionvalue, xk, xk1; //variable to store function value
    Edge *currentEdge;//pointer to keep track of a edge that is iterated 
    while ((currentEdge = vertEdges.next())!=0){//iterating till exhausted
        Face *left = currentEdge->Left();//grabbing the inner face of the edge
        unsigned int innerid = left->getID();//id of inner face
        Edge *prevEdge = currentEdge->Lprev();//edge that is previous on the same left face
        //grabbing j+1 and j-1 vertex
        //Vertex *vertex1 = prevEdge->Org();//j-1
        Vertex *vertex2 = currentEdge->Org();//j
        Vertex *vertex3 = currentEdge->Dest();//j+1
        xk = vertex2->getProjectedXcoordinate(innerid);
        xk1 = vertex3->getProjectedXcoordinate(innerid);
        //calculating the fucntion
        functionvalue = xk*xk + xk*xk1+xk1*xk1;
        //putting function1 value in the array
        this->Function3[innerid] = functionvalue;
    }
}
// ******************************************************* //
void Vertex::setAk(){
  VertexEdgeIterator vertEdges(this);//iterator to iterate through the vertex for outgoign edge
  double functionvalue, xk, xk1, yk, yk1; //variable to store function value
    Edge *currentEdge;//pointer to keep track of a edge that is iterated 
    while ((currentEdge = vertEdges.next())!=0){//iterating till exhausted
        Face *left = currentEdge->Left();//grabbing the inner face of the edge
        unsigned int innerid = left->getID();//id of inner face
        Edge *prevEdge = currentEdge->Lprev();//edge that is previous on the same left face
        //grabbing j+1 and j-1 vertex
        //Vertex *vertex1 = prevEdge->Org();//j-1
        Vertex *vertex2 = currentEdge->Org();//j
        Vertex *vertex3 = currentEdge->Dest();//j+1
        xk = vertex2->getProjectedXcoordinate(innerid);
        xk1 = vertex3->getProjectedXcoordinate(innerid);
        yk = vertex2->getProjectedYcoordinate(innerid);
        yk1 = vertex3->getProjectedYcoordinate(innerid);
        //calculating the fucntion
        functionvalue = xk*yk1 - xk1*yk;
        //putting function1 value in the array
        this->Ak[innerid] = functionvalue;
  }
}
