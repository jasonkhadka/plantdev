
/* ============================================================================
 * p2/cell/edge.cc
 * ========================================================================= */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <map>
#include "edge.hh"
#include "face.hh"
#include "vertex.hh"

/* ----------------------------------------------------------------------------
 * QuadEdge
 * ------------------------------------------------------------------------- */

class QuadEdge
{

  /* -- public instance methods -------------------------------------------- */

  public:

  /*
   * Initialize the edges of this quad edge with no connections.
   */
  QuadEdge()
  {
    edges[0].index = 0;
    edges[1].index = 1;
    edges[2].index = 2;
    edges[3].index = 3;

    edges[0].next = edges+0;
    edges[1].next = edges+3;
    edges[2].next = edges+2;
    edges[3].next = edges+1;

    unsigned int id = Edge::nextID;

    edges[0].id = id+0;
    edges[1].id = id+1;
    edges[2].id = id+2;
    edges[3].id = id+3;

    Edge::nextID = id+4;
  }

  /* -- public instance variables ------------------------------------------ */

  public:

  /*
   * The edges of this quad edge.
   */
  Edge edges[4];

};

/* ----------------------------------------------------------------------------
 * Edge
 * ------------------------------------------------------------------------- */

/* -- public class methods ------------------------------------------------- */

Edge *Edge::make()
{
  return (new QuadEdge())->edges;
}

void Edge::kill(Edge *edge)
{
  assert(edge!=0);

  // detach the edge from its cell

  splice(edge, edge->Oprev());
  splice(edge->Sym(), edge->Sym()->Oprev());

  // free the quad edge that the edge belongs to

  delete (QuadEdge*)(edge-edge->index);
}

void Edge::splice(Edge *a, Edge *b)
{
  assert(a!=0);
  assert(b!=0);

  // see Guibas and Stolfi

  Edge* alpha = a->Onext()->Rot();
  Edge* beta  = b->Onext()->Rot();

  Edge* t1 = b->Onext();
  Edge* t2 = a->Onext();
  Edge* t3 = beta->Onext();
  Edge* t4 = alpha->Onext();

  a->next     = t1;
  b->next     = t2;
  alpha->next = t3;
  beta->next  = t4;
}

/* -- public instance methods ---------------------------------------------- */

void Edge::setID(unsigned int id)
{
  assert(id>0);

  this->id = id;
}

void Edge::setOrg(Vertex *org)
{
  // add this edge to the (vertex) orbit of _org_

  vertex = org;

  org->addEdge(this);
}

void Edge::setDest(Vertex *dest)
{
  // add this edge to the (vertex) orbit of _dest_

  Sym()->vertex = dest;

  dest->addEdge(Sym());
}

void Edge::setLeft(Face *left)
{
  // add this edge to the (face) orbit of _left_

  Rot()->face = left;

  left->addEdge(this);
}

void Edge::setRight(Face *right)
{
  // add this edge to the (face) orbit of _right_

  InvRot()->face = right;

  right->addEdge(Sym());
}
// ******************************************************* //
void Edge::setFunction1(unsigned int faceid, Edge* edge){
	//face id cant be 0
	assert(faceid>0);
	// variable to store function value at this edge
	double functionvalue;
	// vertex1 and vertex2 : Origin and Destination of this edge
	Vertex *vertex1;
	Vertex *vertex2;
	vertex1 = edge->Org();
	vertex2 = edge->Dest();
	// getting the projected Y coordinates of the origin and destionation
	// this yk and yk1 are y_k and y_k+1 : origin and destination 
	double yk = vertex1->getProjectedYcoordinate(faceid);
	double yk1 = vertex2->getProjectedYcoordinate(faceid);
	// calculation of the function value 
	functionvalue = yk*yk + yk*yk1+yk1*yk1;
	// assiging function value to the map function1 with key : face id
	this->Function1[faceid] = functionvalue;
}
// ******************************************************* //
void Edge::setFunction2(unsigned int faceid, Edge* edge){
	//face id cant be 0
	assert(faceid>0);
	// variable to store function value at this edge
	double functionvalue;
	// vertex1 and vertex2 : Origin and Destination of this edge
	Vertex *vertex1;
	Vertex *vertex2;
	vertex1 = edge->Org();
	vertex2 = edge->Dest();
	// getting the projected Y coordinates of the origin and destionation
	// this yk and yk1 are y_k and y_k+1 : origin and destination 
	double yk = vertex1->getProjectedYcoordinate(faceid);
	double yk1 = vertex2->getProjectedYcoordinate(faceid);
	double xk = vertex1->getProjectedXcoordinate(faceid);
	double xk1 = vertex2->getProjectedXcoordinate(faceid);
	// calculation of the function value 
	functionvalue = xk*yk1 + 2*xk*yk + 2*xk1*yk1 + xk1*yk;
	// assiging function value to the map function2 with key : face id
	this->Function2[faceid] = functionvalue;
}
// ******************************************************* //
void Edge::setFunction3(unsigned int faceid, Edge* edge){
	//face id cant be 0
	assert(faceid>0);
	// variable to store function value at this edge
	double functionvalue;
	// vertex1 and vertex2 : Origin and Destination of this edge
	Vertex *vertex1;
	Vertex *vertex2;
	vertex1 = edge->Org();
	vertex2 = edge->Dest();
	// getting the projected Y coordinates of the origin and destionation
	// this yk and yk1 are y_k and y_k+1 : origin and destination 
	double xk = vertex1->getProjectedXcoordinate(faceid);
	double xk1 = vertex2->getProjectedXcoordinate(faceid);
	// calculation of the function value 
	functionvalue = xk*xk + xk*xk1+xk1*xk1;
	// assiging function value to the map function3 with key : face id
	this->Function3[faceid] = functionvalue;
}
// ******************************************************* //
void Edge::setAk(unsigned int faceid, Edge* edge){
	//face id cant be 0
	assert(faceid>0);
	// variable to store function value at this edge
	double functionvalue;
	// vertex1 and vertex2 : Origin and Destination of this edge
	Vertex *vertex1;
	Vertex *vertex2;
	vertex1 = edge->Org();
	vertex2 = edge->Dest();
	// getting the projected Y coordinates of the origin and destionation
	// this yk and yk1 are y_k and y_k+1 : origin and destination 
	double yk = vertex1->getProjectedYcoordinate(faceid);
	double yk1 = vertex2->getProjectedYcoordinate(faceid);
	double xk = vertex1->getProjectedXcoordinate(faceid);
	double xk1 = vertex2->getProjectedXcoordinate(faceid);
	// calculation of the function value 
	functionvalue = xk*yk1 - xk1*yk;
	// assiging function value to the map Ak with key : face id
	this->Ak[faceid] = functionvalue;
}
/* -- protected instance methods ------------------------------------------- */

Edge::Edge()
{
  // _index_ is initialized by QuadEdge
  // _next_ is initialized by QuadEdge
  // _id_ is initialized by QuadEdge

  data   = 0;
  vertex = 0;
  face   = 0;
}

Edge::~Edge()
{
}

/* -- private class variables ---------------------------------------------- */

unsigned int Edge::nextID = 4;

