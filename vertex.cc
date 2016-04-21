
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

