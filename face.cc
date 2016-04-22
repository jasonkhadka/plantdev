
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

void Face::addEdge(Edge *edge)
{
  assert(edge!=0);

  // only keep track of one edge in the orbit--this one is as good as any

  this->edge = edge;
}

void Face::removeEdge(Edge *edge)
{
  assert(edge!=0);

  // replace the arbitrary edge with another edge in the orbit
  // use null if this is the only edge
  // assumes that the edge hasn't been actually removed yet

  Edge *next = edge->Onext();

  this->edge = next!=edge ? next : 0;
}

double Face::getAreaOfFace(){
	return this->areaOfFace;
}

void Face::setAreaOfFace(){
	FaceEdgeIterator edges(this);//iterator to iterate through the face
	Edge *newedge;//a pointer to keep track of current edge
	double areasum;//variable to store the area of the face
	unsigned int faceid = this->id;//getting the faceid
	//looping through all the edges in the face until exhausted
	while ((newedge = edges.next())!=0){
		areasum += newedge->Ak[faceid];//summing up all the Ak values for this face
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

