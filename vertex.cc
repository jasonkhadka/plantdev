
/* ============================================================================
 * p2/cell/vertex.cc
 * ========================================================================= */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <svl/Vec3.h>
#include "svl-1.5/include/svl/SVL.h"
#include "cell.hh"
#include "edge.hh"
#include "vertex.hh"
#include "derivative.hh"

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
void Vertex::insertProjectedZcoordinate(unsigned int faceid, double zcood){
  assert(faceid>0);
  this->ProjectedZcoordinate[faceid] = zcood;
}
///////////////////////////////////////////////////////////////////////
double Vertex::getProjectedXcoordinate(unsigned int faceid){
	return this->ProjectedXcoordinate[faceid];
}
///////////////////////////////////////////////////////////////////////
double Vertex::getProjectedYcoordinate(unsigned int faceid){
	return this->ProjectedYcoordinate[faceid];
}
double Vertex::getProjectedZcoordinate(unsigned int faceid){
  return this->ProjectedZcoordinate[faceid];
}
void Vertex::insertNonCentralisedProjectedXcoordinate(unsigned int faceid, double xcood){
  assert(faceid>0);
  this->NonCentralisedProjectedXcoordinate[faceid] = xcood;
}
///////////////////////////////////////////////////////////////////////
void Vertex::insertNonCentralisedProjectedYcoordinate(unsigned int faceid, double ycood){
  assert(faceid>0);
  this->NonCentralisedProjectedYcoordinate[faceid] = ycood;
}
void Vertex::insertNonCentralisedProjectedZcoordinate(unsigned int faceid, double zcood){
  assert(faceid>0);
  this->NonCentralisedProjectedZcoordinate[faceid] = zcood;
}
/////////////
///////////////////////////////////////////////////////////////////////
double Vertex::getNonCentralisedProjectedXcoordinate(unsigned int faceid){
  return this->NonCentralisedProjectedXcoordinate[faceid];
}
///////////////////////////////////////////////////////////////////////
double Vertex::getNonCentralisedProjectedYcoordinate(unsigned int faceid){
  return this->NonCentralisedProjectedYcoordinate[faceid];
}

///////////////////////////////////////////////////////////////////////
double Vertex::getNonCentralisedProjectedZcoordinate(unsigned int faceid){
  return this->NonCentralisedProjectedZcoordinate[faceid];
}

///////////////////////////////////////////////////////////////////////
double * Vertex::getCartesianForce(){
    double * pnt = this->cartesianForce;
    return pnt;
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
  this->Xcoordinate = 0.0;
  this->Ycoordinate = 0.0;
  this->Zcoordinate = 0.0;
  cell->addVertex(this);
  this->domePosition = true;//seting position to dome as True in default
}

Vertex::~Vertex()
{
  cell->removeVertex(this);
}
// ******************************************************* //
void Vertex::setAk(){
  VertexEdgeIterator vertEdges(this);//iterator to iterate through the vertex for outgoign edge
  double functionvalue, xk, xk1, yk, yk1; //variable to store function value
    Edge *currentEdge;//pointer to keep track of a edge that is iterated 
    while ((currentEdge = vertEdges.next())!=0){//iterating till exhausted
        Face *left = currentEdge->Left();//grabbing the inner face of the edge
        unsigned int innerid = left->getID();//id of inner face
        //Edge *prevEdge = currentEdge->Lprev();//edge that is previous on the same left face
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
//***************************************************************************************//
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
//***************************************************************************************//
double Vertex::getMu1(unsigned int faceid){
    VertexEdgeIterator vertEdges(this);//iterator to iterate through the vertex for outgoing edge
    Edge *currentEdge;
    Face *left;//inner face of the current edge
    unsigned int innerid;//id of inner face
    while ((currentEdge = vertEdges.next())!=0){
        left = currentEdge->Left();//grabing the inner face of the edge
        innerid = left->getID();//innerid of the face
        if (innerid == faceid){
            return left->getMu1();
        }
    }
    return 0.;//if not found returns 0.
}
//***************************************************************************************//
double Vertex::getMu2(unsigned int faceid){
    VertexEdgeIterator vertEdges(this);//iterator to iterate through the vertex for outgoing edge
    Edge *currentEdge;
    Face *left;//inner face of the current edge
    unsigned int innerid;//id of inner face
    while ((currentEdge = vertEdges.next())!=0){
        left = currentEdge->Left();//grabing the inner face of the edge
        innerid = left->getID();//innerid of the face
        if (innerid == faceid){
            return left->getMu2();
        }
    }
  return 0.;//if not found returns 0.
}
//***************************************************************************************//
double Vertex::getMu3(unsigned int faceid){
    VertexEdgeIterator vertEdges(this);//iterator to iterate through the vertex for outgoing edge
    Edge *currentEdge;
    Face *left;//inner face of the current edge
    unsigned int innerid;//id of inner face
    while ((currentEdge = vertEdges.next())!=0){
        left = currentEdge->Left();//grabing the inner face of the edge
        innerid = left->getID();//innerid of the face
        if (innerid == faceid){
            return left->getMu3();
        }
    }
    return 0.;//if not found returns 0.

}
//***************************************************************************************//
double Vertex::getMu4(unsigned int faceid){
    VertexEdgeIterator vertEdges(this);//iterator to iterate through the vertex for outgoing edge
    Edge *currentEdge;
    Face *left;//inner face of the current edge
    unsigned int innerid;//id of inner face
    while ((currentEdge = vertEdges.next())!=0){
        left = currentEdge->Left();//grabing the inner face of the edge
        innerid = left->getID();//innerid of the face
        if (innerid == faceid){
            return left->getMu4();
        }
    }
    return 0.;//if not found returns 0.
}
//***************************************************************************************//
void Vertex::setMuXDerivative(){
    VertexEdgeIterator vertEdges(this);//iterator to iterate through the vertex for outgoign edge
    double temp1, temp3; //variable to store value of mu1squared, m2squared, mu4 sqaured derivative
    double xk_1,yk_1,xk,yk,xk1,yk1; // varibales to store the coordinates of the vertices
    Edge *currentEdge;//pointer to keep track of a edge that is iterated 
    Edge *prevEdge;//edge that is previous on the same left face
    Face *left;//inner face (left) of the currentEdge
    unsigned int innerid; //id of inner face
    while ((currentEdge = vertEdges.next())!=0){//iterating till exhausted
        left = currentEdge->Left();//grabbing the inner face of the edge
        innerid = left->getID();//id of inner face
        prevEdge = currentEdge->Lprev();//edge that is previous on the same left face, ccw previous
        //grabbing j+1 and j-1 vertex
        Vertex *vertex1 = prevEdge->Org();//j-1
        Vertex *vertex2 = currentEdge->Org();//j
        Vertex *vertex3 = currentEdge->Dest();//j+1
        //getting the vertices of the above vertex
        //j-1 vertex
        xk_1 = vertex1->getProjectedXcoordinate(innerid);
        yk_1 = vertex1->getProjectedYcoordinate(innerid);
        //j vertex
        xk = vertex2->getProjectedXcoordinate(innerid);
        yk = vertex2->getProjectedYcoordinate(innerid);
        //j+1 vertex
        xk1 = vertex3->getProjectedXcoordinate(innerid);
        yk1 = vertex3->getProjectedYcoordinate(innerid);
        //calculating xderivative of mu1squared, mu2sq, mu4sq
        temp1 = 1./12.*((vertex2->getFunction1(innerid))*yk1-(vertex1->getFunction1(innerid))*yk_1);
        temp3 =1./12.*((vertex2->getFunction3(innerid))*yk1-(vertex1->getFunction3(innerid))*yk_1 +
                (vertex2->getAk(innerid))*2*xk+(vertex1->getAk(innerid))*xk_1 + 
                (vertex2->getAk(innerid))*xk1+(vertex1->getAk(innerid))*2*xk);
        this->mu1XDerivative[innerid] = temp1;
        this->mu4XDerivative[innerid] = temp3;
    }
}
//***************************************************************************************//
void Vertex::setMuYDerivative(){
    VertexEdgeIterator vertEdges(this);//iterator to iterate through the vertex for outgoign edge
    double temp1, temp3; //variable to store value of mu1squared, m2squared, mu4 sqaured derivative
    double xk_1,yk_1,xk,yk,xk1,yk1; // varibales to store the coordinates of the vertices
    Edge *currentEdge;//pointer to keep track of a edge that is iterated 
    Edge *prevEdge;//edge that is previous on the same left face
    Face *left;//inner face (left) of the currentEdge
    unsigned int innerid; //id of inner face
    while ((currentEdge = vertEdges.next())!=0){//iterating till exhausted
        left = currentEdge->Left();//grabbing the inner face of the edge
        innerid = left->getID();//id of inner face
        prevEdge = currentEdge->Lprev();//edge that is previous on the same left face, ccw previous
        //grabbing j+1 and j-1 vertex
        Vertex *vertex1 = prevEdge->Org();//j-1
        Vertex *vertex2 = currentEdge->Org();//j
        Vertex *vertex3 = currentEdge->Dest();//j+1
        //getting the vertices of the above vertex
        //j-1 vertex
        xk_1 = vertex1->getProjectedXcoordinate(innerid);
        yk_1 = vertex1->getProjectedYcoordinate(innerid);
        //j vertex
        xk = vertex2->getProjectedXcoordinate(innerid);
        yk = vertex2->getProjectedYcoordinate(innerid);
        //j+1 vertex
        xk1 = vertex3->getProjectedXcoordinate(innerid);
        yk1 = vertex3->getProjectedYcoordinate(innerid);
        //calculating xderivative of mu1squared, mu2sq, mu4sq
        temp1 = 1./12.*((vertex1->getFunction1(innerid))*xk_1-(vertex2->getFunction1(innerid))*xk1 +
                (vertex2->getAk(innerid))*2*yk+(vertex1->getAk(innerid))*yk_1 + 
                (vertex2->getAk(innerid))*yk1+(vertex1->getAk(innerid))*2*yk);
        temp3 = 1./12.*((vertex1->getFunction3(innerid))*xk_1-(vertex2->getFunction3(innerid))*xk);
        this->mu1YDerivative[innerid] = temp1;
        this->mu4YDerivative[innerid] = temp3;
    }
}
//***************************************************************************************//
void Vertex::setMuSquaredXDerivative(){
    VertexEdgeIterator vertEdges(this);//iterator to iterate through the vertex for outgoign edge
    double temp1, temp2, temp3; //variable to store value of mu1squared, m2squared, mu4 sqaured derivative
    double xk_1,yk_1,xk,yk,xk1,yk1; // varibales to store the coordinates of the vertices
    Edge *currentEdge;//pointer to keep track of a edge that is iterated 
    Edge *prevEdge;//edge that is previous on the same left face
    Face *left;//inner face (left) of the currentEdge
    unsigned int innerid; //id of inner face
    while ((currentEdge = vertEdges.next())!=0){//iterating till exhausted
        left = currentEdge->Left();//grabbing the inner face of the edge
        innerid = left->getID();//id of inner face
        prevEdge = currentEdge->Lprev();//edge that is previous on the same left face, ccw previous
        //grabbing j+1 and j-1 vertex
        Vertex *vertex1 = prevEdge->Org();//j-1
        Vertex *vertex2 = currentEdge->Org();//j
        Vertex *vertex3 = currentEdge->Dest();//j+1
        //getting the vertices of the above vertex
        //j-1 vertex
        xk_1 = vertex1->getProjectedXcoordinate(innerid);
        yk_1 = vertex1->getProjectedYcoordinate(innerid);
        //j vertex
        xk = vertex2->getProjectedXcoordinate(innerid);
        yk = vertex2->getProjectedYcoordinate(innerid);
        //j+1 vertex
        xk1 = vertex3->getProjectedXcoordinate(innerid);
        yk1 = vertex3->getProjectedYcoordinate(innerid);
        //calculating xderivative of mu1squared, mu2sq, mu4sq
        temp1 = (left->getMu1())*1./6.*((vertex2->getFunction1(innerid))*yk1-(vertex1->getFunction1(innerid))*yk_1);
        temp2 = (left->getMu2())*1./12.*((vertex2->getFunction2(innerid))*yk1-(vertex1->getFunction2(innerid))*yk_1 +
                (vertex2->getAk(innerid))*yk1+(vertex2->getAk(innerid))*2*yk + 
                (vertex1->getAk(innerid))*yk+(vertex1->getAk(innerid))*yk_1);
        temp3 =(left->getMu4())*1./6.*((vertex2->getFunction3(innerid))*yk1-(vertex1->getFunction3(innerid))*yk_1 +
                (vertex2->getAk(innerid))*2*xk+(vertex1->getAk(innerid))*xk_1 + 
                (vertex2->getAk(innerid))*xk1+(vertex1->getAk(innerid))*2*xk);
        this->mu1SquaredXDerivative[innerid] = temp1;
        this->mu2SquaredXDerivative[innerid] = temp2;
        this->mu3SquaredXDerivative[innerid] = temp2;
        this->mu4SquaredXDerivative[innerid] = temp3;
    }
}

void Vertex::setMuSquaredYDerivative(){
    VertexEdgeIterator vertEdges(this);//iterator to iterate through the vertex for outgoign edge
    double temp1, temp2, temp3; //variable to store value of mu1squared, m2squared, mu4 sqaured derivative
    double xk_1,yk_1,xk,yk,xk1,yk1; // varibales to store the coordinates of the vertices
    Edge *currentEdge;//pointer to keep track of a edge that is iterated 
    Edge *prevEdge;//edge that is previous on the same left face
    Face *left;//inner face (left) of the currentEdge
    unsigned int innerid; //id of inner face
    while ((currentEdge = vertEdges.next())!=0){//iterating till exhausted
        left = currentEdge->Left();//grabbing the inner face of the edge
        innerid = left->getID();//id of inner face
        prevEdge = currentEdge->Lprev();//edge that is previous on the same left face, ccw previous
        //grabbing j+1 and j-1 vertex
        Vertex *vertex1 = prevEdge->Org();//j-1
        Vertex *vertex2 = currentEdge->Org();//j
        Vertex *vertex3 = currentEdge->Dest();//j+1
        //getting the vertices of the above vertex
        //j-1 vertex
        xk_1 = vertex1->getProjectedXcoordinate(innerid);
        yk_1 = vertex1->getProjectedYcoordinate(innerid);
        //j vertex
        xk = vertex2->getProjectedXcoordinate(innerid);
        yk = vertex2->getProjectedYcoordinate(innerid);
        //j+1 vertex
        xk1 = vertex3->getProjectedXcoordinate(innerid);
        yk1 = vertex3->getProjectedYcoordinate(innerid);
        //calculating xderivative of mu1squared, mu2sq, mu4sq
        temp1 = (left->getMu1())*1./6.*((vertex1->getFunction1(innerid))*xk_1-(vertex2->getFunction1(innerid))*xk1 +
                (vertex2->getAk(innerid))*2*yk+(vertex1->getAk(innerid))*yk_1 + 
                (vertex2->getAk(innerid))*yk1+(vertex1->getAk(innerid))*2*yk);
        temp2 = (left->getMu2())*1./12.*((vertex1->getFunction2(innerid))*xk_1-(vertex2->getFunction2(innerid))*xk1 +
                (vertex1->getAk(innerid))*xk_1+(vertex2->getAk(innerid))*2*xk + 
                (vertex1->getAk(innerid))*2*xk+(vertex2->getAk(innerid))*xk1);
        temp3 =(left->getMu4())*1./6.*((vertex1->getFunction3(innerid))*xk_1-(vertex2->getFunction3(innerid))*xk);
        this->mu1SquaredYDerivative[innerid] = temp1;
        this->mu2SquaredYDerivative[innerid] = temp2;
        this->mu3SquaredYDerivative[innerid] = temp2;
        this->mu4SquaredYDerivative[innerid] = temp3;
    }
}
// ******************************************************* //
void Vertex::setFunctions(){
  this->setFunction1();
  this->setFunction2();
  this->setFunction3();
}
// ******************************************************* //
void Vertex::setFunction1(){
  VertexEdgeIterator vertEdges(this);//iterator to iterate through the vertex for outgoign edge
  double functionvalue, yk, yk1; //variable to store function value
    Edge *currentEdge;//pointer to keep track of a edge that is iterated 
    while ((currentEdge = vertEdges.next())!=0){//iterating till exhausted
        Face *left = currentEdge->Left();//grabbing the inner face of the edge
        unsigned int innerid = left->getID();//id of inner face
        //Edge *prevEdge = currentEdge->Lprev();//edge that is previous on the same left face
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
        //Edge *prevEdge = currentEdge->Lprev();//edge that is previous on the same left face
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
  double functionvalue(0), xk, xk1; //variable to store function value
    Edge *currentEdge;//pointer to keep track of a edge that is iterated 
    while ((currentEdge = vertEdges.next())!=0){//iterating till exhausted
        Face *left = currentEdge->Left();//grabbing the inner face of the edge
        unsigned int innerid = left->getID();//id of inner face
        //Edge *prevEdge = currentEdge->Lprev();//edge that is previous on the same left face
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
// %%%%%%%%%%%%%%%%%%%%%%%%% X Derivative %%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
// ******************************************************* //
void Vertex::setFirstTermXDerivative(){
    VertexEdgeIterator vertEdges(this);//iterator to iterate through the outgoing edge
    double functionvalue(0), area(0), sumMuSqDerivative(0),sumMuSq(0),areaDerivative(0);//variable to store function value
    unsigned int innerid; //variable to store inner id of face
    Edge *currentEdge;//pointer to the current edge that is iterated
    while ((currentEdge = vertEdges.next())!=0){//iterating till exhausted
        Face *left = currentEdge->Left();//grabbing the inner face of the edge
        innerid = left->getID();//id of inner face
        if (innerid != 1){//skipping the rest of the calculation, if the face is outer face
            area = left->getAreaOfFace();//area of inner face
            sumMuSqDerivative = this->getMu1SquaredXDerivative(innerid) + this->getMu2SquaredXDerivative(innerid) +
                                 this->getMu3SquaredXDerivative(innerid) + this->getMu4SquaredXDerivative(innerid);//sumof mu sq X derivaitves
            sumMuSq = pow(left->getMu1(),2)+ pow(left->getMu2(),2)+pow(left->getMu3(),2)+ pow(left->getMu4(),2);
            areaDerivative = this->getAreaXDerivative(innerid);
            functionvalue += 1./(left->getTraceSquaredTargetFormMatrix())*(area*(sumMuSqDerivative)+sumMuSq*areaDerivative);
            //std::cout<<"vertid : "<<this->getID()<<"face id : "<<innerid<<"functionvalue :"<<functionvalue<<std::endl;
        }
    }
    this->firstTermXDerivative = functionvalue;
    //std::cout<<"written !! "<<this->firstTermXDerivative<<std::endl;
    //std::cout<<"vertid : "<<this->getID()<<"  face id : "<<innerid<<"  functionvalue :"<<functionvalue<<std::endl;
}

// ******************************************************* //
 void Vertex::setSecondTermXDerivative(){
    VertexEdgeIterator vertEdges(this);//iterator to iterate through the outgoing edge
    double functionvalue(0), area, sumMuDerivative,sumMu,areaDerivative;//variable to store function value
    unsigned int innerid; //variable to store inner id of face
    Edge *currentEdge;//pointer to the current edge that is iterated
    while ((currentEdge = vertEdges.next())!=0){//iterating till exhausted
        Face *left = currentEdge->Left();//grabbing the inner face of the edge
        innerid = left->getID();//id of inner 
        if (innerid == 1) continue;
        area = left->getAreaOfFace();//area of inner face
        sumMuDerivative = (left->getMu4())*(this->getMu1XDerivative(innerid)) + 
                          (left->getMu1())*(this->getMu4XDerivative(innerid));//sumof mu sq X derivaitves
        sumMu = (left->getMu1()+left->getMu4());
        areaDerivative = this->getAreaXDerivative(innerid);
        functionvalue += 1./(left->getTraceSquaredTargetFormMatrix())*
                          ((sumMu*sumMu)*areaDerivative+area*2*sumMu*sumMuDerivative);
    }
    this->secondTermXDerivative = functionvalue;
}
// ******************************************************* //
void Vertex::setThirdTermXDerivative(){
    VertexEdgeIterator vertEdges(this);//iterator to iterate through the outgoing edge
    Cell *currentCell = this->getCell();// cell this vertex lies in
    Face *left;//pointer to the left face of the edge
    //double pressure = currentCell->getPressure();//pressure in the cells of this tissue
    double functionvalue(0), areaDerivative;//variable to store function value
    unsigned int innerid; //variable to store inner id of face
    Edge *currentEdge;//pointer to the current edge that is iterated
    //Edge *prevEdge; // pointer to edge previous on the left face of the current edge
    while ((currentEdge = vertEdges.next())!=0){//iterating till exhausted
        left = currentEdge->Left();//grabbing the inner face of the edge
        innerid = left->getID();//id of inner face
        //getting area derivative wrt X
        areaDerivative = this->getAreaXDerivative(innerid);
        //calculating xderivative third term
        functionvalue += areaDerivative;
    }
    this->thirdTermXDerivative =functionvalue;
}
// %%%%%%%%%%%%%%%%%%%%%%%%% Y Derivative %%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
void Vertex::setFirstTermYDerivative(){
    VertexEdgeIterator vertEdges(this);//iterator to iterate through the outgoing edge
    double functionvalue(0), area, sumMuSqDerivative,sumMuSq,areaDerivative;//variable to store function value
    unsigned int innerid; //variable to store inner id of face
    Edge *currentEdge;//pointer to the current edge that is iterated
    while ((currentEdge = vertEdges.next())!=0){//iterating till exhausted
        Face *left = currentEdge->Left();//grabbing the inner face of the edge
        innerid = left->getID();//id of inner face
        if (innerid == 1) continue; //skipping if face is the outer face
        area = left->getAreaOfFace();//area of inner face
        sumMuSqDerivative = this->getMu1SquaredYDerivative(innerid) + this->getMu2SquaredYDerivative(innerid) +
                             this->getMu3SquaredYDerivative(innerid) + this->getMu4SquaredYDerivative(innerid);//sumof mu sq X derivaitves
        sumMuSq = (left->getMu1())*(left->getMu1()) + (left->getMu2())*(left->getMu2()) +
                  (left->getMu3())*(left->getMu3()) + (left->getMu4())*(left->getMu4());
        areaDerivative = this->getAreaYDerivative(innerid);
        functionvalue += 1./(left->getTraceSquaredTargetFormMatrix())*(area*(sumMuSqDerivative)+sumMuSq*areaDerivative);
    }
    this->firstTermYDerivative = functionvalue;
}
// ******************************************************* //
 void Vertex::setSecondTermYDerivative(){
    VertexEdgeIterator vertEdges(this);//iterator to iterate through the outgoing edge
    double functionvalue(0), area, sumMuDerivative,sumMu,areaDerivative;//variable to store function value
    unsigned int innerid; //variable to store inner id of face
    Edge *currentEdge;//pointer to the current edge that is iterated
    while ((currentEdge = vertEdges.next())!=0){//iterating till exhausted
        Face *left = currentEdge->Left();//grabbing the inner face of the edge
        innerid = left->getID();//id of inner face
        if (innerid == 1) continue; //skipping if face is the outer face
        area = left->getAreaOfFace();//area of inner face
        sumMuDerivative = (left->getMu4())*(this->getMu1YDerivative(innerid)) + 
                          (left->getMu1())*(this->getMu4YDerivative(innerid));//sumof mu sq X derivaitves
        sumMu = (left->getMu1()+left->getMu4());
        areaDerivative = this->getAreaYDerivative(innerid);
        functionvalue += 1./(left->getTraceSquaredTargetFormMatrix())*
                          ((sumMu*sumMu)*areaDerivative+area*2*sumMu*sumMuDerivative);
    }
    this->secondTermYDerivative = functionvalue;
}
// ******************************************************* //
void Vertex::setThirdTermYDerivative(){
    VertexEdgeIterator vertEdges(this);//iterator to iterate through the outgoing edge
    Cell *currentCell = this->getCell();// cell this vertex lies in
    Face *left;//pointer to the left face of the edge
    //double pressure = currentCell->getPressure();//pressure in the cells of this tissue
    double functionvalue(0), areaDerivative;//variable to store function value
    unsigned int innerid; //variable to store inner id of face
    Edge *currentEdge;//pointer to the current edge that is iterated
    //Edge *prevEdge; // pointer to edge previous on the left face of the current edge
    while ((currentEdge = vertEdges.next())!=0){//iterating till exhausted
        left = currentEdge->Left();//grabbing the inner face of the edge
        innerid = left->getID();//id of inner face
        //getting area derivative wrt X
        areaDerivative = this->getAreaYDerivative(innerid);
        //calculating xderivative third term
        functionvalue += areaDerivative;
    }
    this->thirdTermYDerivative =functionvalue;
}
// ******************************************************* //
void Vertex::setAlphaBetaGamma(){
    VertexEdgeIterator vertEdges(this);//iterator to iterate through the outgoing edge
    Face *left; //pointer to left face
    double xmean, ymean, zmean, tempalpha, tempgamma, tempbeta;
    double x1, y1, z1, x2, y2, z2;
    unsigned int innerid;// ID of inner Face
    Edge *currentEdge; // pointer to current edge
    Vertex *vertex1, *vertex2;
    while ((currentEdge = vertEdges.next())!=0){
        left = currentEdge->Left();//grabbing the inner face of the edge
        innerid = left->getID();
        vertex1 = currentEdge->Org();
        vertex2 = currentEdge->Dest();
        x1 = vertex1->getXcoordinate();
        y1 = vertex1->getYcoordinate();
        z1 = vertex1->getZcoordinate();
        x2 = vertex2->getXcoordinate();
        y2 = vertex2->getYcoordinate();
        z2 = vertex2->getZcoordinate();
        // getting the mean x, y and z from the left 
        xmean = left->getXCentralised();
        ymean = left->getYCentralised();
        zmean = left->getZCentralised();
        //calculating alpha beta and gamma
        tempalpha = (y1-ymean)*(z2-zmean)-(z1-zmean)*(y2-ymean);
        tempbeta = (z1-zmean)*(x2-xmean)-(x1-xmean)*(z2-zmean);
        tempgamma = (x1-xmean)*(y2-ymean)-(y1-ymean)*(x2-xmean);
        //now storing the values
        this->alpha[innerid] = tempalpha;
        this->beta[innerid] = tempbeta;
        this->gamma[innerid] = tempgamma;
    }
}
// ******************************************************* //
// ******************************************************* //
void Vertex::calculateCartesianForce(){
    // getting a derivative instance
    CentralisedDerivative derivative;
    //calculating the derivative now
    this->cartesianForce[0] = derivative.numericalEnergyXDerivative(this);
    this->cartesianForce[1] = derivative.numericalEnergyYDerivative(this);
    this->cartesianForce[2] = derivative.numericalEnergyZDerivative(this);
}
// ******************************************************* //
void Vertex::setCylindrical(){
    this->rCoordinate = sqrt(pow(this->Xcoordinate,2)+pow(this->Ycoordinate,2));
    this->thetaCoordinate = atan2((this->Ycoordinate),(this->Xcoordinate));//arctan2 is needed for this conversion
}

// ******************************************************* //
void Vertex::setCartesian(){
    this->Xcoordinate = (this->rCoordinate)*cos(this->thetaCoordinate);
    this->Ycoordinate = (this->rCoordinate)*sin(this->thetaCoordinate);
}
// ******************************************************* //
void Vertex::setparameters(){
    //setting alpha beta gamma
    this->setAlphaBetaGamma();
    //Here on All need Projected coordinates <<>>
    //setting function 1, 2 & 3
    this->setFunctions();
    //set A_k
    this->setAk();
}
// ******************************************************** //
void Vertex::setDerivatives(){
    //setting A_k derivative
    this->setAkDerivative();
    //setting area derivative
    this->setAreaDerivative();
    //setting Mux derivative
    this->setMuXDerivative();
    //setting Muy derivative
    this->setMuYDerivative();
    //setting Musq X derivative
    this->setMuSquaredXDerivative();
    //setting Musq Y derivative
    this->setMuSquaredYDerivative();
    //setting first x derivative 
    this->setFirstTermXDerivative();
    //setting second x deriv
    this->setSecondTermXDerivative();
    //setting third xderiv
    this->setThirdTermXDerivative();
    //setting first Y derivative 
    this->setFirstTermYDerivative();
    //setting second Y deriv
    this->setSecondTermYDerivative();
    //setting third Y deriv
    this->setThirdTermYDerivative();
}



