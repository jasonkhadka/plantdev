/* ============================================================================
 * p2/cell/cell.cc
 * ========================================================================= */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h> //has abs and others
#include <math.h> //sqrt and others
#include <vector>// for vector 
#include <sys/time.h>//to get time of the day

//to get Maximum double value 
#include <limits> //std::numeric_limits

//random number generating
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>

// NLOPT for relaxation
#include <nlopt.hpp>

//QuadEdge libraries
#include "cell.hh"
#include "face.hh"
#include "vertex.hh"

//Macro to change 2d index to linear index [i.j] = i*numOfRow + j
#define getIndex(i,j) (i*3 + j)
// ***************************************************************** //
// This is to seed the random number generator, for "truly" getting random number
unsigned long int random_seed()
{

 unsigned int seed;
 struct timeval tv;
 FILE *devrandom;

 if ((devrandom = fopen("/dev/urandom","r")) == NULL) {
   gettimeofday(&tv,0);
   seed = tv.tv_sec + tv.tv_usec;
   //printf("Got seed %u from gettimeofday()\n",seed);
 } else {
   fread(&seed,sizeof(seed),1,devrandom);
   //printf("Got seed %u from /dev/random\n",seed);
   fclose(devrandom);
 }

 return(seed);

}
// ********************************************************************** //
/* ----------------------------------------------------------------------------
 * Cell
 * ------------------------------------------------------------------------- */
// added features //
 void Cell::setPressure(double inputpressure){
      this->pressure = inputpressure;
 }
 void Cell::setAlpha(double inputAlpha){
      this->alpha = inputAlpha;
 }
 void Cell::setBeta(double inputBeta){
      this->beta = inputBeta;
 }
 /**
  * Energy of cell calculator,
  * bascially go over faces and sum up the energies
  */
 double Cell::getEnergy(){
      //iterating the faces
      CellFaceIterator faces(this);
      Face * face;
      double totalenergy = 0.;
      while ((face = faces.next())!= 0){
          if (face->getID()==1) continue;
          totalenergy += face->getEnergy();
      }
      totalenergy -= this->getGamma()*this->getVolume();//subtracting the volume energy
      totalenergy +=  this->getZeta()*this->getFourthTerm();//subtracting the fourth term : z_proj penalty
      return totalenergy ;
 }
//******************************************************************************* //
  /**
  * Energy of cell calculator,
  * bascially go over faces and sum up the energies with VOlume -> cartesian volume without centroid
  */
 double Cell::getEnergyCartesianVolumeWOCentroid(){
      //iterating the faces
      CellFaceIterator faces(this);
      Face * face;
      double totalenergy = 0.;
      while ((face = faces.next())!= 0){
          if (face->getID()==1) continue;
          totalenergy += face->getEnergy();
      }
      totalenergy -= this->getGamma()*this->getCartesianVolumeWOCentroid();
      totalenergy +=  this->getZeta()*this->getFourthTerm();//subtracting the fourth term : z_proj penalty
      return totalenergy ;
 }

//******************************************************************************* //
  /**
  * Energy of cell calculator,
  * bascially go over faces and sum up the energies ->
  * Volume : cartesian volume with centroid
  */
 double Cell::getEnergyCartesianVolume(){
      //iterating the faces
      // First : check if the cell is convex or not
      if (!(this->isConvex())){ //if the cell is not convex then return HIGHVALUE
          std::cout<<"Convesity Reached ! : "<<this->isConvex()<<std::endl;
          return std::numeric_limits<double>::max();
      }
      double fourthterm = this->getFourthTerm();
      // Second : check if the cell has bent than threshold
      if (fourthterm > (this->getBendingThreshold())){
        std::cout<<"Bending Threshold Reached"<<std::endl;
        return std::numeric_limits<double>::max();
      }
      // If not continue with calculation of energy
      CellFaceIterator faces(this);
      Face * face;
      double totalenergy = 0.;
      while ((face = faces.next())!= 0){
          if (face->getID()==1) continue;
          totalenergy += face->getEnergy();
      }
      totalenergy -= this->getGamma()*this->getCartesianVolume();
      totalenergy +=  this->getZeta()*fourthterm;//subtracting the fourth term : z_proj penalty

      return totalenergy ;
 }
 //******************************************************************************* //
   /**
  * Energy of cell calculator,
  * bascially go over faces and sum up the energies ->
  * Volume : cartesian volume with centroid
  */
 double Cell::getDerivativeEnergyCartesianVolume(){
      //iterating the faces
      // First : check if the cell is convex or not
      /*
      if (!(this->isConvex())){ //if the cell is not convex then return HIGHVALUE
          std::cout<<"Convesity Reached ! : "<<this->isConvex()<<std::endl;
          return std::numeric_limits<double>::max();
      }
      */
      double fourthterm = this->getFourthTerm();
      /*
      // Second : check if the cell has bent than threshold
      if (fourthterm > (this->getBendingThreshold())){
        std::cout<<"Bending Threshold Reached"<<std::endl;
        return std::numeric_limits<double>::max();
      }
      */
      // If not continue with calculation of energy
      CellFaceIterator faces(this);
      Face * face;
      double totalenergy = 0.;
      while ((face = faces.next())!= 0){
          if (face->getID()==1) continue;
          totalenergy += face->getEnergy();
      }
      totalenergy -= this->getGamma()*this->getCartesianVolume();
      totalenergy +=  this->getZeta()*fourthterm;//subtracting the fourth term : z_proj penalty

      return totalenergy ;
 }
//******************************************************************************* //
 /**
  * caculating the volume without Centroid
  */
 double Cell::getCartesianVolumeWOCentroid(){
    //iterating the faces
  CellFaceIterator faces(this);
  Face * face;
  Edge * edge;
  Vertex * first;
  Vertex * second;
  Vertex * third;
  double tempVolume = 0.;
  double productVector[3];
  while ((face= faces.next()) != 0){
    if (face->getID()==1) continue;
    //getting the centroid of face
    FaceEdgeIterator edges(face);
    edge = edges.next();
    first = edge->Org();//grabbing the first vertex in the face
    double firstCoordinate[3] = {first->getXcoordinate(),first->getYcoordinate(),first->getZcoordinate()};
    //double centroid[3] = {face->getXCentralised(),face->getYCentralised(),face->getZCentralised()};
    // iterating the faces
    while ((edge = edges.next()) != 0){
          second = edge->Org();
          third = edge->Dest();
          double secondCoordinate[3] = {second->getXcoordinate(),second->getYcoordinate(),second->getZcoordinate()};
          double thirdCoordinate[3] = {third->getXcoordinate(),third->getYcoordinate(),third->getZcoordinate()};
          //calculating cross product of two vectors first and second
          productVector[0] = secondCoordinate[1]*thirdCoordinate[2] - secondCoordinate[2]*thirdCoordinate[1];//also definition of alpha
          productVector[1] = secondCoordinate[2]*thirdCoordinate[0] - secondCoordinate[0]*thirdCoordinate[2];//also definition of beta
          productVector[2] = secondCoordinate[0]*thirdCoordinate[1] - secondCoordinate[1]*thirdCoordinate[0];//also definition of gamma
          //calculating dot product and resulting volumne
          //std::cout<< " centroid : "<< centroid[0] <<centroid[1]<<centroid[2]<<std::endl;
          tempVolume += (1./6.)*abs(productVector[0]*firstCoordinate[0]+
                                    productVector[1]*firstCoordinate[1]+
                                    productVector[2]*firstCoordinate[2]);
        }
        //std::cout<<"Volume " << tempVolume<<std::endl;
  }
  return tempVolume;
 }
  /**
  * caculating the volume
  */
  double Cell::getVolume(){
  //iterating the faces
  CellFaceIterator faces(this);
  Face * face;
  Edge * edge;
  Vertex * second;
  Vertex * third;
  double tempVolume = 0.;
  double productVector[3];
  while ((face= faces.next()) != 0){
    if (face->getID()==1) continue;
    //getting the centroid of face
    FaceEdgeIterator edges(face);
    double centroid[3] = {face->getXCentralised(),face->getYCentralised(),face->getZCentralised()};
    //double centroid[3] = {face->getXCentralised(),face->getYCentralised(),face->getZCentralised()};
    // iterating the faces
    while ((edge = edges.next()) != 0){
          second = edge->Org();
          third = edge->Dest();
          double secondCoordinate[3] = {second->getNonCentralisedProjectedXcoordinate(face->getID()),second->getNonCentralisedProjectedYcoordinate(face->getID()),second->getNonCentralisedProjectedZcoordinate(face->getID())};
          double thirdCoordinate[3] = {third->getNonCentralisedProjectedXcoordinate(face->getID()),third->getNonCentralisedProjectedYcoordinate(face->getID()),third->getNonCentralisedProjectedZcoordinate(face->getID())};
          //calculating cross product of two vectors first and second
          productVector[0] = secondCoordinate[1]*thirdCoordinate[2] - secondCoordinate[2]*thirdCoordinate[1];//also definition of alpha
          productVector[1] = secondCoordinate[2]*thirdCoordinate[0] - secondCoordinate[0]*thirdCoordinate[2];//also definition of beta
          productVector[2] = secondCoordinate[0]*thirdCoordinate[1] - secondCoordinate[1]*thirdCoordinate[0];//also definition of gamma
          //calculating dot product and resulting volumne
          //std::cout<< " centroid : "<< centroid[0] <<centroid[1]<<centroid[2]<<std::endl;
          tempVolume += (1./6.)*abs(productVector[0]*centroid[0]+
                                    productVector[1]*centroid[1]+
                                    productVector[2]*centroid[2]);
        }
        //std::cout<<"Volume " << tempVolume<<std::endl;
  }
  return tempVolume;
 }
 /**
  * caculating the volume
  */
  double Cell::getCartesianVolume(){
  //iterating the faces
  CellFaceIterator faces(this);
  Face * face;
  Edge * edge;
  Vertex * second;
  Vertex * third;
  double tempVolume = 0.;
  double productVector[3];
  while ((face= faces.next()) != 0){
    if (face->getID()==1) continue;
    //getting the centroid of face
    FaceEdgeIterator edges(face);
    double centroid[3] = {face->getXCentralised(),face->getYCentralised(),face->getZCentralised()};
    //double centroid[3] = {face->getXCentralised(),face->getYCentralised(),face->getZCentralised()};
    // iterating the faces
    while ((edge = edges.next()) != 0){
          second = edge->Org();
          third = edge->Dest();
          double secondCoordinate[3] = {second->getXcoordinate(),second->getYcoordinate(),second->getZcoordinate()};
          double thirdCoordinate[3] = {third->getXcoordinate(),third->getYcoordinate(),third->getZcoordinate()};
          //calculating cross product of two vectors first and second
          productVector[0] = secondCoordinate[1]*thirdCoordinate[2] - secondCoordinate[2]*thirdCoordinate[1];//also definition of alpha
          productVector[1] = secondCoordinate[2]*thirdCoordinate[0] - secondCoordinate[0]*thirdCoordinate[2];//also definition of beta
          productVector[2] = secondCoordinate[0]*thirdCoordinate[1] - secondCoordinate[1]*thirdCoordinate[0];//also definition of gamma
          //calculating dot product and resulting volumne
          //std::cout<< " centroid : "<< centroid[0] <<centroid[1]<<centroid[2]<<std::endl;
          tempVolume += (1./6.)*abs(productVector[0]*centroid[0]+
                                    productVector[1]*centroid[1]+
                                    productVector[2]*centroid[2]);
        }
        //std::cout<<"Volume " << tempVolume<<std::endl;
  }
  return tempVolume;
 }
//******************************************************************************* //
 /**
  * First term of Energy of cell calculator,
  * bascially go over faces and sum up the energies
  */
 double Cell::getFirstTerm(){
      //iterating the faces
      CellFaceIterator faces(this);
      Face * face;
      double totalenergy = 0.;
      while ((face = faces.next())!= 0){
          if (face->getID()==1) continue;
          totalenergy += face->getFirstTerm();
      }
      return totalenergy;
 }
 /**
  * Second term of Energy of cell calculator,
  * bascially go over faces and sum up the energies
  */
 double Cell::getSecondTerm(){
      //iterating the faces
      CellFaceIterator faces(this);
      Face * face;
      double totalenergy = 0.;
      while ((face = faces.next())!= 0){
          if (face->getID()==1) continue;
          totalenergy += face->getSecondTerm();
      }
      return totalenergy;
 }
/**
  * Third term of Energy of cell calculator,
  * bascially go over faces and sum up the energies
  */
 double Cell::getThirdTerm(){
      //iterating the faces
      CellFaceIterator faces(this);
      Face * face;
      double totalenergy = 0.;
      while ((face = faces.next())!= 0){
          if (face->getID()==1) continue;
          totalenergy += face->getThirdTerm();
      }
      return totalenergy;
 }
 /**
  * Fourth term of the cell energy: 
  * this calculates the deformation off the plane
  * sum of the z_projected coordinate or the perpendicular distance of vertices from the plane of projection
  */
 double Cell::getFourthTerm(){
    CellFaceIterator faces(this);
    Face * face;
    Edge * edge;
    Vertex * vertex;
    double totalenergy = 0.;
    while((face = faces.next()) != 0){
       if (face->getID() == 1) continue;
      FaceEdgeIterator edges(face);
      while ((edge = edges.next())!= 0){
          vertex = edge->Dest();
          totalenergy += pow(vertex->getProjectedZcoordinate(face->getID()),2);
      }
    }
  return totalenergy;
 }
 //********************************************************************************* //
 void Cell::setCylindrical(){
    CellVertexIterator vertices(this);
    Vertex * vertex;
    while ((vertex = vertices.next())!= 0){
      vertex->setCylindrical();
    }
 }
//********************************************************************************* //
 void Cell::setCartesian(){
    CellVertexIterator vertices(this);
    Vertex * vertex;
    while ((vertex = vertices.next())!= 0){
      vertex->setCartesian();
    }
 }

 //********************************************************************************* //
 void Cell::setInitialParameters(){
  
  Face * face;
  //Edge * edge;
  Vertex * vertex;
  {
    CellFaceIterator faces(this);
    while((face = faces.next())!= 0){
          if (face->getID() == 1) continue;
          face->setProjectedCoordinate();
    }
  }
  /////////////////////////////////////////
  {
    CellVertexIterator vertices(this);
    while ((vertex = vertices.next())!= 0){
      vertex->setparameters();
    }
  }
  /////////////////////////////////////////
  {
    CellFaceIterator faces(this);
    while((face = faces.next())!= 0){
          if (face->getID() == 1) continue;
          face->setMu();
    }
  }
  /////////////////////////////////////////
  {
    CellFaceIterator faces(this);
    while((face = faces.next())!= 0){
          if (face->getID() == 1) continue;
          face->setTempTargetFormMatrixCurrent();
    }
  }
  /////////////////////////////////////////
  {
    CellFaceIterator faces(this);
    while((face = faces.next())!= 0){
          if (face->getID() == 1) continue;
          face->setMu();
    }
  }
  /////////////////////////////////////////
  {
    CellFaceIterator faces(this);
    while((face = faces.next())!= 0){
          if (face->getID() == 1) continue;
          face->setAreaOfFace();
    }
  }
  /////////////////////////////////////////
  {
    CellFaceIterator faces(this);
    while((face = faces.next())!= 0){
          if (face->getID() == 1) continue;
          face->setEnergyTerms();
    }
  }
  /////////////////////////////////////////
  {
    CellFaceIterator faces(this);
    while((face = faces.next())!= 0){
          if (face->getID() == 1) continue;
          face->setDivisionThreshold();
    }
  }
  /////////////////////////////////////////
  // Setting the bendingThreshold to initial bending energy value
  this->setBendingThreshold(this->getFourthTerm());
 }
 //********************************************************************************* //
void Cell::setParameters(){
  Face * face;
  //Edge * edge;
  Vertex * vertex;
  {
    CellFaceIterator faces(this);
    
    while((face = faces.next())!= 0){
          if (face->getID() == 1) continue;
          face->setProjectedCoordinate();
    }
  }
  //////////////////////////////////
  {
    CellVertexIterator vertices(this);
    while ((vertex = vertices.next())!= 0){
      vertex->setparameters();
    } 
  }
  //////////////////////////////////
  {
    CellFaceIterator faces(this);
    while((face = faces.next())!= 0){
          if (face->getID() == 1) continue;
          face->setMu();
    }
  }
  //////////////////////////////////
  {
    CellFaceIterator faces(this);
    while((face = faces.next())!= 0){
          if (face->getID() == 1) continue;
          face->setAreaOfFace();
    }
  }
  //////////////////////////////////
  {
    CellFaceIterator faces(this);
    while((face = faces.next())!= 0){
          if (face->getID() == 1) continue;
          face->setEnergyTerms();
    }
  }
  //////////////////////////////////
 }
//********************************************************************************* //
 void Cell::setDivisionFactor(double newfactor){
  this->divisionFactor = newfactor;
  {
    Face * face;
    CellFaceIterator faces(this);
    while((face = faces.next())!= 0){
          face->setDivisionFactor(newfactor);
    }
  }
}
//********************************************************************************* //
bool Cell::isConvex(){
  CellFaceIterator faces(this);
  Face * face;
  while ((face = faces.next())!=0){
        if (face->getID() == 1) continue;
        if (!(face->isConvex())){// if face->isConvex() condition returns False
            return false;  
        }
      }
  return true; // if non of the face return isConvex() == False, then all are convex, hence cell is Convex
}
 //********************************************************************************* //


 /*------------------------------------------------------------------------- */
/* -- public class methods ------------------------------------------------- */

Cell *Cell::make()
{
  // create a looping edge that connects to itself at a single vertex
  // the edge delimits two faces
  // this is the smallest cell that is consistent with our invariants

  Cell *cell = new Cell();

  Vertex *vertex = Vertex::make(cell);
  Face   *left   = Face::make(cell);
  Face   *right  = Face::make(cell);
  Edge   *edge   = Edge::make()->InvRot();

  edge->setOrg(vertex);
  edge->setDest(vertex);
  edge->setLeft(left);
  edge->setRight(right);
//********************* added features *****************************************//
  //adding vertex to left and right face
  //left->addVertex(vertex);
  //right->addVertex(vertex);
  //initialising the thresholdMatrix
  //cell->thresholdMatrix[0][0] = 0.;
  //cell->thresholdMatrix[0][1] = 0.;
  //cell->thresholdMatrix[1][0] = 0.;
  //cell->thresholdMatrix[1][1] = 0.;
//********************* end added features *****************************************//

  
  return cell;
}

Cell *Cell::makeTetrahedron()
{
  // use the Euler operators to construct a tetrahedron

  Cell *cell = make();

  // grab the initial vertex

  Vertex *vertex1;

  {
    CellVertexIterator iterator(cell);

    vertex1 = iterator.next();

    assert(vertex1!=0);
  }

  // grab the initial edge and the initial faces

  Edge *edge1 = vertex1->getEdge();
  Face *left  = edge1->Left();
  Face *right = edge1->Right();

  // drop in four vertices along the initial edge

  Vertex *vertex2 = cell->makeVertexEdge(vertex1, left, right)->Dest();
  Vertex *vertex3 = cell->makeVertexEdge(vertex2, left, right)->Dest();
  Vertex *vertex4 = cell->makeVertexEdge(vertex3, left, right)->Dest();

  // cut each of the faces in half from complementary vertices

  Face *front  = cell->makeFaceEdge(left, vertex2, vertex4)->Right();
  Face *bottom = cell->makeFaceEdge(right, vertex1, vertex3)->Right();

  return cell;
}

void Cell::kill(Cell *cell)
{
  assert(cell!=0);

  delete cell;
}

/* -- public instance methods (Euler operators) ---------------------------- */

Edge *Cell::makeVertexEdge(Vertex *vertex, Face *left, Face *right)
{
  assert(vertex!=0);
  assert(left!=0);
  assert(right!=0);

  // locate the edges to the right of each of the faces in the orbit of the
  // vertex

  Edge *edge  = vertex->getEdge();
  Edge *edge1 = getOrbitLeft(edge, right);
  Edge *edge2 = getOrbitLeft(edge, left);

  if (edge1==0)
  {
    (void)fprintf(stderr, "Cell::makeVertexEdge: unable to locate right face %u on vertex %u\n",
		 right->getID(), vertex->getID());
    abort();
  }

  if (edge2==0)
  {
    (void)fprintf(stderr, "Cell::makeVertexEdge: unable to locate left face %u on vertex %u\n",
		 left->getID(), vertex->getID());
    abort();
  }

  // create a new vertex and copy the position of the vertex of origin

  Vertex *vertexNew = Vertex::make(this);

  vertexNew->pos = vertex->pos;

  // create a new edge and rotate it to make a clockwise loop

  Edge *edgeNew = Edge::make()->Rot();

  // connect the origin (and destination) of the new edge to _vertex_ so that
  // the left face of the edge is _left_
  // this makes a loop on the inside of _left_

  Edge::splice(edge2, edgeNew);

  // split the origin and destination of the loop so that the right face of the
  // edge is now _right_
  // this results in a non-loop edge dividing _left_ from _right_

  Edge::splice(edge1, edgeNew->Sym());

  // initialize the secondary attributes of the new edge

  edgeNew->setOrg(edge1->Org());
  edgeNew->setLeft(edge2->Left());
  edgeNew->setRight(edge1->Left());

  // all edges leaving the destination orbit of the new edge now have the new
  // vertex as their vertex of origin
//********************* added features *****************************************//
  //adding vertex to left and right face
  //left->addVertex(vertexNew);
  //right->addVertex(vertexNew);
//********************* end added features *****************************************//
  setOrbitOrg(edgeNew->Sym(), vertexNew);

  return edgeNew;
}

void Cell::killVertexEdge(Edge *edge)
{
  assert(edge!=0);

  // locate _edge1_ and _edge2_ as in _makeVertexEdge_

  Edge *edge1 = edge->Oprev();
  Edge *edge2 = edge->Lnext();

  // use _edge1_ for _edge2_ if the destination vertex is isolated

  if (edge2==edge->Sym())
    edge2 = edge1;

  // inverse of _makeVertexEdge_

  Edge::splice(edge1, edge->Sym());
  Edge::splice(edge2, edge);

  // all edges leaving the destination orbit of the deleted edge now have its
  // origin vertex as their vertex of origin

  setOrbitOrg(edge2, edge1->Org());

  // don't use the deleted edge as a reference edge any more

  edge1->Org()->addEdge(edge1);
  edge1->Left()->addEdge(edge1);
  edge2->Left()->addEdge(edge2);

  // reclaim the vertex and the edge

  Vertex::kill(edge->Dest());
  Edge::kill(edge);
}

Edge *Cell::makeFaceEdge(Face *face, Vertex *org, Vertex *dest)
{
  assert(face!=0);
  assert(org!=0);
  assert(dest!=0);

  // locate the edges leaving each of the vertices in the orbit of the face

  Edge *edge  = face->getEdge();
  Edge *edge1 = getOrbitOrg(edge, org);
  Edge *edge2 = getOrbitOrg(edge, dest);

  if (edge1==0)
  {
    (void)fprintf(stderr, "Cell::makeFaceEdge: unable to locate origin vertex %u on face %u\n",
		 org->getID(), face->getID());
    abort();
  }

  if (edge2==0)
  {
    (void)fprintf(stderr, "Cell::makeFaceEdge: unable to locate destination vertex %u on face %u\n",
		 dest->getID(), face->getID());
    abort();
  }

  // create a new face

  Face *faceNew = Face::make(this);

  // create a new (non-loop) edge

  Edge *edgeNew = Edge::make();

  // connect the destination of the new edge to the origin of _edge2_
  // both faces of the edge are now _face_

  Edge::splice(edge2, edgeNew->Sym());

  // connect the origin of the new edge to _edge1_
  // _face_ is split in half along the new edge, with the new face introduced
  // on the right

  Edge::splice(edge1, edgeNew);

  // initialize the secondary attributes of the new edge

  edgeNew->setOrg(edge1->Org());
  edgeNew->setDest(edge2->Org());
  edgeNew->setLeft(edge2->Left());

  // all edges in the right orbit of the new edge (i.e. the left orbit of its
  // Sym) now have the new face as their left face

  setOrbitLeft(edgeNew->Sym(), faceNew);

  return edgeNew;
}

void Cell::killFaceEdge(Edge *edge)
{
  assert(edge!=0);

  // locate _edge1_ and _edge2_ as in _makeFaceEdge_

  Edge *edge1 = edge->Oprev();
  Edge *edge2 = edge->Lnext();

  // use _edge2_ for _edge1_ if the right face is inside a loop

  if (edge1==edge->Sym())
    edge1 = edge2;

  // inverse of _makeFaceEdge_

  Edge::splice(edge2, edge->Sym());
  Edge::splice(edge1, edge);

  // all edges in the right orbit of the deleted edge now have its left face
  // as their left face

  setOrbitLeft(edge1, edge2->Left());

  // don't use the deleted edge as a reference edge any more

  edge1->Org()->addEdge(edge1);
  edge2->Org()->addEdge(edge2);
  edge2->Left()->addEdge(edge2);

  // reclaim the face and the edge

  Face::kill(edge->Right());
  Edge::kill(edge);
}

/* -- public instance methods ---------------------------------------------- */

void Cell::addVertex(Vertex *vertex)
{
  assert(vertex!=0);

  // expand the vertex array, if necessary

  if (vertexCount>=vertexSize)
  {
    unsigned int vertexSizeNew = vertexSize*2;
    Vertex **verticesNew   = new Vertex*[vertexSizeNew];

    for (unsigned int i = 0; i<vertexCount; i++)
      verticesNew[i] = vertices[i];

    delete[] vertices;

    vertices   = verticesNew;
    vertexSize = vertexSizeNew;
  }

  // add the vertex to the end of the array

  vertices[vertexCount++] = vertex;
}


void Cell::removeVertex(Vertex *vertex)
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


void Cell::addFace(Face *face)
{
  assert(face!=0);

  // expand the face array, if necessary

  if (faceCount>=faceSize)
  {
    unsigned int faceSizeNew = faceSize*2;
    Face       **facesNew    = new Face*[faceSizeNew];

    for (unsigned int i = 0; i<faceCount; i++)
      facesNew[i] = faces[i];

    delete[] faces;

    faces    = facesNew;
    faceSize = faceSizeNew;
  }

  // add the face to the end of the array

  faces[faceCount++] = face;
}

void Cell::removeFace(Face *face)
{
  assert(face!=0);

  // locate the face in the array and replace it with the current last face
  // if already the last face, just overwrite it
  // slow: should make this a doubly-linked list ???

  for (unsigned int i = faceCount; i>0; i--)
    if (faces[i-1]==face)
    {
      faces[i-1] = faces[--faceCount];
      return;
    }

  assert(0);
}

/* -- protected instance methods ------------------------------------------- */

Cell::Cell():gaussianWidth(1.0), //initialising the gaussian Widht to 1
 randomNumberGeneratorType(gsl_rng_default) // initialising the random Num generator to be Mersenne Twister
{
  // preallocate enough elements for a cube

  vertices    = new Vertex*[8];
  vertexCount = 0;
  vertexSize  = 8;
  vertexID    = 1;

  faces     = new Face*[6];
  faceCount = 0;
  faceSize  = 6;
  faceID    = 1;
  divisionCounter = 0;
  divisionFactor = 1.1;
  convexAngleThreshold = 180.;
  //setting the random number generator
  // intialised in Initialising list :-> randomNumberGeneratorType = gsl_rng_default;//this is Mersenne Twister algorithm
  randomNumberGenerator = gsl_rng_alloc(randomNumberGeneratorType);
  gsl_rng_set(randomNumberGenerator,38270);//some number as seed-> this can be set with another random number/time or /dev/random /urandom

  // RANDOM NUMBER GENERATOR FOR CELL DIVISION
  cellDivisionRandomNumberGenerator = gsl_rng_alloc(randomNumberGeneratorType);
  gsl_rng_set(cellDivisionRandomNumberGenerator, random_seed());//using /dev/urandom to seed this generator
  //calculating the square root of epsilon
  sqrtEpsilon = sqrt(std::numeric_limits<double>::epsilon());
}

Cell::~Cell()
{
  // reclaim each of the vertices and faces still owned by the cell
  // go in backwards order so that when the elements try to remove themselves,
  // it will be linear time

  {
    for (unsigned int i = vertexCount; i>0; i--)
      Vertex::kill(vertices[i-1]);
  }

  {
    for (unsigned int i = faceCount; i>0; i--)
      Face::kill(faces[i-1]);
  }
  
    // reclaim the vertex and face arrays

  delete[] vertices;
  delete[] faces;
  //release the randomNumberGenerator
  gsl_rng_free(randomNumberGenerator);
  gsl_rng_free(cellDivisionRandomNumberGenerator);

}

/* -- private instance methods --------------------------------------------- */

Edge *Cell::getOrbitOrg(Edge *edge, Vertex *org)
{
  assert(edge!=0);
  assert(org!=0);

  // traverse the Lnext orbit of _edge_ looking for an edge whose origin is
  // _org_

  Edge *scan = edge;

  do
  {
    if (scan->Org()==org)
      return scan;

    scan = scan->Lnext();
  }
  while (scan!=edge);

  return 0;
}

void Cell::setOrbitOrg(Edge *edge, Vertex *org)
{
  assert(edge!=0);
  assert(org!=0);

  // traverse the Onext orbit of _edge_, setting the origin of each edge to
  // _org_

  Edge *scan = edge;

  do
  {
    scan->setOrg(org);

    scan = scan->Onext();
  }
  while (scan!=edge);
}

Edge *Cell::getOrbitLeft(Edge *edge, Face *left)
{
  assert(edge!=0);
  assert(left!=0);

  // traverse the Onext orbit of _edge_ looking for an edge whose left face is
  // _left

  Edge *scan = edge;

  do
  {
    if (scan->Left()==left)
      return scan;

    scan = scan->Onext();
  }
  while (scan!=edge);

  return 0;
}

void Cell::setOrbitLeft(Edge *edge, Face *left)
{
  assert(edge!=0);
  assert(left!=0);

  // traverse the Lnext orbit of _edge_, setting the left face of each edge to
  // _left_

  Edge *scan = edge;

  do
  {
    scan->setLeft(left);

    scan = scan->Lnext();
  }
  while (scan!=edge);
}
// **************************************************************************** //
//structure to store X,Y,Z coordinate of vertex

struct vertex_coordinate {
    unsigned int id;
    double x,y,z;
};
// *********************************************************************************************** //

// Function to relax the Tissue
int Cell::relax()
{
  /* // Under Construction still
   this->functionCallCounter += 1;
   // Step 1. Gathering all the coordinates of the vertices
      int numVertices = (int) this->countVertices();
      int numOfLayer = this->getLayer();
      double radius =  (numOfLayer>1)*(np.sqrt(3.)*(numOfLayer-1)-Length)+Length;
      int perimeterVertexNum = 6+12*(numOfLayer-1);//number of verts in permeter
      //making the array to store all the vertices
      std::vector<double> coordinates(numVertices*3);
      //iterating the vertices of cell to gather coordinates
      CellVertexIterator vertices(this);
      Vertex * vertex;
      int counter = 0;
      while ((vertex = vertices.next())!= 0){
          coordinates[counter]= vertex.getXcoordinate();
          counter += 1;
          coordinates[counter]= vertex.getYcoordinate();
          counter += 1;
          coordinates[counter]= vertex.getZcoordinate();
          counter += 1;
      }
      // *************************************************************** //
      //BOUNDS FOR OPTIMIZATION //
      // *************************************************************** //
        double factor(4.), bound(0.);
        double epsbasefreedom = (1./10.)*2.*M_PI*radius/perimeterVertexNum;
        //array to store the bounds
        std::vector<double> low_bound(numVertices*3);
        std::vector<double> up_bound(numVertices*3);
        //bounds for a Dome //
        CellVertexIterator vertices(this);
        int vertcounter = 0;
        while ((vertex = vertices.next()!= 0))
            {
                if ((vertex->getID() > perimeterVertexNum) and (vertex->getID() <= numVertices))
                { //DOME VERTEX
                  bound = factor*radius;
                  up_bound[getIndex(vertcounter,0)] = bound;
                  up_bound[getIndex(vertcounter,1)] = bound;
                  up_bound[getIndex(vertcounter,2)] = bound;
                  low_bound[getIndex(vertcounter,0)] = -1.* bound;
                  low_bound[getIndex(vertcounter,1)] = -1.* bound;
                  low_bound[getIndex(vertcounter,2)] = 0.;
                }
                else
                {
                  // CYLINDRICAL FLANK VERTEX
                  up_bound[getIndex(vertcounter,0)] = coordinates[getIndex(vertcounter,0)]+epsbasefreedom;
                  up_bound[getIndex(vertcounter,1)] = coordinates[getIndex(vertcounter,1)]+epsbasefreedom;
                  up_bound[getIndex(vertcounter,2)] = factor*radius;
                  low_bound[getIndex(vertcounter,0)] = coordinates[getIndex(vertcounter,0)]-epsbasefreedom;
                  low_bound[getIndex(vertcounter,1)] = coordinates[getIndex(vertcounter,1)]-epsbasefreedom;
                  low_bound[getIndex(vertcounter,2)] = 0.;
                }
                vertcounter++;
            }
        // *************************************************************** //
        //BOUNDS FOR FLOOR VERTICES//
        // *************************************************************** //
        CellVertexIterator vertices(this);
        int vertcounter = 0;
        while ((vertex = vertices.next()!= 0))
            {
                if (checkExternalVertex(vertex))//If the vertex is external vertex then new bounds applied
                {
                  up_bound[getIndex(vertcounter,0)] = coordinates[getIndex(vertcounter,0)]+epsbasefreedom;
                  up_bound[getIndex(vertcounter,1)] = coordinates[getIndex(vertcounter,1)]+epsbasefreedom;
                  up_bound[getIndex(vertcounter,2)] = coordinates[getIndex(vertcounter,2)]+epsbasefreedom;
                  low_bound[getIndex(vertcounter,0)] = coordinates[getIndex(vertcounter,0)]-epsbasefreedom;
                  low_bound[getIndex(vertcounter,1)] = coordinates[getIndex(vertcounter,1)]-epsbasefreedom;
                  low_bound[getIndex(vertcounter,2)] = 0.;
                }
                vertcounter++;
            }
      // *************************************************************** //
      //                    NLOPT CONFIGURATION                         //
      // *************************************************************** //
      //Initialising NLOPT optimizer -> Here using SBPLX
      nlopt:opt opt(nlopt::LN_SBPLX,numVertices*3);
      //setting the objective function
      opt.set_min_objective(objectiveEnergyFunction,NULL);
      //setting the tolerance 
      opt.set_xtol_abs(this->tolerance);
      //setting the initial step
      opt.set_initial_step(this->initialStep);
      //setting the bounds for calculation
      opt.set_lower_bounds(low_bound);
      opt.set_upper_bounds(up_bound);
      // now optimizing the bounds
      nlopt::result result = opt.optimize(coordinates)
      //checking optimisation result 
      */
      return 0;
      // *************************************************************** //
      //                     Releasing the Memory                       //
      // *************************************************************** //
}

// ********************************************************************************** //
//this is the function to use as objective function in calculation of the energy for
//NLOPT relaxation
/*
double objectiveEnergyFunction(const double * inputcoordinates, double * grad, Cell * cell)
{
    CellVertexIterator vertices(cell);
    Vertex * vertex;
    int counter(0);
    //Updating the coordinates of the Cell
    while((vertex = vertices.next())!= 0)
    {
        vertex->setXcoordinate(inputcoordinates[getIndex(counter,0)]);
        vertex->setYcoordinate(inputcoordinates[getIndex(counter,1)]);
        vertex->setZcoordinate(inputcoordinates[getIndex(counter,2)]);
        counter += 1;
    }
    //now setting the parameters with respect to new Coordinates
    cell->setParameters();
    //now calculatng the Energy with Cartesian Volume and return
    return cell->getEnergyCartesianVolume();
}

*/
// ********************************************************************************** //
bool checkExternalVertex(Vertex * vertex)
{
    VertexEdgeIterator edges(vertex);
    Edge * edge;
    Face * face;
    while ((edge = edges.next())!= 0)
    {
          face = edge->Left();
          if (face->getID() == 1)
          {
              return true;
          }
    }
    return false;
}




