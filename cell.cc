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

//for eigenvalue conputation
#include "./eigen/Eigen/Dense"
#include "./eigen/Eigen/Eigenvalues"

//define 2pi
#define M_2PI 6.283185307179586476925286766559005768394338798750211641949
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


// ********************************************************************************** //
bool checkExternalFace(Face * face)
{
    if (face->getID() == 1){return true;};
    FaceEdgeIterator edges(face);
    Edge * edge;
    Face * faceright;
    while ((edge = edges.next())!= 0)
    {
          faceright = edge->Right();
          if (faceright->getID() == 1)
          {
              return true;
          }
    }
    return false;
}
// ********************************************************************************** //


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
      //totalenergy +=  this->getZeta()*this->getFourthTerm();//subtracting the fourth term : z_proj penalty
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
      //double fourthterm = this->getFourthTerm();
      // First : check if the cell is convex or not
      /*
      if (!(this->isConvex())){ //if the cell is not convex then return HIGHVALUE
          //std::cout<<"Cell Is not Convex ! isConvex -> "<<this->isConvex()<<std::endl;
          return std::numeric_limits<double>::max();
      }
      // Second : check if the cell has bent than threshold
      if (fourthterm > (this->getBendingThreshold())){
        //std::cout<<"Bending Threshold Reached"<<std::endl;
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
      //totalenergy -= (this->getGamma()/this->initialVolume)*this->getCartesianVolume();
      //totalenergy +=  (this->getZeta()/this->initialFourthTerm)*fourthterm;//subtracting the fourth term : z_proj penalty
      totalenergy -= (this->getGamma())*this->getCartesianVolume();
      //totalenergy +=  (this->getZeta())*fourthterm;//subtracting the fourth term : z_proj penalty
      //totalenergy += this->getSigma()*this->getSumEdgeLength();
      totalenergy += this->getBendingEnergy();//adding the Helfrich bending energy
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
  /*
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
  */
    // VOLUME USING DIVERGENCE THEOREM
  CellFaceIterator faces(this);
  Face * face;
  Edge * edge;
  Vertex * second;
  Vertex * third;
  double * normal;
  double tempVolume = 0.;
  double productVector[3];
  double dotproduct,area;
  while ((face= faces.next()) != 0){
    //getting the centroid of face
    double centroid[3] = {face->getXCentralised(),face->getYCentralised(),face->getZCentralised()};
    normal = face->getNormal();
    //-1 is there to correct for normal: As we need outward facing normal but the calculated normal is inward
    dotproduct = (-1*normal[0])*centroid[0] + (-1*normal[1])*centroid[1] +(-1*normal[2])*centroid[2];
    area = face->getAreaOfFace();
    //calculating and summing over temp volume
    tempVolume += dotproduct*area;
  }
  return abs(tempVolume)/3.;
 }
 /**
  * caculating the volume
  */
  double Cell::getCartesianVolume(){
    return this->getNonConvexVolume();//volume computed using triangular faces
  //Divergence theorem  ::
  /*
    * Important Change :: 
              The Volume calculation is changed for removing Volume Pressure dependence on Area of Face
               -> Thus, Volume = <A> * SUM_c (V_c/A_c)
               volume is normalised by cell area and calcualted with respect to average area
  */
    /*
  CellFaceIterator faces(this);
  Face * face;
  //Edge * edge;
  //Vertex * second;
  //Vertex * third;
  double tempVolume = 0.;
  double dotproduct;
  double *normalofface;
  //double areaofface;
  while ((face= faces.next()) != 0){
    if (face->getID()==1) continue;
    //getting the centroid of face
    double centroid[3] = {face->getXCentralised(),face->getYCentralised(),face->getZCentralised()};
    //areaofface = face->getAreaOfFace();
    normalofface = face->getNormal();
    dotproduct = normalofface[0]*centroid[0]+
                 normalofface[1]*centroid[1]+
                 normalofface[2]*centroid[2];
    //tempVolume+= dotproduct*areaofface;
    tempVolume += dotproduct;
    }
  //return 1./3.*abs(tempVolume)*(this->averageFaceArea);
  return 1./3.*abs(tempVolume);
  */
  /*
  //Summing over the triangles
  CellFaceIterator faces(this);
  Face * face;
  Edge * edge;
  Vertex * second;
  Vertex * third;
  double localVolume = 0.;
  double productVector[3];
  double triCen[3];
  double vector1[3], vector2[3];
  double sumArea = 0;
  double counter = 0;
  double triArea = 0;
  double norm =0.;
  double trivec1[3], trivec2[3], trivec3[3];
  double trivecarea[3],sumtrivecarea;
  double faceVolume;
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
          double secondCoordinate[3] = {second->getXcoordinate(),
                                        second->getYcoordinate(),
                                        second->getZcoordinate()};
          double thirdCoordinate[3] =  {third->getXcoordinate(),
                                        third->getYcoordinate(),
                                        third->getZcoordinate()};
          //now calculating the Center of triangle : Weighted Centroid
          triCen[0] = 1./3.*(centroid[0]+secondCoordinate[0]+thirdCoordinate[0]);
          triCen[1] = 1./3.*(centroid[1]+secondCoordinate[1]+thirdCoordinate[1]);
          triCen[2] = 1./3.*(centroid[2]+secondCoordinate[2]+thirdCoordinate[2]);
          /*
          // centers of smaller triangles
          double firstcen[3] = {1./3.*(secondCoordinate[0]+thirdCoordinate[0]+triCen[0]),
                                1./3.*(secondCoordinate[1]+thirdCoordinate[1]+triCen[1]),
                                1./3.*(secondCoordinate[2]+thirdCoordinate[2]+triCen[2])};
          double secondcen[3] = {1./3.*(centroid[0]+thirdCoordinate[0]+triCen[0]),
                                 1./3.*(centroid[1]+thirdCoordinate[1]+triCen[1]),
                                 1./3.*(centroid[2]+thirdCoordinate[2]+triCen[2])};
          double thirdcen[3] = {1./3.*(centroid[0]+secondCoordinate[0]+triCen[0]),
                                1./3.*(centroid[1]+secondCoordinate[1]+triCen[1]),
                                1./3.*(centroid[2]+secondCoordinate[2]+triCen[2])};
          
          //trivec1 
          trivec1[0] = secondCoordinate[0]-triCen[0];
          trivec1[1] = secondCoordinate[1]-triCen[1];
          trivec1[2] = secondCoordinate[2]-triCen[2];
          //trivec2
          trivec2[0] = thirdCoordinate[0]-triCen[0];
          trivec2[1] = thirdCoordinate[1]-triCen[1];
          trivec2[2] = thirdCoordinate[2]-triCen[2];
          // trivec3
          trivec3[0] = centroid[0]-triCen[0];
          trivec3[1] = centroid[1]-triCen[1];
          trivec3[2] = centroid[2]-triCen[2];
          // area of small triangles
          //mini triangle 1
          productVector[0] = trivec1[1]*trivec2[2] - trivec1[2]*trivec2[1];
          productVector[1] = trivec1[2]*trivec2[0] - trivec1[0]*trivec2[2];
          productVector[2] = trivec1[0]*trivec2[1] - trivec1[1]*trivec2[0];
          trivecarea[0] = 0.5*sqrt(pow(productVector[0],2)+pow(productVector[1],2)+pow(productVector[2],2));
          //mini triangle 2
          productVector[0] = trivec2[1]*trivec3[2] - trivec2[2]*trivec3[1];
          productVector[1] = trivec2[2]*trivec3[0] - trivec2[0]*trivec3[2];
          productVector[2] = trivec2[0]*trivec3[1] - trivec2[1]*trivec3[0];
          trivecarea[1] = 0.5*sqrt(pow(productVector[0],2)+pow(productVector[1],2)+pow(productVector[2],2));
          //mini triangle 3
          productVector[0] = trivec3[1]*trivec1[2] - trivec3[2]*trivec1[1];
          productVector[1] = trivec3[2]*trivec1[0] - trivec3[0]*trivec1[2];
          productVector[2] = trivec3[0]*trivec1[1] - trivec3[1]*trivec1[0];
          trivecarea[2] = 0.5*sqrt(pow(productVector[0],2)+pow(productVector[1],2)+pow(productVector[2],2));
          sumtrivecarea = trivecarea[0]+trivecarea[1]+trivecarea[2];
          // Weighted triCen
          triCen[0] = 1./sumtrivecarea*(firstcen[0]*trivecarea[0]+
                                        secondcen[0]*trivecarea[1]+
                                        thirdcen[0]*trivecarea[2]);
          triCen[1] = 1./sumtrivecarea*(firstcen[1]*trivecarea[0]+
                                        secondcen[1]*trivecarea[1]+
                                        thirdcen[1]*trivecarea[2]);
          triCen[2] = 1./sumtrivecarea*(firstcen[2]*trivecarea[0]+
                                        secondcen[2]*trivecarea[1]+
                                        thirdcen[2]*trivecarea[2]);
          
          // Get the vectors on the Triangle
          //Vector 1
          vector1[0] = secondCoordinate[0]-centroid[0];
          vector1[1] = secondCoordinate[1]-centroid[1];
          vector1[2] = secondCoordinate[2]-centroid[2];
          //Vector 2
          vector2[0] = thirdCoordinate[0]-centroid[0];
          vector2[1] = thirdCoordinate[1]-centroid[1];
          vector2[2] = thirdCoordinate[2]-centroid[2];
          //cross product of the two vectors
          productVector[0] = vector1[1]*vector2[2] - vector1[2]*vector2[1];
          productVector[1] = vector1[2]*vector2[0] - vector1[0]*vector2[2];
          productVector[2] = vector1[0]*vector2[1] - vector1[1]*vector2[0];
          //maginitude of cross product :: 2*Area Of this Triangle
          norm = sqrt(pow(productVector[0],2)+pow(productVector[1],2)+pow(productVector[2],2));
          triArea= 1./2.*norm;
          // Normal of Triangle
          productVector[0] /= norm;
          productVector[1] /= norm;
          productVector[2] /= norm;
          // Now calculating the height and then Volume of this tetrahedron by projecting the triCen to Normal
          faceVolume =  1./3.*(triCen[0]*productVector[0]+
                               triCen[1]*productVector[1]+
                               triCen[2]*productVector[2]);
           localVolume += faceVolume;
           face->setFaceVolume(faceVolume);
          // Summing the Areas of triangle and counting the number of triangles
          sumArea += triArea;
          counter += 1;
        }
        //std::cout<<"Volume " << tempVolume<<std::endl;
  }
*/
//return localVolume*(sumArea/counter);
//return localVolume;
    /*
  // Summing the distance to origin of all vertices
  CellVertexIterator vertices(this);
  Vertex *vertex;
  //Edge * edge;
  //Vertex * second;
  //Vertex * third;
  double tempVolume = 0.;
  double dotproduct;
  double *normalofface;
  //double areaofface;
  while ((vertex= vertices.next()) != 0){
    //if (face->getID()==1) continue;
    //getting the centroid of face
    double coordinates[3] = {vertex->getXcoordinate(),
                              vertex->getYcoordinate(),
                              vertex->getZcoordinate()};
    tempVolume += sqrt(pow(coordinates[0],2)+pow(coordinates[1],2)+pow(coordinates[2],2));
    //tempVolume+= dotproduct*areaofface;
    }
  //return 1./3.*abs(tempVolume)*(this->averageFaceArea);
  return abs(tempVolume);
  */
 }

 /**
  * Volume calculation using the non-convex generalised method
  * this assumes faces have arranged triangular faces, which is true 
  * if faces are iterated as triangle with respect to centroid of face
  */
 double Cell::getNonConvexVolume(){
  //Summing over the triangles
  CellFaceIterator faces(this);
  Face * face;
  Edge * edge;
  Vertex * second;
  Vertex * third;
  double determinantsum = 0.;
  Eigen::Matrix3f trimat;
  while ((face= faces.next()) != 0){
    // iterating the faces
    FaceEdgeIterator edges(face);
    while ((edge = edges.next()) != 0){
        // getting the 2 vertex to form triangle with the centroid
        // as the edge iteration is always anti-clockwise
        // the repeated edges are always oriented in opposite direction
        // when they are counted twice
        second = edge->Org();
        third = edge->Dest();
        //Forming the matrix to calcualted determinant
        // 3x3 matrix is of form : {A,B,C}, with A,B,C as column vector for the vertex of this triangle
        trimat<< second->getXcoordinate(), third->getXcoordinate(), face->getXCentralised(),
                second->getYcoordinate(), third->getYcoordinate(), face->getYCentralised(),
                second->getZcoordinate(), third->getZcoordinate(), face->getZCentralised();
        determinantsum += trimat.determinant();// calculated the derterminant and summing it up
      }
    }
  // after summing over all the triangular faces on the surface of this polyhedron returning the volume
    return abs(determinantsum/6.);
 }

double Cell::getProjectedCoordinateVolume(){

  //Divergence theorem  ::
  /*
    * Important Change :: 
              The Volume calculation is changed for removing Volume Pressure dependence on Area of Face
               -> Thus, Volume = <A> * SUM_c (V_c/A_c)
               volume is normalised by cell area and calcualted with respect to average area
  */
  CellFaceIterator faces(this);
  Face * face;
  Edge * edge;
  Vertex * second;
  Vertex * third;
  double localVolume = 0.;
  double productVector[3];
  double triCen[3];
  double vector1[3], vector2[3];
  double sumArea = 0;
  double counter = 0;
  double triArea = 0;
  double norm =0.;
  double trivec1[3], trivec2[3], trivec3[3];
  double trivecarea[3],sumtrivecarea;
  double faceVolume;
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
          double secondCoordinate[3] = {second->getXcoordinate(),
                                        second->getYcoordinate(),
                                        second->getZcoordinate()};
          double thirdCoordinate[3] =  {third->getXcoordinate(),
                                        third->getYcoordinate(),
                                        third->getZcoordinate()};
          //now calculating the Center of triangle : Weighted Centroid
          triCen[0] = 1./3.*(centroid[0]+secondCoordinate[0]+thirdCoordinate[0]);
          triCen[1] = 1./3.*(centroid[1]+secondCoordinate[1]+thirdCoordinate[1]);
          triCen[2] = 1./3.*(centroid[2]+secondCoordinate[2]+thirdCoordinate[2]);
          
          vector1[0] = secondCoordinate[0]-centroid[0];
          vector1[1] = secondCoordinate[1]-centroid[1];
          vector1[2] = secondCoordinate[2]-centroid[2];
          //Vector 2
          vector2[0] = thirdCoordinate[0]-centroid[0];
          vector2[1] = thirdCoordinate[1]-centroid[1];
          vector2[2] = thirdCoordinate[2]-centroid[2];
          //cross product of the two vectors
          productVector[0] = vector1[1]*vector2[2] - vector1[2]*vector2[1];
          productVector[1] = vector1[2]*vector2[0] - vector1[0]*vector2[2];
          productVector[2] = vector1[0]*vector2[1] - vector1[1]*vector2[0];
          //maginitude of cross product :: 2*Area Of this Triangle
          norm = sqrt(pow(productVector[0],2)+pow(productVector[1],2)+pow(productVector[2],2));
          triArea= 1./2.*norm;
          // Normal of Triangle
          productVector[0] /= norm;
          productVector[1] /= norm;
          productVector[2] /= norm;
          // Now calculating the height and then Volume of this tetrahedron by projecting the triCen to Normal
          faceVolume =  1./3.*(triCen[0]*productVector[0]+
                               triCen[1]*productVector[1]+
                               triCen[2]*productVector[2]);
           localVolume += faceVolume;
           face->setFaceVolume(faceVolume);
          // Summing the Areas of triangle and counting the number of triangles
          sumArea += triArea;
          counter += 1;
        }
  }
return localVolume*(sumArea/counter);
 }

 //******************************************************************************* //
 void Cell::setAverageFaceArea(){
  CellFaceIterator faces(this);
  Face * face;
  double areaofface = 0.;
  double  facecount = (double) this->countFaces();
  while ((face = faces.next())!= 0){
      if (face->getID()==1) continue;
      areaofface += face->getAreaOfFace();
  }
  this->surfaceArea = areaofface;
  this->averageFaceArea = areaofface/facecount;
 }
//******************************************************************************* //
 void Cell::calculateVertexForce(){
  CellVertexIterator vertices(this);
  Vertex * vertex;
  //std::cout<<"Start of Computation "<<std::endl;
  while ((vertex = vertices.next())!= 0){
    //std::cout<<"vertex ID : "<< vertex->getID()<<std::endl;
    vertex->calculateCartesianForce();
  }
 }
 // ****************************************************************************** //
 void Cell::calculateStressStrain(){
  CellFaceIterator faces(this);
  Face * face;
  while ((face = faces.next())!= 0){
      face->calculateStrain();
      face->calculateStress();
      
  }
 }
 // ****************************************************************************** //
 void Cell::calculateStrain(){
  CellFaceIterator faces(this);
  Face * face;
  while ((face = faces.next())!= 0){
      face->calculateStress();
      face->calculateStrain();
  }
 }
 //******************************************************************************* //
 double Cell::getMeanStrainDeterminant(){
  CellFaceIterator faces(this);
  Face * face;
  double  facecount = (double) this->countFaces();
  double meanstrain(0.);
  while ((face = faces.next())!= 0){
      meanstrain += face->getStrainDeterminant();
  }
  return meanstrain/facecount;
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
 //   Calculating the Helfirch Energy
 //   The mean curvature should be already calculated
void Cell::calculateBendingEnergy(){
    // E = 2(k_b)*sum[Ai_mix(H-H0)^2] 
    // here, omega = k_b
    CellFaceIterator faces(this);
    Face * face;
    Vertex * vertex;
    double bending = 0.;
    double faceomega;
    double vertexomega;
    {
      CellFaceIterator faces(this);
      while((face = faces.next())!= 0){
            if (face->getID() == 1) continue;
            if (face->getOmega() == 0.){
                faceomega = this->getOmega();
             }else{
               faceomega = face->getOmega();
             }
            bending += 2.*faceomega*(face->getAreaMixed())*pow((face->getMeanCurvature() - face->getInitialMeanCurvature()),2);
        }
    }
    // on all the vertices
    {
      CellVertexIterator vertices(this);
      while ((vertex = vertices.next())!= 0){
             if (vertex->getOmega() == 0.){
                vertexomega = this->getOmega();
             }else{
               vertexomega = vertex->getOmega();
             }
            bending += 2.*vertexomega*(vertex->getAreaMixed())*pow((vertex->getMeanCurvature() - vertex->getInitialMeanCurvature()),2);
        }
    }
    //setting bending
    this->bendingEnergy = bending;
}
//********************************************************************************* //
// function to return norm of a vector
double getVectorNorm(double A[]){
    return sqrt(pow(A[0],2)+
                pow(A[1],2)+
                pow(A[2],2)
                );
}
//********************************************************************************* //
// function to return norm of a vector
double getVectorDotProduct(double A[],double B[]){
    return A[0]*B[0]+
           A[1]*B[1]+
           A[2]*B[2];
}
//********************************************************************************* //
// function to return cot(angle)
double getCot(double angle){
  return tan(M_PI_2 - angle);//as cot(angle) = tan(pi/2 - angle)
}
//********************************************************************************* //
//Function to check the triangle is obtuse or not
// return : True if obtuse
//          False if not obtuse
bool checkObtuseTriangle(double A[],double B[], double C[]){ 
    // getting the vector
    double AC[3] = {C[0]- A[0],C[1]- A[1],C[2]- A[2]};
    double AB[3] = {B[0]- A[0],B[1]- A[1],B[2]- A[2]};
    double BC[3] = {C[0]- B[0],C[1]- B[1],C[2]- B[2]};
    double BA[3] = {-1.*AB[0],-1.*AB[1],-1.*AB[2]};
    double CA[3] =  {-1.*AC[0],-1.*AC[1],-1.*AC[2]};
    double CB[3] =  {-1.*BC[0],-1.*BC[1],-1.*BC[2]};

    // dot products
    double d1 = getVectorDotProduct(AC,AB);
    double d2 = getVectorDotProduct(BC,BA);
    double d3 = getVectorDotProduct(CA,CB);
    if (d1*d2*d3 < 0.) {//meaning one of the dot product is negative
                      // or the angle between any two of the vectors should be more than 90 degrees
      return true;//this triangle is obtuse
    }
    return false;// this triangle is not obtuse
}
//********************************************************************************* //
//Function to get angle between centroid (reference point) and other two vertex
double getAngle(double A[],double B[], double C[]){ 
    // getting the vector
    double AC[3] = {C[0]- A[0],C[1]- A[1],C[2]- A[2]};
    double AB[3] = {B[0]- A[0],B[1]- A[1],B[2]- A[2]};
    //Angle at A
    return acos(getVectorDotProduct(AC,AB)/(getVectorNorm(AC)*getVectorNorm(AB))); 
}
//********************************************************************************* //
//Function to check if angle between centroid (reference point) and other two vertex is obtuse or not 
bool checkAngleObtuse(double A[],double B[], double C[]){ 
    // getting the vector
    double AC[3] = {C[0]- A[0],C[1]- A[1],C[2]- A[2]};
    double AB[3] = {B[0]- A[0],B[1]- A[1],B[2]- A[2]};
    //Angle at A
    double angle = acos(getVectorDotProduct(AC,AB)/(getVectorNorm(AC)*getVectorNorm(AB)));
    if (angle > (M_PI_2)) {//meaning angle is greater than 90 degrees or is Obtuse
      return true;//this angle at reference point is obtuse
    }
    return false;// this angle at reference point is not obtuse
}
//********************************************************************************* //
 // Function to calculate the area of a given triangle
double getAreaOfTriangle(double A[],double B[], double C[]){ 
    // getting the vector
    double AC[3] = {C[0]- A[0],C[1]- A[1],C[2]- A[2]};
    double AB[3] = {B[0]- A[0],B[1]- A[1],B[2]- A[2]};
    double BC[3] = {C[0]- B[0],C[1]- B[1],C[2]- B[2]};
    // Area = 1/2*(Norm(CrossProduct of Two vectors of triangle))
    //calculating cross product of two vectors ABxAC
    double productVector[3];
      productVector[0] = AB[1]*AC[2] - AB[2]*AC[1];
      productVector[1] = AB[2]*AC[0] - AB[0]*AC[2];
      productVector[2] = AB[0]*AC[1] - AB[1]*AC[0];
    //Calculating area and returning
    return 0.5*(sqrt(pow(productVector[0],2)+
                      pow(productVector[1],2)+
                      pow(productVector[2],2)
                    )
                );
}
//********************************************************************************* //
 // Function to calculate non-obtuse voronoi area for a given vertex
 //   IMPORTANT :: VERTEX A (passed below) is the reference for calculating the voronoi component
//     (reference) A ------- B
//                  \       /
//                   \     /
//                    \   /
//                     \ /
//                      C
//********************************************************************************* //
double getAreaVoronoiOfTriangle(double A[],double B[], double C[]){ 
    // getting the vector
    double CA[3] = {A[0]- C[0],A[1]- C[1],A[2]- C[2]};
    double CB[3] = {B[0]- C[0],B[1]- C[1],B[2]- C[2]};
    double BA[3] = {A[0]- B[0],A[1]- B[1],A[2]- B[2]};
    double BC[3] = {-1*CB[0],-1*CB[1],-1*CB[2]};
    // Getting anlge on the either side of triangle (angle B, angle C)
    double angleC = acos(getVectorDotProduct(CB,CA)/(getVectorNorm(CB)*getVectorNorm(CA)));
    double angleB = acos(getVectorDotProduct(BA,BC)/(getVectorNorm(BA)*getVectorNorm(BC)));
    // now calculating the voronoi area
    return 1./8. * (pow(getVectorNorm(BA),2)*getCot(angleC)+
                    pow(getVectorNorm(CA),2)*getCot(angleB)
                   );
}
//********************************************************************************* //
//Function to calculate area Mixed with three given vertices, centroid needs to be the reference vertex !!
// check if the triangle is obtuse
  double calculateAreaMixedComponent(double centroid[], double vertex1coordinate[],double vertex2coordinate[]) {
          double areacomponent  = 0.;
          if (checkObtuseTriangle(centroid, vertex1coordinate,vertex2coordinate)){
            //std::cout<< " Check Obtuse True" ;
                //True == obtuse triangle, so special calculation
                // Check if the angle at centroid (reference point) is obtuse or not
                if (checkAngleObtuse(centroid,vertex1coordinate,vertex2coordinate)){
                      //True == angle at centroid for this triangle is Obtuse, hence areaMixed += area(T)/2
                     areacomponent= 0.5*getAreaOfTriangle(centroid, vertex1coordinate,vertex2coordinate);
                     //std::cout<< ":: angle at centroid Obtuse" << areacomponent;
                  }
                else{
                      //False == angle at centroid is not obtuse, and hence areaMixed += area(T)/4
                      areacomponent= 0.25*getAreaOfTriangle(centroid, vertex1coordinate,vertex2coordinate);
                      //std::cout<< ":: angle at centroid not Obtuse" << areacomponent;
                }
            }
          else{
             //std::cout<< " Check Obtuse False" ;
                //False == this triangle is not obtuse triangle, so voronoi component for this vertex is calculated
                areacomponent=  getAreaVoronoiOfTriangle(centroid, vertex1coordinate, vertex2coordinate);
            }
          return areacomponent;
  }


 

//********************************************************************************* //
 // Function to calculate Laplace-Beltrami Operator for a given vertex
 //   IMPORTANT :: VERTEX A (passed below) is the reference for calculating the voronoi component
//     (reference) A ------- B
//                  \       /
//                   \     /
//                    \   /
//                     \ /
//                      C
//  LBOperator[3] would be the array that will be edited with the calculated LB values
//********************************************************************************* //
void getLaplaceBeltrami(double A[],double B[], double C[], double LBOperator[]){ 
    // getting the vector
    double CA[3] = {A[0]- C[0],A[1]- C[1],A[2]- C[2]};
    double CB[3] = {B[0]- C[0],B[1]- C[1],B[2]- C[2]};
    double BA[3] = {A[0]- B[0],A[1]- B[1],A[2]- B[2]};
    double BC[3] = {-1*CB[0],-1*CB[1],-1*CB[2]};
    // Getting anlge on the either side of triangle (angle B, angle C)
    double angleC = acos(getVectorDotProduct(CB,CA)/(getVectorNorm(CB)*getVectorNorm(CA)));
    double angleB = acos(getVectorDotProduct(BA,BC)/(getVectorNorm(BA)*getVectorNorm(BC)));
    // Calculating half of cot(angles)
    double cotC = getCot(angleC);
    double cotB = getCot(angleB);
    // calculating the Laplace-Beltrami Operator
    LBOperator[0] = cotB*CA[0] + cotC*BA[0];
    LBOperator[1] = cotB*CA[1] + cotC*BA[1];
    LBOperator[2] = cotB*CA[2] + cotC*BA[2];  
};
//********************************************************************************* //
 // Function to calculate Laplace-Beltrami Operator on Mean Curvature for a given vertex
 //   IMPORTANT :: VERTEX A (passed below) is the reference for calculating the voronoi component
//     (reference) A ------- B
//                  \       /
//                   \     /
//                    \   /
//                     \ /
//                      C
//  this would return a scalar
//********************************************************************************* //
double getLaplaceBeltramiOnMeanCurvature(double A[],double B[], double C[], double curvatureDiffBA,double curvatureDiffCA){ 
    // getting the vector
    double CA[3] = {A[0]- C[0],A[1]- C[1],A[2]- C[2]};
    double CB[3] = {B[0]- C[0],B[1]- C[1],B[2]- C[2]};
    double BA[3] = {A[0]- B[0],A[1]- B[1],A[2]- B[2]};
    double BC[3] = {-1*CB[0],-1*CB[1],-1*CB[2]};
    // Getting anlge on the either side of triangle (angle B, angle C)
    double angleC = acos(getVectorDotProduct(CB,CA)/(getVectorNorm(CB)*getVectorNorm(CA)));
    double angleB = acos(getVectorDotProduct(BA,BC)/(getVectorNorm(BA)*getVectorNorm(BC)));
    // Calculating half of cot(angles)
    double cotC = getCot(angleC);
    double cotB = getCot(angleB);
    // calculating the Laplace-Beltrami Operator
    return cotB*curvatureDiffCA + cotC*curvatureDiffBA;
};
//********************************************************************************* //
void Cell::setMeanCurvature(){
  Face * face;
  Edge * edge1,*edge2;
  Vertex * vertex1,*vertex2, *vertex3; 
  double totalAreaMixed = 0.;
  // First iterating through all the faces to set areaMixed for the centroid of faces
  CellFaceIterator faces(this);
  while((face = faces.next()) != 0){
      if (face->getID() == 1) continue; // not calculating for the external face
                                        // imagine this is constant
      FaceEdgeIterator edges(face);
      double areaMixed = 0.;
      double sumAngle = 0.;
      double centroid[3] = {face->getXCentralised(),
                            face->getYCentralised(),
                            face->getZCentralised()};
      double LBOperator[3] = {0.,0.,0.};
      while ((edge1 = edges.next())!= 0){
              //getting vertices
              /*          centroid
                          /     \ 
                         /       \
                      v1 --------v2

              */
              vertex1 = edge1->Org();
              vertex2 = edge1->Dest();
              double vertex1coordinate[3] =  {vertex1->getXcoordinate(),
                                              vertex1->getYcoordinate(),
                                              vertex1->getZcoordinate()};
              double vertex2coordinate[3] =  {vertex2->getXcoordinate(),
                                              vertex2->getYcoordinate(),
                                              vertex2->getZcoordinate()};
              //caculate AreaMixedComponent
              areaMixed += calculateAreaMixedComponent(centroid,vertex1coordinate,vertex2coordinate);
              sumAngle += getAngle(centroid,vertex1coordinate,vertex2coordinate);
              //std::cout<<"Face Sum Angle :"<<sumAngle<<std::endl;
              /*if (face->getID() == 113){
                std::cout<<std::endl<<"113  areaMixed component  :v1 =  "<< vertex1->getID()<<" v2 = "<< vertex2->getID()<<" AM = "<< calculateAreaMixedComponent(centroid,vertex1coordinate,vertex2coordinate) << std::endl;
              }*/
              // calculating the Laplace Beltrami operator for this vector on this triangle
              double tempLBOperator[3];
              getLaplaceBeltrami(centroid, vertex1coordinate, vertex2coordinate,tempLBOperator);

              // value is stored in LBO operator
              LBOperator[0] += tempLBOperator[0];
              LBOperator[1] += tempLBOperator[1];
              LBOperator[2] += tempLBOperator[2];
            }
      totalAreaMixed += areaMixed;
      // dividing the LB-Operator by 2*(AreaMixed of this vertex) as final step of its calculation
      LBOperator[0] /= (2.*areaMixed);
      LBOperator[1] /= (2.*areaMixed);
      LBOperator[2] /= (2.*areaMixed);
      //calculating guassian curvature
      //std::cout<<"Gaussian curvature :"<<(M_2PI - sumAngle)/(areaMixed)<<std::endl;
      face->setGaussianCurvature((M_2PI - sumAngle)/(areaMixed));
      //now setting the LBOperator value for this vertex (which is face for this case)
      face->setLBOperator(LBOperator);
      //Calcuating Mean curvature = 0.5*norm(LBO)
      face->setMeanCurvature((0.5*getVectorNorm(LBOperator)));
      //saving the areaMixed for the centroid of this face
      face->setAreaMixed(areaMixed);
      }
  // with this all the faces are iterated over and the areaMixed for all the centroid of faces are summed on areaMixed
  //  =============================================================================================================
  // Now iterating over all the vertices on the Structure
  CellVertexIterator vertices(this);
  Vertex * referenceVertex;
  Face * right, *left;
  while((referenceVertex = vertices.next()) != 0){
      // iterating around all the edges from this vertex
      VertexEdgeIterator edges(referenceVertex);
      double areaMixed = 0.;
      double sumAngle = 0.;
      //referenceVertex is the centroid of its vertex orbit, similar to centroid for face above
      double centroid[3] = {referenceVertex->getXcoordinate(),
                            referenceVertex->getYcoordinate(),
                            referenceVertex->getZcoordinate()};
      double LBOperator[3] = {0.,0.,0.};
      //std::cout<<"Vertex ID :" << referenceVertex->getID()<<" Sum Angle :"<<sumAngle<<std::endl;
      while ((edge1 = edges.next())!= 0){
          //getting vertices
          /*           refVert
                      /  |   \ 
                left /   |e   \ right
                    /    |     \
                  v1 ----v2-----v3
                  where: v1,v3 = face Centroid left and right, respectively
                          e : iterated edge
                          v2 : e->Dest()

          */
          right = edge1->Right();
          left = edge1->Left();
          vertex2 = edge1->Dest();
          double vertex1coordinate[3] =   {left->getXCentralised(),
                                           left->getYCentralised(),
                                           left->getZCentralised()};
          double vertex2coordinate[3] =  {vertex2->getXcoordinate(),
                                          vertex2->getYcoordinate(),
                                          vertex2->getZcoordinate()};
          double vertex3coordinate[3] =  {right->getXCentralised(),
                                          right->getYCentralised(),
                                          right->getZCentralised()};

          // first areaMixedComponent for refVert-v1-v2
          if (left->getID() != 1){// Keeping a check to not calculate the outer face
                areaMixed += calculateAreaMixedComponent(centroid,vertex1coordinate,vertex2coordinate);
                sumAngle += getAngle(centroid,vertex1coordinate,vertex2coordinate);
                //std::cout<<"Vertex Sum Angle :"<<sumAngle<<std::endl;
                // calculating the Laplace Beltrami operator for this vector on this triangle
                double tempLBOperator[3];
                getLaplaceBeltrami(centroid, vertex1coordinate, vertex2coordinate,tempLBOperator);
                // value is stored in LBO operator
                LBOperator[0] += tempLBOperator[0];
                LBOperator[1] += tempLBOperator[1];
                LBOperator[2] += tempLBOperator[2];
                /*std::cout<< std::endl<<"Ref "<<referenceVertex->getID()<<"  Face ::"<<left->getID()<<" Vertex ::"<<vertex2->getID()
           <<"LBO OPERATOR "<< tempLBOperator[0] <<tempLBOperator[1] <<tempLBOperator[2] <<std::endl;
           */
              }
          
          // second areaMixedComponent for refVert-v2-v3
          if (right->getID() != 1){
                areaMixed += calculateAreaMixedComponent(centroid,vertex2coordinate,vertex3coordinate);
                sumAngle += getAngle(centroid,vertex2coordinate,vertex3coordinate);
                //std::cout<<"Vertex Sum Angle :"<<sumAngle<<std::endl;
                // calculating the Laplace Beltrami operator for this vector on this triangle
                double tempLBOperator[3];
                getLaplaceBeltrami(centroid, vertex2coordinate, vertex3coordinate,tempLBOperator);
                // value is stored in LBO operator
                LBOperator[0] += tempLBOperator[0];
                LBOperator[1] += tempLBOperator[1];
                LBOperator[2] += tempLBOperator[2];
                /*
                std::cout<< std::endl<<"Ref "<<referenceVertex->getID()<<"  Face ::"<<right->getID()<<" Vertex ::"<<vertex2->getID()
           <<"LBO OPERATOR "<< tempLBOperator[0] <<tempLBOperator[1] <<tempLBOperator[2] <<std::endl;
           */
              }
          
            }  
          totalAreaMixed += areaMixed;
          // dividing the LB-Operator by 2*(AreaMixed of this vertex) as final step of its calculation
          LBOperator[0] /= (2.*areaMixed);
          LBOperator[1] /= (2.*areaMixed);
          LBOperator[2] /= (2.*areaMixed);
          //calculating guassian curvature
          double tempNorm[3];
          double LBOnorm = getVectorNorm(LBOperator);
          tempNorm[0] = 1./LBOnorm*LBOperator[0];
          tempNorm[1] = 1./LBOnorm*LBOperator[1];
          tempNorm[2] = 1./LBOnorm*LBOperator[2];
          //std::cout<<"area mixed :"<<(areaMixed)<<std::endl;
          //std::cout<<"Gaussian curvature :"<<(M_2PI - sumAngle)/(areaMixed)<<std::endl;
          referenceVertex->setGaussianCurvature((M_2PI - sumAngle)/(areaMixed));
          //now setting the LBOperator value for this vertex (which is face for this case)
          referenceVertex->setLBOperator(LBOperator);
          //Calcuating Mean curvature = 0.5*norm(LBO)
          referenceVertex->setMeanCurvature((0.5*LBOnorm));
          //saving the areaMixed for the centroid of this face
          referenceVertex->setAreaMixed(areaMixed);
          // saving the normal to the vertex
          double * pntnorm = tempNorm;
          referenceVertex->setNormal(pntnorm);
      }
  this->totalAreaMixed = totalAreaMixed;
  // Now Calculate Bending Energy :: on newly calculated mean curvature
  this->calculateBendingEnergy();
}
//********************************************************************************* //
void Cell::setRandomInitialMeanCurvature(){
  CellFaceIterator faces(this);
  Face * face;
  Vertex * vertex;
  while((face= faces.next()) != 0 ){
      face->setRandomInitialMeanCurvature();
  }
  CellVertexIterator vertices(this);
  while((vertex = vertices.next())!= 0){
      vertex->setRandomInitialMeanCurvature();
  }
}
//********************************************************************************* //
void Cell::calculateBendingForce(){
  Face * face;
  Edge * edge1,*edge2;
  Vertex * vertex1,*vertex2, *vertex3; 
  // First iterating through all the faces to set areaMixed for the centroid of faces
  CellFaceIterator faces(this);
  while((face = faces.next()) != 0){
      if (face->getID() == 1) continue; // not calculating for the external face
                                        // imagine this is constant
      FaceEdgeIterator edges(face);
      double centroid[3] = {face->getXCentralised(),
                            face->getYCentralised(),
                            face->getZCentralised()};
      double LBOperator(0.);
      while ((edge1 = edges.next())!= 0){
              //getting vertices
              /*          centroid
                          /     \ 
                         /       \
                      v1 --------v2

              */
              vertex1 = edge1->Org();
              vertex2 = edge1->Dest();
              double vertex1coordinate[3] =  {vertex1->getXcoordinate(),
                                              vertex1->getYcoordinate(),
                                              vertex1->getZcoordinate()};
              double vertex2coordinate[3] =  {vertex2->getXcoordinate(),
                                              vertex2->getYcoordinate(),
                                              vertex2->getZcoordinate()};
              //std::cout<<"Face Sum Angle :"<<sumAngle<<std::endl;
              /*if (face->getID() == 113){
                std::cout<<std::endl<<"113  areaMixed component  :v1 =  "<< vertex1->getID()<<" v2 = "<< vertex2->getID()<<" AM = "<< calculateAreaMixedComponent(centroid,vertex1coordinate,vertex2coordinate) << std::endl;
              }*/
              // calculating the Laplace Beltrami operator for this vector on this triangle
              double tempLBOperator;
              // curvatureDiff = H_i - H_j 
              double curvatureDiffBA = (face->getMeanCurvature() - vertex1->getMeanCurvature());
              double curvatureDiffCA = (face->getMeanCurvature() - vertex2->getMeanCurvature());
              tempLBOperator = getLaplaceBeltramiOnMeanCurvature(centroid, vertex1coordinate, vertex2coordinate,curvatureDiffBA, curvatureDiffCA);
              // value is stored in LBO operator
              LBOperator += tempLBOperator;
            }
      // dividing the LB-Operator by 2*(AreaMixed of this vertex) as final step of its calculation
      double areaMixed = face->getAreaMixed();
      // Computing the LBOperator for this face
      LBOperator /= (2.*areaMixed);
      // saving hte LBOpertor on Mean Curvature
      face->setLBOperatorOnMeanCurvature(LBOperator);
      face->calculateBendingForce();
      }
  // with this all the faces are iterated over and the areaMixed for all the centroid of faces are summed on areaMixed
  //  =============================================================================================================
  // Now iterating over all the vertices on the Structure
  CellVertexIterator vertices(this);
  Vertex * referenceVertex;
  Face * right, *left;
  while((referenceVertex = vertices.next()) != 0){
      // iterating around all the edges from this vertex
      VertexEdgeIterator edges(referenceVertex);
      //referenceVertex is the centroid of its vertex orbit, similar to centroid for face above
      double centroid[3] = {referenceVertex->getXcoordinate(),
                            referenceVertex->getYcoordinate(),
                            referenceVertex->getZcoordinate()};
      double LBOperator(0.);
      //std::cout<<"Vertex ID :" << referenceVertex->getID()<<" Sum Angle :"<<sumAngle<<std::endl;
      while ((edge1 = edges.next())!= 0){
          //getting vertices
          /*           refVert
                      /  |   \ 
                left /   |e   \ right
                    /    |     \
                  v1 ----v2-----v3
                  where: v1,v3 = face Centroid left and right, respectively
                          e : iterated edge
                          v2 : e->Dest()

          */
          right = edge1->Right();
          left = edge1->Left();
          vertex2 = edge1->Dest();
          double vertex1coordinate[3] =   {left->getXCentralised(),
                                           left->getYCentralised(),
                                           left->getZCentralised()};
          double vertex2coordinate[3] =  {vertex2->getXcoordinate(),
                                          vertex2->getYcoordinate(),
                                          vertex2->getZcoordinate()};
          double vertex3coordinate[3] =  {right->getXCentralised(),
                                          right->getYCentralised(),
                                          right->getZCentralised()};

          // first areaMixedComponent for refVert-v1-v2
          if (left->getID() != 1){// Keeping a check to not calculate the outer face
                // curvatureDiff = H_i - H_j 
                double curvatureDiffBA = (referenceVertex->getMeanCurvature() - left->getMeanCurvature());
                double curvatureDiffCA = (referenceVertex->getMeanCurvature() - vertex2->getMeanCurvature());
                double tempLBOperator = getLaplaceBeltramiOnMeanCurvature(centroid, vertex1coordinate, vertex2coordinate,curvatureDiffBA, curvatureDiffCA);
                // value is stored in LBO operator
                LBOperator += tempLBOperator;
              }
          
          // second areaMixedComponent for refVert-v2-v3
          if (right->getID() != 1){
                // curvatureDiff = H_i - H_j 
                double curvatureDiffBA = (referenceVertex->getMeanCurvature() - vertex2->getMeanCurvature());
                double curvatureDiffCA = (referenceVertex->getMeanCurvature() - right->getMeanCurvature());
                double tempLBOperator = getLaplaceBeltramiOnMeanCurvature(centroid, vertex2coordinate, vertex3coordinate,curvatureDiffBA, curvatureDiffCA);
                // value is stored in LBO operator
                LBOperator += tempLBOperator;
              }
          
            }  
          double areaMixed = referenceVertex->getAreaMixed();
          // dividing the LB-Operator by 2*(AreaMixed of this vertex) as final step of its calculation
          LBOperator /= (2.*areaMixed);
          // saving LBOperator on Mixed Area
          referenceVertex->setLBOperatorOnMeanCurvature(LBOperator);
          referenceVertex->calculateBendingForce();
      }
}
 //********************************************************************************* //
 double Cell::getSumEdgeLength(){
  CellFaceIterator faces(this);
    Face * face;
    double sumlength = 0.;
    while((face = faces.next()) != 0){
       sumlength += face->getSumEdgeLength();
    }
  sumlength /= 2.;//As each edge length is calculated for twice
  return sumlength;
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
 //******************************************************************************* //
 void Cell::setGrowthVar(double tempvar){
  growthvar = tempvar;
  this->setFaceGrowthVar();
}
 //********************************************************************************* //
 void Cell::setKappa(double newkappa){
  this->kappa = newkappa;
  this->setFaceKappa();
}
//********************************************************************************* //
 void Cell::setFaceKappa(){
  Face * face;
  {
    CellFaceIterator faces(this);
    while((face = faces.next())!= 0){
          face->setKappa(this->kappa);
    }
  }
 }
 //********************************************************************************* //
 void Cell::calculateAverageTFM(){
  Face * face;
  double facecount = ((double) this->countFaces())-1;//number of faces in this structure (minus face 1)
  double tfm[2][2] = {{0,0},{0,0}};
  {
    CellFaceIterator faces(this);
    while((face = faces.next())!= 0){
      if (face->getID() == 1 ) continue;
          tfm[0][0] += face->targetFormMatrix[0][0];
          tfm[0][1] += face->targetFormMatrix[0][1];
          tfm[1][0] += face->targetFormMatrix[1][0];
          tfm[1][1] += face->targetFormMatrix[1][1];
    }
  }
  tfm[0][0] /= facecount;
  tfm[0][1] /= facecount;
  tfm[1][0] /= facecount;
  tfm[1][1] /= facecount;
  //setting average tfm
  this->averageTargetFormMatrix[0][0] = tfm[0][0];
  this->averageTargetFormMatrix[0][1] = tfm[0][1];
  this->averageTargetFormMatrix[1][0] = tfm[1][0];
  this->averageTargetFormMatrix[1][1] = tfm[1][1];
 }
  //********************************************************************************* //
 void Cell::setAverageTFM(){
  // Calculating the Average TFM for all the faces in this Structure
  this->calculateAverageTFM();
  // Now Setting Average TFM to all the faces
  Face * face;
    {
    CellFaceIterator faces(this);
    while((face = faces.next())!= 0){
      if (face->getID() == 1 ) continue;
      // setting Average Target Form Matrix as  TFM for the faces
      face->targetFormMatrix[0][0] = this->averageTargetFormMatrix[0][0];
      face->targetFormMatrix[0][1] = this->averageTargetFormMatrix[0][1];
      face->targetFormMatrix[1][0] = this->averageTargetFormMatrix[1][0];
      face->targetFormMatrix[1][1] = this->averageTargetFormMatrix[1][1];
      // setting respectively proportional constantMatrix for growth
      face->constantGrowthMatrix[0][0] = 0.01*face->targetFormMatrix[0][0];
      face->constantGrowthMatrix[0][1] = 0.01*face->targetFormMatrix[0][1];
      face->constantGrowthMatrix[1][0] = 0.01*face->targetFormMatrix[1][0];
      face->constantGrowthMatrix[1][1] = 0.01*face->targetFormMatrix[1][1];
      }
    }
  // Now recalculating the energy terms
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
          face->setEnergyTerms();
      }
    }
  /////////////////////////////////
 }
  //********************************************************************************* //
 void Cell::setInitialMeanCurvature(){
    //first calculating the mean curvature
    this->setMeanCurvature();
    // now setting it as initial mean curvature on the structure
    Face * face;
    //Edge * edge;
    Vertex * vertex;
    // on all the faces
  {
    CellFaceIterator faces(this);
    while((face = faces.next())!= 0){
          if (face->getID() == 1) continue;
          face->setInitialMeanCurvature(face->getMeanCurvature());
    }
  }
    // on all the vertices
  {
    CellVertexIterator vertices(this);
    while ((vertex = vertices.next())!= 0){
      vertex->setInitialMeanCurvature(vertex->getMeanCurvature());
    }
  }
 }
 //********************************************************************************* //
 void Cell::setInitialMeanCurvature(double initialCurvature){
    // now setting it as initial mean curvature on the structure
    Face * face;
    //Edge * edge;
    Vertex * vertex;
    // on all the faces
  {
    CellFaceIterator faces(this);
    while((face = faces.next())!= 0){
          if (face->getID() == 1) continue;
          face->setInitialMeanCurvature(initialCurvature);
    }
  }
    // on all the vertices
  {
    CellVertexIterator vertices(this);
    while ((vertex = vertices.next())!= 0){
      vertex->setInitialMeanCurvature(initialCurvature);
    }
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
          //if (face->getID() == 1) continue;
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
  this->initialFourthTerm = this->getFourthTerm();// setting the initial fourth term
  /////////////////////////////////////////
  {
    CellFaceIterator faces(this);
    while((face = faces.next())!= 0){
          //if (face->getID() == 1) continue;
          face->setSumEdgeLength();
    }
  }
  /////////////////////////////////////////
  this->setMeanCurvature();// First calculating the mean curvature
  //this->setInitialMeanCurvature(0.);
  /////////////////////////////////////////
  {
    CellFaceIterator faces(this);
    while((face = faces.next())!= 0){
          //if (face->getID() == 1) continue;
          face->setMu();
    }
  }
  /////////////////////////////////////////
  {
    CellFaceIterator faces(this);
    while((face = faces.next())!= 0){
          //if (face->getID() == 1) continue;
          face->setInitialTargetFormMatrixCurrent();
    }
  }
  /////////////////////////////////////////
  this->calculateAverageTFM(); //calculatign averag TFM for the face
  /////////////////////////////////////////
  {
    CellFaceIterator faces(this);
    while((face = faces.next())!= 0){
          //if (face->getID() == 1) continue;
          face->setConstantGrowthMatrix();
    }
  }
  /////////////////////////////////////////
  
  {
    CellFaceIterator faces(this);
    while((face = faces.next())!= 0){
          //if (face->getID() == 1) continue;
          face->setMu();
    }
  }
  /////////////////////////////////////////
  {
    CellFaceIterator faces(this);
    while((face = faces.next())!= 0){
          //if (face->getID() == 1) continue;
          face->setAreaOfFace();
    }
  }
  /////////////////////////////////////////
  {
    CellFaceIterator faces(this);
    while((face = faces.next())!= 0){
          //if (face->getID() == 1) continue;
          face->setEnergyTerms();
    }
  }
  /////////////////////////////////////////
  {
    CellFaceIterator faces(this);
    while((face = faces.next())!= 0){
          //if (face->getID() == 1) continue;
          face->setDivisionThreshold();
    }
  }
  /////////////////////////////////////////
  // setting external position true for all the faces that are external faces
  {
    CellFaceIterator faces(this);
    while((face = faces.next())!= 0){
          //if (face->getID() == 1) continue;
          face->setExternalPosition(checkExternalFace(face));
    }
  }
  /////////////////////////////////////////
  //seting average face area
  this->setAverageFaceArea();
  // calculating Average TFM for faces
  this->calculateAverageTFM();
  // Setting the bendingThreshold to initial bending energy value
  double fourthterm = this->getFourthTerm();
  this->setBendingThreshold(fourthterm + 0.05*fourthterm);
  // Calculate and set the initial volume of this cell
  this->initialVolume = this->getCartesianVolume();
 }


 //********************************************************************************* //
void Cell::setFaceGrowthVar(){
  Face * face;
  //Edge * edge;
  Vertex * vertex;
  {
    CellFaceIterator faces(this);    
    while((face = faces.next())!= 0){
          //if (face->getID() == 1) continue;
          face->setGrowthVar(this->getGrowthVar());
    }
  }
}
//********************************************************************************* //
void Cell::setRadialOrthoradialVector(Face * primordialFace){ 
  Face * face;
  {
    CellFaceIterator faces(this);
    while((face = faces.next())!= 0){
          face->setRadialOrthoradialVector(primordialFace);
    }
  }
}

void Cell::setRadialOrthoradialStress(){
Face * face;
  {
    CellFaceIterator faces(this);
    while((face = faces.next())!= 0){
          face->setRadialOrthoradialStress();
    }
  }
}

void Cell::setRadialOrthoradialFeedbackCorrection(){
Face * face;
  {
    CellFaceIterator faces(this);
    while((face = faces.next())!= 0){
          face->setRadialOrthoradialFeedbackCorrection();
    }
  }
}
//********************************************************************************* //
 
void Cell::setParameters(){
  Face * face;
  //Edge * edge;
  Vertex * vertex;
  {
    CellFaceIterator faces(this);
    
    while((face = faces.next())!= 0){
          //if (face->getID() == 1) continue;
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
  /*
  ////////////////////////////////////////
  {
    CellFaceIterator faces(this);
    while((face = faces.next())!= 0){
          //if (face->getID() == 1) continue;
          face->setSumEdgeLength();
    }
  }
  ////////////////////////////////////////
  */
  this->setMeanCurvature();
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
         //if (face->getID() == 1) continue;
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
  //seting average face area
  this->setAverageFaceArea();
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
          face->setDivisionThreshold();
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
void Cell::setDivisionThreshold(double newthreshold){
  CellFaceIterator faces(this);
  Face * face;
  while ((face = faces.next())!=0){
        if (face->getID() == 1) continue;
        face->setDivisionThreshold(newthreshold);
      }
}

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

Cell::Cell():gaussianWidth(0.125), //initialising the Standard Deviaton of Gaussian Noise = 0.125
 randomNumberGeneratorType(gsl_rng_default), // initialising the random Num generator to be Mersenne Twister
 seedNumberGeneratorType(gsl_rng_default)
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
  divisionFactor = 1.5;
  convexAngleThreshold = 360.;//2 pi
  this->sigma = 0.;
  this->omega = 0.;
  this->meanCurvatureWidth = 0.1;
  this->lambda = 0.0;
  //setting the random number generator
  // intialised in Initialising list :-> randomNumberGeneratorType = gsl_rng_default;//this is Mersenne Twister algorithm
  randomNumberGenerator = gsl_rng_alloc(randomNumberGeneratorType);
  gsl_rng_set(randomNumberGenerator,482023);//some number as seed-> this can be set with another random number/time or /dev/random /urandom
  // RandomNumberGenrator for Seed in growth for faces
  seedRandomNumberGenerator = gsl_rng_alloc(seedNumberGeneratorType);
  gsl_rng_set(seedRandomNumberGenerator,76500);//some number as seed-> this can be set with another random number/time or /dev/random /urandom

  // RANDOM NUMBER GENERATOR FOR CELL DIVISION
  cellDivisionRandomNumberGenerator = gsl_rng_alloc(randomNumberGeneratorType);
  //gsl_rng_set(cellDivisionRandomNumberGenerator, random_seed());//using /dev/urandom to seed this generator
  gsl_rng_set(cellDivisionRandomNumberGenerator, 123194);

  // RANDOM number generator for spontaneous (initial mean) curvature
  meanCurvatureRandomNumberGenerator = gsl_rng_alloc(randomNumberGeneratorType);
  gsl_rng_set(meanCurvatureRandomNumberGenerator, 902913);
  
  // RANDOM Number for RANDOM ANGLE GENERATOR OF FACE
  randomAngleGenerator = gsl_rng_alloc(randomNumberGeneratorType);
  gsl_rng_set(randomAngleGenerator, 33210);

  // Random Number (GAUSSIAN) for varying Angle of growth for FACE
  randomAngleGaussianVarianceGenerator = gsl_rng_alloc(randomNumberGeneratorType);
  gsl_rng_set(randomAngleGaussianVarianceGenerator, 64201);

  //calculating the square root of epsilon
  sqrtEpsilon = sqrt(std::numeric_limits<double>::epsilon());
  //sqrtEpsilon = pow(10.,-12.);
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
  gsl_rng_free(seedRandomNumberGenerator);
  gsl_rng_free(meanCurvatureRandomNumberGenerator);
  gsl_rng_free(randomAngleGenerator);
  gsl_rng_free(randomAngleGaussianVarianceGenerator);

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