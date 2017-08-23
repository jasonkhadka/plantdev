
/* ============================================================================
 * p2/cell/cell.hh
 * ========================================================================= */

#ifndef cellINCLUDED
#define cellINCLUDED

#include <iostream>
#include <set>

//random number generating
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>//for gaussian distribution

//for eigenvalue conputation
#include "./eigen/Eigen/Dense"
#include "./eigen/Eigen/Eigenvalues"

//#include "svl-1.5/include/svl/SVL.h"

#include "edge.hh"
#include "face.hh"
#include "vertex.hh"

class CellVertexIterator;
class CellFaceIterator;

/* ----------------------------------------------------------------------------
 * Cell
 * ------------------------------------------------------------------------- */

/*
 * An enclosed volume, bounded by a set of vertices and faces.
 *
 * Vertices   the vertices of the cell;
 *            all are nonnull
 * VertexIDs  an increasing sequence of positive integers used to number
 *            distinct vertices;
 *            all are positive
 * Faces      the faces of the cell;
 *            all are nonnull
 * FaceIDs    an increasing sequence of positive integers used to number
 *            distinct faces;
 *            all are positive
 */
class Cell
{

  /* -- public class methods ----------------------------------------------- */

  public:

  /*
   * Return a new, degenerate cell consisting of a single closed edge (loop),
   * a single vertex at the origin, and a pair of faces.
   * <- the new cell;
   *    will be nonnull
   */
  static Cell *make();

  /*
   * Return a new cell with the topology of a tetrahedron and all vertices at
   * the origin.
   * <- the new tetrahedron;
   *    will be nonnull
   */
  static Cell *makeTetrahedron();

  /*
   * Release the storage occupied by a given cell.
   * cell -> the cell to kill;
   *         must be nonnull
   */
  static void kill(Cell *cell);

  /* -- public instance methods (Euler operators) -------------------------- */

  /*
   * Use these methods to construct cells with guaranteed consistent topology.
   * Other means of modifying a cell can potentially produce bizarre results.
   */

  public:

  /*
   * Return a new edge formed by splitting a given vertex between a given pair
   * of faces.
   * A new vertex is introduced at the destination of the new edge.
   * The new edge has _left_ along its left and _right_ along its right.
   * vertex      -> the vertex to split to make the new edge;
   *                must be nonnull;
   *                must share an edge with both _left_ and _right_
   * left, right -> the faces adjacent to the new edge;
   *                must be nonnull;
   *                must share an edge with _vertex_
   * <- the new edge;
   *    will be nonnull
   */
  Edge *makeVertexEdge(Vertex *vertex, Face *left, Face *right);

  /*
   * Delete a given edge from this cell, along with its destination vertex.
   * edge -> the edge to kill;
   *         must be nonnull
   */
  void killVertexEdge(Edge *edge);

  /*
   * Return a new edge formed by splitting a given face through a given pair
   * of vertices.
   * A new face is introduced to the right of the new edge.
   * The new edge has _org_ as its origin and _dest_ as its destination.
   * face      -> the face to divide to make the new edge;
   *              must be nonnull;
   *              must have both _org_ and _dest_ on its perimiter
   * org, dest -> the vertices for the endpoints of the new edge;
   *              must be nonnull;
   *              must be located on the perimiter of _face_
   * <- the new edge;
   *    will be nonnull
   */
  Edge *makeFaceEdge(Face *face, Vertex *org, Vertex *dest);

  /*
   * Delete a given edge from this cell, along with its right face.
   * edge -> the edge to kill;
   *         must be nonnull
   */
  void killFaceEdge(Edge *edge);

  /* -- public instance methods -------------------------------------------- */

  public:

  /*
   * Return the number of vertices in this cell.
   * <- the number of vertices
   */
  unsigned int countVertices();
  /**
   * Return the ith vertex of this face. 
   * <- ith vertex on this face
   */
   Vertex *getIthVertex(unsigned int position);

  /*
   * Add a given vertex to this cell.
   * vertex -> the vertex to add;
   *           must be nonnull;
   *           must not be in the cell
   */
  void addVertex(Vertex *vertex);

  /*
   * Remove a given vertex from this cell.
   * vertex -> the vertex to remove;
   *           must be nonnull;
   *           must be in the cell
   */
  void removeVertex(Vertex *vertex);
  /**
   * Return the total number of vertex size in this cell
   * This is the max number of vertex cell can hold at current
   * Does not limit the further addition as vertexSize is increased with new addition
   */
  unsigned int totalVertexSize();
  /*
   * Return a new vertex ID.
   * <- a new vertex ID;
   *    will be positive
   */
  unsigned int makeVertexID();

  /*
   * Return the number of faces in this cell.
   * <- the number of faces
   */
  unsigned int countFaces();

  /*
   * Add a given face to this cell.
   * face -> the face to add;
   *         must be nonnull
   *         must not be in the cell
   */
  void addFace(Face *face);

  /*
   * Remove a given face from this cell.
   * face -> the face to remove;
   *         must be nonnull;
   *         must be in the cell
   */
  void removeFace(Face *face);

  /*
   * Return a new face ID.
   * <- a new face ID;
   *    will be positive
   */
  unsigned int makeFaceID();
  /**
   * calculate the volume of the structure
   */
  double getVolume();
  /**
   * calculate the volume of the structure
   * without Centroid
   */
  double getCartesianVolumeWOCentroid();
  /**
   * calculate the volume with cartesian (or actual) coordiantes
   */
  double getCartesianVolume();
  /**
   * Function to Initially set all the parameters to all the Faces & Vertices
   * the targetFormMatrix for faces is set to be current Form Matrix
   * to prepare for the further calculation of Energy and other quantities
   */
  void setInitialParameters();
  /**
   * Function to set all the parameters to all the Faces & Vertices
   * to prepare for the further calculation of Energy and other quantities
   */
  void setParameters();
  /**
   * Function to calculate the forces on the vertices
   */
  void calculateVertexForce();
  /**
   * Function to calculate stress and strain matrices in all the faces
   */
  void calculateStressStrain();
  /* -- protected instance methods ----------------------------------------- */

  protected:

  /*
   * Initialize this cell consisting of no vertices and no faces.
   */
  Cell();

  /*
   * Release the storage occupied by the contents of this cell.
   */
  ~Cell();

  /* -- private instance methods ------------------------------------------- */

  private:

  /*
   * Return the edge with a given origin vertex in the face orbit of a given
   * edge.
   * edge -> an edge of the orbit to look for the vertex in;
   *         must be nonnull
   * org  -> the origin vertex to look for;
   *         must be nonnull
   * <- the edge in the same face orbit as _edge_ with origin vertex _org_;
   *    null if not found
   */
  Edge *getOrbitOrg(Edge *edge, Vertex *org);

  /*
   * Set the origin of the vertex orbit of a given edge to a given vertex.
   * edge -> an edge of the orbit to set the origin vertex of;
   *         must be nonnull
   * org  -> the new origin vertex;
   *         must be nonnull
   */
  void setOrbitOrg(Edge *edge, Vertex *org);

  /*
   * Return the edge with a given left face in the vertex orbit of a given
   * edge.
   * edge -> an edge of the orbit to look for the face in;
   *         must be nonnull
   * left -> the left face to look for;
   *         must be nonnull
   * <- the edge in the same vertex orbit as _edge_ with left face _left_;
   *    null if not found
   */
  Edge *getOrbitLeft(Edge *edge, Face *left);

  /*
   * Set the left face of the face orbit of a given edge to a given face.
   * edge -> an edge of the orbit to set the left face of;
   *         must be nonnull
   * left -> the new left face;
   *         must be nonnull
   */
  void setOrbitLeft(Edge *edge, Face *left);
  /* -- private instance variables ----------------------------------------- */
  private:
  /**
   * Layer of cells on the Dome of this cell
   */
  int layer;
  /*
   * The vertices in this cell.
   * Nonnull.
   */
  Vertex **vertices;

  /*
   * The number of vertices in this cell.
   * Less than or equal to _vertexSize_.
   */
  unsigned int vertexCount;

  /*
   * The number of vertices allocated for this cell.
   * Greater than or equal to _vertexCount_.
   */
  unsigned int vertexSize;

  /*
   * The next unused vertex ID.
   */
  unsigned int vertexID;

  /*
   * The faces in this cell.
   * Nonnull.
   */
  Face **faces;

  /*
   * The number of faces in this cell.
   * Less than or equal to _faceSize_.
   */
  unsigned int faceCount;

  /*
   * The number of faces allocated for this cell.
   * Greater than or equal to _faceCount_.
   */
  unsigned int faceSize;

  /*
   * The next unused face ID.
   */
  unsigned int faceID;

  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
  // Parameters of Simulation 
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
  /**
   * Pressure of cells (faces) from inside the cells
   */
  double pressure = 1.;
  /**
   * Alpha : coefficient of first term of Energy
   */
   double alpha = 1.;
   /**
   * Beta : coefficient of second term of Energy
   */
   double beta = 1.;
   /**
    * Kappa : the growth rate of faces
    */
   double kappa = 1.;
   /**
    * Gamma : the pressure from underneath the L1 layer
    */
   double gamma = 1.;
   /**
    * Zeta : penalty for the bending of cells off the projection plane
    */
   double zeta = 1.;
   /**
    * Growth Variation in Face growth rate
    */
   double growthvar = 0.5;
   /**
    * Square root of Epsilon to use for finite diference stepsize
    * as recommended step size for finite difference is = sqrt(eps)*x for x != 0
    * assigned in constructor
    */
   double sqrtEpsilon;
   /**
    * Tolerance for the optimizer
    */
   double tolerance = 0.0001;
   /**
    * Initial step for simulation
    */
   double initialStep = 0.01;
   /**
    * function call counter
    */
   int functionCallCounter = 0;
   /**
    * cell division counter 
    */
   unsigned int divisionCounter;
   /**
    * division factor for cell division Threshold
    */
   double divisionFactor;
   /**
    * Angle threshold for decision for Cell is Convex or Not
    */
   double convexAngleThreshold;
   /**
    * Bending Threshold for cell for total Tissue
    */
   double bendingThreshold;
   // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
  //    Random Number Generator : seeded with *some* seed 
  //      Right now it is just a number i chose for testing*
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
   gsl_rng * randomNumberGenerator; //Random number generator
   const gsl_rng_type * randomNumberGeneratorType; //type of random number generator
   gsl_rng * cellDivisionRandomNumberGenerator; //random number generator for cell division
   const double gaussianWidth;

public:
  /**
   *                    Public members
   */
  /**
   * Function to set Bending threshold to find cell is convex or not
   */
  void setBendingThreshold(double);
  /**
   * Function to get back the Bending threshold to find cell is convex or not
   */
  double getBendingThreshold();
  /**
   * Function to set Angle threshold to find cell is convex or not
   */
  void setConvexAngleThreshold(double);
  /**
   * Function to get back the Angle threshold to find cell is convex or not
   */
  double getConvexAngleThreshold();
  /**
   * change division Factor for cell division threshold in this cell
   */
  void setDivisionFactor(double);
  /**
   * get the division factor
   */
  double getDivisionFactor();
  /**
   * register cell division
   */
  void countCellDivision();
  /**
   * return divisionCounter()
   */
  unsigned int getDivisionCounter();
   /**
    * Growth threshold for faces
    * it is initialised with constructor in Cell
    * */
   double thresholdMatrix[2][2] = {{0.,0.},{0.,0.}};
  /**
    * time step of calculation
    */
   double hstepsize = 1.;
   /**
    * Energy of Cell calculator : ADDED - P*V in the energy
    * where Volume is calculated by projected coordiantes
    */
   double getEnergy();
   /**
    * Energy of Cell calculator : -> with -P*V in the energy
    * where VOlume is calculated by cartesian coordiates by triangulating the vertices
    * without the Centroid
    */
   double getEnergyCartesianVolumeWOCentroid();
   /**
    * Energy of Cell calculator : -> with -P*V in the energy 
    * where VOlume is calcualted by Cartesian coordinates with Centroid
    */
   double getEnergyCartesianVolume();
   /**
    * Energy Function to calculate for Energy for use in numerical derivative
    */
   double getDerivativeEnergyCartesianVolume();
   
   /**
    * 
    */
   /**
    * First term of energy for this  cell
    */
   double getFirstTerm();
   /**
    * Second term of energy for this  cell
    */
   double getSecondTerm();
   /**
    * Third term of energy for this  cell
    */
   double getThirdTerm();
   /**
    * Fourth term of energy for this cell
    */
   double getFourthTerm();

  // **************************************************************** //
  // Public instance method //
public:
  /**
   * set the number of layer
   */
  void setLayer(int);
  /**
   * get the number of layer
   */
  int getLayer();

  /**
   * setting the pressure acting on the cells of this tissue
   * in the faces in quadedge 
   */
  void setPressure(double );
  /**
   * function to get the pressure of the cell
   */
  double getPressure();
  /**
   * setting the alpha parameter of energy
   */
  void setAlpha(double );
  /**
   * function to get the alpha of the cell
   */
  double getAlpha();

  /**
   * function to Set gamma value
   */
  void setGamma(double );
  /**
   * function to get Gamma value
   */
  double getGamma();
  /**
   * setting the beta parameter of energy
   */
  void setBeta(double );
  /**
   * function to get the beta of the cell
   */
  double getBeta();
  /**
   * setKappa()
   */
  void setKappa(double);
  /**
   * getKappa()
   */
  double getKappa();
  /**
   * get Zeta()
   */
  double getZeta();
  /**
   * set Zeta()
   */
  void setZeta(double);
  /**
   * Get a UNIFORM RANDOM NUMBER from the random number generator of this cell
   * in range [0,1) exclusive 1
   */
  double getRandomNumber();
  /**
   * get RANDOM NUMBER for CELL DIVISION
   */
  double getCellDivisionRandomNumber();
  /**
   * get the growth variation
   */
  double getGrowthVar();
  /**
   * set the growth variation
   */
  void setGrowthVar(double);
  /**
   * get Sqrt of Epsilon
   */
  double getSqrtEpsilon();
  /**
   * set the tolerance value
   */
  void setTolerance(double);
  /**
   * get the tolerance value
   */
  double getTolerance();
  /**
   * set the initial step
   */
  void setInitialStep(double);
  /**
   * get the initial step
   */
  double getInitialStep();
  /**
   * set function call counter
   */
  void setFunctionCallCounter(int);
  /**
   * get function call counter
   */
  int getFunctionCallCounter();
  /**
   * Set Cylindrical Coordinates from the Cartesian coordinates for all vertices
   */
  void setCylindrical();
  /**
   * Set Cartesian Coordinates from the Cylindrical coordinates for all vertices
   */
  void setCartesian();
  /**
   * Function to check if all the polygons (faces) in this Cell is convex or not
   * Return : 
   *        True : all are convex
   *        False : if one is not-convex
   */
  bool isConvex();
  // ********************************************************************* //
  /**
   * Function to relax the tissue
   * This uses NLOPT optmizer to relax the complete tissue
   * For more details see its implementation in cell.cc
   */
  int relax();
  // ********************************************************************* //
  /* -- friend classes ----------------------------------------------------- */

  friend class CellVertexIterator;
  friend class CellFaceIterator;
  friend class CentralisedDerivative;

};

/* -- inline instance methods ---------------------------------------------- */
inline double Cell::getBendingThreshold(){
  return this->bendingThreshold;
}
inline void Cell::setBendingThreshold(double newthreshold){
  this->bendingThreshold = newthreshold;
}
inline double Cell::getConvexAngleThreshold(){
  return this->convexAngleThreshold;
}
inline void Cell::setConvexAngleThreshold(double newangle){
  this->convexAngleThreshold = newangle;
}
inline double Cell::getDivisionFactor(){
  return this->divisionFactor;
}
inline void Cell::countCellDivision(){
  this->divisionCounter = this->divisionCounter + 1;
}
inline unsigned int Cell::getDivisionCounter(){
  return divisionCounter;
}
inline int Cell::getFunctionCallCounter(){
  return this->functionCallCounter;
}
inline void Cell::setFunctionCallCounter(int newvalue){
  this->functionCallCounter = newvalue;
}
inline void Cell::setInitialStep(double newinitial){
  this->initialStep = newinitial;
}
inline double Cell::getInitialStep(){
  return initialStep;
}
inline void Cell::setTolerance(double newtol){
  this->tolerance = newtol;
}
inline double Cell::getTolerance(){
  return tolerance;
}
inline void Cell::setLayer(int newlayer){
  this->layer = newlayer;
}
inline int Cell::getLayer(){
  return layer;
}
inline double Cell::getSqrtEpsilon(){
  return sqrtEpsilon;
}
inline double Cell::getGrowthVar(){
  return growthvar;
}
inline void Cell::setGrowthVar(double tempvar){
  growthvar = tempvar;
}
inline double Cell::getRandomNumber(){
      //using 
      return gsl_ran_gaussian_ziggurat(randomNumberGenerator, this->gaussianWidth);
}
inline double Cell::getCellDivisionRandomNumber(){
      return gsl_ran_gaussian_ziggurat(cellDivisionRandomNumberGenerator, this->gaussianWidth);
}
inline double Cell::getKappa(){
  return kappa;
}
inline void Cell::setKappa(double newkappa){
  this->kappa = newkappa;
}
inline double Cell::getPressure()
{
  return pressure;
}
inline double Cell::getAlpha()
{
  return alpha;
}
inline double Cell::getGamma()
{
  return gamma;
}
inline void Cell::setGamma(double newgamma)
{
  this->gamma = newgamma;
}
inline double Cell::getBeta()
{
  return beta;
}
inline void Cell::setZeta(double newzeta){
  this->zeta = newzeta;
}
inline double Cell::getZeta(){
  return zeta;
}
inline unsigned int Cell::countVertices()
{
  return vertexCount;
}

inline unsigned int Cell::makeVertexID()
{
  return vertexID++;
}

inline unsigned int Cell::countFaces()
{
  return faceCount;
}

inline unsigned int Cell::makeFaceID()
{
  return faceID++;
}
inline unsigned int Cell::totalVertexSize()
{
  return vertexSize;
}
inline Vertex *Cell::getIthVertex(unsigned int position)
{
  return vertices[position];
}
/* ----------------------------------------------------------------------------
 * CellVertexIterator
 * ------------------------------------------------------------------------- */

/*
 * Enumerates the vertices of a given cell in arbitrary order.
 */
class CellVertexIterator
{

  /* -- public instance methods -------------------------------------------- */

  public:

  /*
   * Initialize this vertex iterator over a given cell.
   * cell -> the cell to iterate the vertices of;
   *         must be nonnull
   */
  CellVertexIterator(Cell *cell)
  {
    this->cell  = cell;
    this->count = cell->vertexCount;
  }

  /*
   * Release the storage occupied by this vertex iterator.
   */
  ~CellVertexIterator()
  {
  }

  /*
   * Return the next vertex of this vertex iterator, if any.
   * <- the next vertex of this vertex iterator;
   *    null if none
   */
  Vertex *next()
  {
    // iterate the array in reverse order so that the current vertex can be
    // removed during iteration

    if (count<1)
      return 0;

    return cell->vertices[--count];
  }

  /* -- private instance variables ----------------------------------------- */

  private:

  /*
   * The cell whose vertices are being iterated.
   * Nonnull.
   */
  Cell *cell;

  /*
   * The number of vertices left to iterate.
   */
  unsigned int count;

};

/* ----------------------------------------------------------------------------
 * CellFaceIterator
 * ------------------------------------------------------------------------- */

/*
 * Enumerates the faces of a given cell in arbitrary order.
 */
class CellFaceIterator
{

  /* -- public instance methods -------------------------------------------- */

  public:

  /*
   * Initialize this face iterator over a given cell.
   * cell -> the cell to iterate the faces of;
   *         must be nonnull
   */
  CellFaceIterator(Cell *cell)
  {
    this->cell  = cell;
    this->count = cell->faceCount;
  }

  /*
   * Release the storage occupied by this face iterator.
   */
  ~CellFaceIterator()
  {
  }

  /*
   * Return the next face of this face iterator, if any.
   * <- the next face of this face iterator;
   *    null if none
   */
  Face *next()
  {
    // iterate the array in reverse order so that the current face can be
    // removed during iteration

    if (count<1)
      return 0;

    return cell->faces[--count];
  }

  /* -- private instance variables ----------------------------------------- */

  private:

  /*
   * The cell whose faces are being iterated.
   * Nonnull.
   */
  Cell *cell;

  /*
   * The number of faces left to iterate.
   */
  unsigned int count;

};

#endif /* #ifndef cellINCLUDED */

