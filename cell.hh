
/* ============================================================================
 * p2/cell/cell.hh
 * ========================================================================= */

#ifndef cellINCLUDED
#define cellINCLUDED

#include <iostream>
#include <set>
#include <math.h>
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

// SQRT_M_PI_5degree  = sqrt(pi/36)
#define SQRT_M_PI_5degree 0.295408975150919337883027913890190863799591576020397854702


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
   * calculating the volume of the structure 
   * using the non-convex generalised method
   */
  double getNonConvexVolume();
  /**
   * calculate the volume with cartesian (or actual) coordiantes
   */
  double getCartesianVolume();
  /**
   * calculate the volume with projected coordinate
   */
  double getProjectedCoordinateVolume();
  /**
   * Function to Initially set all the parameters to all the Faces & Vertices
   * the targetFormMatrix for faces is set to be current Form Matrix
   * to prepare for the further calculation of Energy and other quantities
   */
  void setInitialParameters();
  
  /**
   * Avearge Target Form matrix for the faces of this cell
   */
  double averageTargetFormMatrix[2][2];

  /**
   * Function to calculate average Target Form Matrix for all the cells
   * in this structure
   */
  void calculateAverageTFM();
  /**
   * Function to set average TFM of this structure on all cells
   * to be used to set initial conidition
   */
  void setAverageTFM();
  /**
   * Function to set all the parameters to all the Faces & Vertices
   * to prepare for the further calculation of Energy and other quantities
   */
  void setParameters();
  /**
   * Function to set the growth rate of all faces
   */
  void setFaceKappa();
  /**
   * Function to calculate the forces on the vertices
   */
  void calculateVertexForce();
  /**
   * Function to calculate stress and strain matrices in all the faces
   */
  void calculateStressStrain();
   /**
   * Function to calculate strain matrices in all the faces
   */
  void calculateStrain();
  /**
   * Function to get mean strain determinant
   */
  double getMeanStrainDeterminant();
  /**
   * Function to set areaMixed and Laplace-Beltrami operator for all vertices and the centroid of all faces
   * Also calculates the mean curvature for all centroid of faces and vertices
   */
  void setMeanCurvature();
  /**
   * Functioin to set the intiail mean curvature on all the vertex and face of this structure
   */
  void setInitialMeanCurvature();
    /**
   * Function to set the intiail mean curvature to a given value 
   * on all the vertex and face of this structure
   */
  void setInitialMeanCurvature(double);
  /**
   * get a random number to feed to face/vertex initial meancurvature
   * return : uniform random number in [-meanCurvatureWidth,meanCurvatureWidth)
   */
  double getRandomInitialMeanCurvature();
  /**
   * function to set random number for initial mean curvature for 
   * all vertices and faces
   */
  void setRandomInitialMeanCurvature();
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
  double pressure = 0.;
  /**
   * Alpha : coefficient of first term of Energy
   */
   double alpha = 1.;
   /**
   * Beta : coefficient of second term of Energy
   */
   double beta = 0.;
   /**
    * Kappa : the growth rate of faces
    */
   double kappa = 0.1;
   /**
    * Gamma : the pressure from underneath the L1 layer
    */
   double gamma = 1.;
   /**
    * Zeta : penalty for the bending of cells off the projection plane
    */
   double zeta = 0.;
   /**
    * Sigma : line tension on the walls
    */
   double sigma;
   /**
    * Growth Variation in Face growth rate
    */
   double growthvar = 0.5;
   /**
    * Eta : The feedback strength of stress in FACE GROWTH
    */
   double eta = 0.0;
   /**
    * Omega : the bending stiffness for Helfrich bending energy
    */
   double omega;
   /**
    * lambda : the magnitude of major eigenvalue for random matrix 
    * for growth reorganisation
    */
   double lambda;
   /**
    * initialStrain : The strain to assign as Initial Condition
    *                 In fraction of intial shape 
    */
   double initialStrain = 0.;
   /*
      Average Area of face in this Cell
    */
   double averageFaceArea;
   /*
      Surface Area of total tissue
    */
   double surfaceArea;
   /*
    * total Area Mixed calcualted
    */
   double totalAreaMixed;
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
   /**
    * initial volume
    */
   double initialVolume;
   /**
    * initial fourth term
    */
   double initialFourthTerm;
   /**
    * Helfrich bending energy
    */
   double bendingEnergy;
   // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
  //    Random Number Generator : seeded with *some* seed 
  //      Right now it is just a number i chose for testing*
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
   gsl_rng * randomNumberGenerator; //Random number generator
   const gsl_rng_type * randomNumberGeneratorType; //type of random number generator
   
   gsl_rng * cellDivisionRandomNumberGenerator; //random number generator for cell division
   const double gaussianWidth;
   
   // For the seed of faces
   gsl_rng * seedRandomNumberGenerator;
   const gsl_rng_type * seedNumberGeneratorType;

   // for spontaneous (intial) Mean curvature of surface
   double meanCurvatureWidth;
   gsl_rng * meanCurvatureRandomNumberGenerator;

    // randomized angle for random growth of face
   gsl_rng * randomAngleGenerator;

   // random number for changing angle of Growth of face
   gsl_rng * randomAngleGaussianVarianceGenerator;

public:
  /**getting/setting lambda */
  double getLambda();
  void setLambda(double);
  /**
   * Set/Get mean curvature width
   */
  double getMeanCurvatureWidth();
  void setMeanCurvatureWidth(double);

  /**
   * get random number (guassian) for varying 
   * random angle of growth for face
   */
  double getRandomAngleGaussianVariance();
  /**
   * get randomAngle for growth
   */
  double getRandomGrowthDirectionAngle();
  /*
    * get total area mixed
    */
  double getTotalAreaMixed();
  /**
   *                    Public members
   */
  // function to get initialVolume
  double getInitialVolume();
  // function to get initial Fourth term
  double getInitialFourthTerm();
  //Functions to calculate average area of face in this cell and get  it
  void setAverageFaceArea();
  double getAverageFaceArea();
  // functions to get the surface area of tissue
  double getSurfaceArea();
  /**
   * Getting guassianwidth of growthRandomNumber
   */
  double getGaussianWidth();
  /**
   * Setting & Getting intialStrain
   */
  void setInitialStrain(double);
  double getInitialStrain();
  /*
    *get Bending Energy
  */
  double getBendingEnergy();
  /*
    * calculate bending energy
  */
  void calculateBendingEnergy();
  /*
    * Calcualte bending force
  */
  void calculateBendingForce();
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
  void setDivisionThreshold(double);
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
   /**
    * Get Sum of all edge length of this structure
    */
    double getSumEdgeLength();
    /*
      * set radial orthoradial vector on faces
    */
    void setRadialOrthoradialVector(Face * primordialFace);
    /*
      * set radial and orthoradial stresses on faces
    */
    void setRadialOrthoradialStress();
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
   * function to set and get Omega
   */
  double getOmega();
  void setOmega(double);
  /*
    *set & get Eta()
   */
  double getEta();
  void setEta(double);
  // set and get sigma
  double getSigma();
  void setSigma(double);
  /**
   * Get uniform random number [0,1)
   */
  double getRandomNumber();
  /**
   * get RANDOM NUMBER for CELL DIVISION
   */
  double getCellDivisionRandomNumber();
  //Function to get Random Integer for seeds
  unsigned long int getSeedRandomInteger();

  //Function to set random number generator for seeding the face growth rate
  void setRandomNumberSeed(double newseed);
  /**
   * get the growth variation
   */
  double getGrowthVar();
  /**
   * set growth var for all the faces, copying from the cell->growthvar
   */
  void setFaceGrowthVar();
  /**
   * set the growth variation for this cell
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
inline double Cell::getLambda(){
  return this->lambda;
};
inline void Cell::setLambda(double tl){
  this->lambda = tl;
};
inline double Cell::getRandomAngleGaussianVariance(){//(generator, STDEV)
    return gsl_ran_gaussian_ziggurat(randomAngleGaussianVarianceGenerator, SQRT_M_PI_5degree);
};
inline double Cell::getRandomGrowthDirectionAngle(){
    return gsl_ran_flat(randomAngleGenerator,0.,M_PI);
};

inline double Cell::getMeanCurvatureWidth(){
    return meanCurvatureWidth;
};
inline void Cell::setMeanCurvatureWidth(double newcurve){
    this->meanCurvatureWidth = newcurve;
};
inline double Cell::getRandomInitialMeanCurvature(){
  return gsl_ran_flat(meanCurvatureRandomNumberGenerator,-1.*(this->meanCurvatureWidth),this->meanCurvatureWidth);
}
inline double Cell::getBendingEnergy(){
  return bendingEnergy;
}
inline double Cell::getTotalAreaMixed(){
  return this->totalAreaMixed;
}
inline double Cell::getInitialVolume(){
  return initialVolume;
}
inline double Cell::getInitialFourthTerm(){
  return initialFourthTerm;
}
inline double Cell::getSigma(){
  return this->sigma;
}
inline void Cell::setSigma(double newsigma){
  this->sigma = newsigma;
}
inline double Cell::getAverageFaceArea(){
  return this->averageFaceArea;
}
inline double Cell::getGaussianWidth(){
  return this->gaussianWidth;
}
inline double Cell::getSurfaceArea(){
  return this->surfaceArea;
}
inline void Cell::setInitialStrain(double initial){
  this->initialStrain = initial;
}
inline double Cell::getInitialStrain(){
  return this->initialStrain;
}
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

inline double Cell::getRandomNumber(){
      //using 
     //return gsl_ran_gaussian_ziggurat(randomNumberGenerator, this->gaussianWidth);
    return gsl_rng_uniform(randomNumberGenerator);
}
inline void Cell::setRandomNumberSeed(double newseed){
    gsl_rng_set(this->seedRandomNumberGenerator,newseed);
}
inline unsigned long int Cell::getSeedRandomInteger(){//get a Random integer (to act as a seed)
  return gsl_rng_get(seedRandomNumberGenerator);
}
inline double Cell::getCellDivisionRandomNumber(){
      return gsl_rng_uniform(cellDivisionRandomNumberGenerator);
}
inline double Cell::getKappa(){
  return kappa;
}
inline double Cell::getPressure()
{
  return pressure;
}
inline void Cell::setOmega(double tempo){
  this->omega = tempo;
}
inline double Cell::getOmega(){
  return this->omega;
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
inline double Cell::getEta(){
  return eta;
}
inline void Cell::setEta(double neweta){
    this->eta = neweta;
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

