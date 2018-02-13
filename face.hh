
/* ============================================================================
 * p2/cell/face.hh
 * ========================================================================= */

#ifndef faceINCLUDED
#define faceINCLUDED

#include "edge.hh"
#include "vertex.hh"
#include "./eigen/Eigen/Dense"
#include "./eigen/Eigen/Eigenvalues"
#include <algorithm>//for std::max
//random number generating
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>//for gaussian distribution and flat distribution
class Cell;

/* ----------------------------------------------------------------------------
 * Face
 * ------------------------------------------------------------------------- */

/*
 * A face of a cell, bounded by a set of directed edges.
 *
 * Cell   the cell that the face belongs to;
 *        nonnull
 * ID     the ID number assigned to the face by its cell (or the client);
 *        positive
 * data   generic data attached to the face by the client
 * Edges  the edges in the orbit of the face;
 *        all are nonnull
 */
class Face
{

  /* -- public class methods ----------------------------------------------- */

  public:

  /*
   * Return a new face no adjacent edges.
   * cell -> the cell that the face belongs to;
   *         must be nonnull
   * <- the new degenerate face;
   *    will be nonnull
   */
  static Face *make(Cell *cell);

  /*
   * Release the storage occupied by a given face.
   * face -> the face to kill;
   *         must be nonnull
   */
  static void kill(Face *face);

  /* -- public instance variables ------------------------------------------ */

  /*
   * The data associated with this edge.
   */
  const void *data;

  /* -- public instance methods -------------------------------------------- */

  public:

  /*
   * Return the cell for this face.
   * <- the cell that the face belongs to;
   *    will be nonnull
   */
  Cell *getCell();

  /*
   * Return the ID of this face.
   * <- the ID of this face;
   *    will be positive
   */
  unsigned int getID();

  /*
   * Change the ID of this face.
   * id -> the new id for this face;
   *       must be positive
   */
  void setID(unsigned int id);

  /*
   * Return an arbitrary adjacent edge for this face.
   * <- an edge that is adjacent to this face;
   *    null if degenerate
   */
  Edge *getEdge();

  /*
   * Add a given adjacent edge to this face.
   * edge -> an edge that is adjacent to this face;
   *         must be nonnull
   */
  void addEdge(Edge *edge);
    
  /*
   * Remove a given adjacent from this face.
   * edge -> an edge that is no longer adjacent to this face;
   *         must be nonnull
   */
  void removeEdge(Edge *edge);
  
  //********************added features ****************************//
  /*
   * Return the first vertex of this face.
   * <- a vertex that is adjacent to this face;
   *    null if degenerate
   */
  //Vertex *getVertex();
  /**
   * Return the ith vertex of this face. 
   * <- ith vertex on this face
   */
  //Vertex *getIthVertex(unsigned int position);
  /**
   * Add a given vertex to this face
   * vertex -> a vertex to be added to this face
   *           is one of the vertex forming the face
   *           must be nonnull
   */
   //void addVertex(Vertex *vertex);

   /**
   * Remove a given vertex from this face
   * vertex -> a vertex to be removed to this face
   *           is one of the vertex forming the face
   *           must be nonnull
   */
   
   //void removeVertex(Vertex *vertex);
   //********************end added features ****************************//
  /* -- protected instance methods ----------------------------------------- */

  protected:

  /*
   * Initialize this face with no adjacent edges.
   * cell -> the cell that this face belongs to;
   *         must be nonnull
   */
  Face(Cell *cell);

  /*
   * Release the storage occupied by the contents of this face.
   */
  ~Face();

  /* -- private instance variables ----------------------------------------- */

  private:

  /*
   * The cell that this face belongs to.
   * Nonnull.
   */
  Cell *cell;
  

  /*
   * The ID of this face.
   * Positive.
   */
  unsigned int id;
  /**
   * Rank of the Face in Tissue
   * 1 is the top most cell, 
   * N is the bottom most cells
   */
  unsigned int faceRank;
  /*
   * An arbitrary adjacent edge to this face.
   * Null if degenerate.
   */
  Edge *edge;
  /**
   * x centroid coordinate 
   * of this face or (cell in biology) 
   */
  double xCentralised;
  /**
   * y centroid coordinate 
   * of this face or (cell in biology) 
   */
  double yCentralised;
    /**
   * Z centroid coordinate 
   * of this face or (cell in biology) 
   */
  double zCentralised;
  /**
   * Angle between the projected average normal and the cartesian X- axis. 
   */
   double angleOfTilt;
  /*
   * Area of face
   * Positive
   * stored as Area = 1/2* \sum_{edges} A_k
   */
   double areaOfFace;

  //***************added features*******************************************//
   //defing the type for 2d double array that would be used as matrix for Form Matrix
   /**
    * The vertices in this face
    * nonnull
    */
    Vertex **vertices;
    /*
     * The number of vertices in this face.
     * Less than or equal to _vertexSize_.
     */
    int vertexCount;

    /*
     * The number of vertices allocated for this cell.
     * Greater than or equal to _vertexCount_.
     */
    unsigned int vertexSize;
    /**
     * normal of this face [n^c_x,n^c_y,n^c_z]
     * array stores x,y,z component of the normal
     */
     double normal[3];
     /**
       * normal before normalisation
       */
      double normaltilde[3];
     /**
      * pi_vector of this face (refer to the derivative of coordinate derivation)
      *[ pi_x^p, pi_y^p, pi_z^p]
      */
      double pivector[3]; 
    /**
     * unit vector in X direction
     */
     double unitx[3];
     /**
     * unit vector in Y direction
     */
     double unity[3];
     /**
     * unit vector in Z direction
     */
     double unitz[3]; 
     /**
      * Unit Vector with respect to cartesian coordinate
      */
     double nonCentralisedUnitx[3];
     double nonCentralisedUnity[3];
     double nonCentralisedUnitz[3];
     /**
      * Stress Eigen Directions & Eigen values
      */
     double stressEigenVector1[3];
     double stressEigenVector2[3];
     double stressEigenValue1, stressEigenValue2;
    /**
      * Strain Eigen Directions & Eigen values
      */
     double strainEigenVector1[3];
     double strainEigenVector2[3];
     double strainEigenValue1, strainEigenValue2;
     /**
      * Normal Force on this Face
      */
     double normalForce[3];

     // ****************** ENERGY VALUES ************************************** //
     /**
      * First term value of Energy
      */
     double firstTerm = 0;
     /**
      * Second term value of Energy
      */
     double secondTerm = 0;
     /**
      * third term value of Energy
      */
     double thirdTerm = 0;
     /**
      * Total energy of The face 
      */
     double energy = 0;
       /**
       * Variable that notes if the Face is in Dome or in Cylinder
       */
       bool domePosition;
       /**
        * Alpha for this face
        */
       double alpha =0.;
      /**
       * Determinant of Strain
      */
      double strainDeterminant;
      /**
       * Trace of Strain
      */
      double strainTrace;
      
      /*
        * The growth of of this Face
      */
     double kappa;
     /*
        * Growth Variance for this face
      */
     double growthVar = 0.5;
     /* 
      * Variable to store last growth rate : For plotting the growth rate
     */
     double lastGrowthRate;
     //Volume of this face
     double faceVolume;
     // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
      //    Random Number Generator : seeded with *some* seed 
     //      Right now it is just a number i chose for testing*
     // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
   gsl_rng * randomNumberGenerator; //Random number generator
   const gsl_rng_type * randomNumberGeneratorType; //type of random number generator
   gsl_rng * cellDivisionRandomNumberGenerator; //random number generator for cell division
   const double gaussianWidth;
    //***************end added features****************************************//
  /* ----------Public instance methods ---------*/
 /* -- public instance methods ----------------------------------------- */
    public:
   //Function to get faceVoume
    double getFaceVolume();
    void setFaceVolume(double);
  // Get the random Number //
  double getGrowthRandomNumber();
  /**
   * get RANDOM NUMBER for CELL DIVISION
   */
  double getCellDivisionRandomNumber();
      /**
       * set/get faceRank
       */
    unsigned int getFaceRank();
    void setFaceRank(unsigned int);

    /**
     * Functions to manipulate Kappa
     */
    double getKappa();
    void setKappa(double);

    /**setting and getting growth var**/
     double getGrowthVar();
     void setGrowthVar(double);
    /**
     * get Last growthrate
     */
    double getLastGrowthRate();
    
      /**function to manupulate alpha**/
    double getAlpha();
    void setAlpha(double);
    /**
     * function to get stess Eigen Vectors & values
     */
    double * getStressEigenVector1();
    double * getStressEigenVector2();
    double * getNormalForce();
    double getStressEigenValue1();
    double getStressEigenValue2();    
    /**
     * function to get strain Eigen Vectors & values
     */
    
    double * getStrainEigenVector1();
    double * getStrainEigenVector2();
    double getStrainEigenValue1();
    double getStrainEigenValue2();    
    double getStrainValue(int, int);
    /**
     * setting the number of vertices in this face
     */
    void setVertexCount();
    /**
     * getting the number of vertices in this face
     */
     int getVertexCount();
    /**
     *function to set the normal
     */
     void setNormal(double *);
    /**
     * function to return the normal
     */
     double * getNormal();
     /**
      * function to set the normal tilde
      */
      void setNormalTilde(double *);
      /**
       * function to get the normal tilde
       */
      double * getNormalTilde();
      
     /**
     *function to set the normal
     */
     void setPivector(double *);
    /**
     * function to return the normal
     */
     double * getPivector();
     /**
      * function to set unit vector in X direction 
      */
      void setUnitx(double *);
      /**
       * function to get unit vector in X direction
       */
       double * getUnitx();
       /**
      * function to set unit vector in Y direction 
      */
      void setUnity(double *);
      /**
       * function to get unit vector in Y direction
       */
       double * getUnity();
       /**
      * function to set unit vector in Z direction 
      */
      void setUnitz(double *);
      /**
       * function to get unit vector in Z direction
       */
       double * getUnitz();
        /**
         * function to get non-centralised unitvector in x y z
         */
        double * getNoncentralisedUnitx();
        double * getNoncentralisedUnity();
        double * getNoncentralisedUnitz();
        void setNoncentralisedunitx(double*);
        void setNoncentralisedunity(double*);
        void setNoncentralisedunitz(double*);
  // ***************************************** //
  // TARGET FORM MATRIX //
  // ***************************************** //
  /**
   * target form matrix growth function
   * grows the matrix with 
   * stocahstic and stress driven growth
   */
   void grow();
   /**
    * target form matrix grows proportionally to the itself
    * without any mehcanics just pure inflate growth
    */
   void inflatedGrow();
   /**
   * target form matrix growth function with base grow proportional
   * to the targetformmatrix
   * grows with the Feedback from the stress acting on the cell
   * (mimicing microtubles re-enforcing)
   */
   void feedbackInflatedGrow();
   /**
   * target form matrix growth function with base grow proportional
   * to the strain
   * grows with the Feedback from the stress acting on the cell
   * (mimicing microtubles re-enforcing)
   */
   void feedbackStrainGrow();
   /**
    * The threshold for division
   */
   double divisionThreshold;
   /**
    * the factor for division threshold 
    * divisionThreshold = divisionFactor*IntialArea;
    * initial set as 10%  extra growth
    */
   double divisionFactor;
   /**
    * Function to set divisionFactor
    */
   void setDivisionFactor(double);
   double getDivisionThreshold();
   /**
    * set the threshold for division as compared to Current Area of Face
   */
   void setDivisionThreshold();
   /**
    * this function is to set the division threshold to the given value
    */
   void setDivisionThreshold(double);

  /**
     * function to divide this face in a random axis
     */
     void divideRandom();
  /**
   * function to set the Form Matrix after the division 
   * need to set the form matrix of 
   * 1. the Daugther cells of division 
   * 2. Neighbouring cells as those also got new vertex
   */
    void setDivideFormMatrix();
    /**
   * Function to change the position of Face to Cylinder
   */
  void positionToCylinder();
  /**
   * Function  to change the position of Face to Dome
   */
  void positionToDome();
  /**
   * Function to return Dome Position
   *  True : Face is in Dome
   *  False : Face is in Cylindrical Flanks
   */
  bool getDomePosition();
  /**
   * Function to test if this polygon is Convex or not
   * Convex :- A convex polygon is a simple polygon (not self-intersecting) 
   *           in which no line segment between two points on the boundary 
   *           ever goes outside the polygon. Equivalently, it is a simple 
   *           polygon whose interior is a convex set.[1] In a convex polygon,
   *           all interior angles are less than or equal to 180 degrees, 
   *           while in a strictly convex polygon all interior angles are 
   *           strictly less than 180 degrees.
   *           
   * return :- True : if this polygon is Convex
   *           False : if this polygon is non-Convex
   */
  bool isConvex();

   /**
    * targetFormMatrix of this face (or cell in biological terms)
    * **for now this is jsut the form matrix of current shape
    */
    double targetFormMatrix[2][2] = {{0.,0.},{0.,0.}};
    /**
     * target Area: the ideal area of this cell (area from target Form)
     */
    double targetArea = 1.;
    /**
     * current form Matrix : the matrix that represents the 
     *  current form of the cell
     */
    double currentFormMatrix[2][2]= {{0.,0.},{0.,0.}};
    /**
     * Intrinsic Stress Matrix of Face 
     */
    Eigen::Matrix2d stress;
    /**
     * Intrinsic Strain Matrix of Face
     */
    Eigen::Matrix2d strain;
    /**
     * Function to calculate intrinsic stress matrix
     */
    void calculateStress();
    /**
     * Function to calculate the intrinsic strain matrix
     */
    void calculateStrain();
    /**
     * get Strain Determinant and trace
     */
    double getStrainDeterminant();
    double getStrainTrace();
    /**
     * After cell division,  To update the form matrix after the division
     * Old Mu Matrix : the matrix that represents the 
     *  Mu matrix (difference of currentForm and targetForm) of the cell before cell division
     */
    double oldMuMatrix[2][2]= {{0.,0.},{0.,0.}};
    /**
     * Simular to oldMuMatrix we need trace of current form matrix
     * for targetFormMatrix for cells after divison
     */
    double oldcurrentFormMatrixTrace = 0.;
    
    /**
     * function that prints the targetform matrix
     * this function ideally should be replaced by somethign that 
     * just returnss the target fomr matrix
     */
    void printTargetFormMatrix();
    /**
     * trace of targetFormMatrix : will be updated as soon as
     *  new targetFormMatrix is calculated
     */
     double traceSquaredTargetFormMatrix = 0.;
    /**
     * values of mu1, mu2, mu3, mu4
     * Needed to calculate the analytic jacobian
     */
     double mu1= 0.;
     double mu2 = 0.;
     double mu3 = 0.;
     double mu4 = 0.;
  /**
   * function to set the value of mu1, mu2, mu3, mu4
   */
   void setMu();
   /**
    * get mu1
    */
  double getMu1();
    /**
    * get mu2
    */
  double getMu2();
    /**
    * get mu3
    */
  double getMu3();
    /**
    * get mu4
    */
  double getMu4();
  /**
   * function to set the angle of tilt
   */
   void setAngleOfTilt(double);
   /**
    * function to get the angle of tilt
    */
  double getAngleOfTilt();
  /**
   * the function to set the targetFormMatrix
   * right now just the form matrix of current shape
   * in Future : this should include the time update at each growth step
   */
  void setTargetFormMatrix();
  /**
   * set the target form matrix to be near to the Identity form
   */
  void setTempTargetFormMatrixIdentity();
  /**
   * set the target form matrix to be near to the current form
   */
  void setInitialTargetFormMatrixCurrent();
  /**
   * set the target form matrix to be bit narrower in cylindrical direction
   */
  void setSkewedTargetFormMatrixCurrent();

  /**function to set the trace of targetformmatrix
  */
  void setTraceSquaredTargetFormMatrix();
  /**
   * getTraceTargetFormMatrix
   */
  double getTraceSquaredTargetFormMatrix();
  /**
   * set centralised Coordinates of this face
   */
   void setCentralisedCoordinate(double, double,double);
   /**
    * Function to set new projected X & Y coordinate of vertices of this face
    */
    void setProjectedCoordinate();
    /**
     * get Target Area for calculation of Area Deformation;
     */
    double getTargetArea();
    /**
     * set target area; for calculation of Area Deformation
     */
    void setTargetArea(double);
    /**
     * get target form matrix determinant
     */
    double getTargetFormMatrixDeterminant();

    /**
     * Calcuate Projected Coordinate
     * Function that calcualtes the projection of given vertex
     * on the given face
     * Input : Vertex, Face
     * return : X, Y, Z projection
     */
    
  /**
   * get x coordiante of centroid of this face
   */
   double getXCentralised();
   /**
   * set y coordiante of centroid of this face
   */
   double getYCentralised();
    /**
   * set z coordiante of centroid of this face
   */
   double getZCentralised();

  /*
   * set area of face 
   * areaOfFace
   * no need for input as this function directly goes through all edges
   * and sums up the Ak values
   * and assigns the area value to areaOfFace variable
   */
  void setAreaOfFace();
   /*
    * get area of face
    * this function returns the value of area
    */
  double getAreaOfFace();
  /*
   * Return the number of vertices in this cell.
   * <- the number of vertices
   */
  unsigned int countVertices();
  // ****************** ENERGY VALUES Setters/ Getters ************************************** //
     /**
      * get First term value of Energy
      */
     double getFirstTerm();
     /**
      * get Second term value of Energy
      */
     double getSecondTerm();
     /**
      *  get third term value of Energy
      */
     double getThirdTerm();
     /**
      * get Energy values
      */
     double getEnergy();
      /**
      * Set the values of Energy terms :- First, Second, Third term value of Energy 
      */
     void setEnergyTerms();


  // **********************************************************************// 
  /* -- friend classes ----------------------------------------------------- */

  friend class CentralisedDerivative;
};
/* ----- inline instance methods --------------------------------*/
inline double Face::getFaceVolume(){
  return this->faceVolume;
}
inline void  Face::setFaceVolume(double newvol){
  this->faceVolume = newvol;
};
inline double Face::getGrowthVar(){
  return growthVar;
}
inline void Face::setGrowthVar(double newvar){
  this->growthVar = newvar;
}
inline double Face::getLastGrowthRate(){
  return lastGrowthRate;
}
inline double Face::getGrowthRandomNumber(){//getting random number between kappa-growthvar to kappa+growthvar
    return gsl_ran_flat(randomNumberGenerator, std::max(0.,kappa-growthVar),kappa+growthVar);
}
inline double Face::getCellDivisionRandomNumber(){
      return gsl_rng_uniform(cellDivisionRandomNumberGenerator);
}
inline double Face::getKappa(){
  return kappa;
}
inline double Face::getStrainValue(int a, int b){
  return this->strain(a,b);
}
inline double Face::getStrainDeterminant(){
  return this->strainDeterminant;  
};
inline double Face::getStrainTrace(){
  return this->strainTrace;
}
inline unsigned int Face::getFaceRank(){
      return this->faceRank;
}
inline void Face::setFaceRank(unsigned int newrank){
      this->faceRank = newrank;
}
inline double Face::getAlpha(){
  return this->alpha;
}
inline void Face::setAlpha(double newalpha){
  this->alpha = newalpha;
}
inline double * Face::getNormalForce(){
    double * pnt = this->normalForce;
    return pnt;
}
inline double Face::getStressEigenValue1(){
  return this->stressEigenValue1;
};
inline double Face::getStressEigenValue2(){
  return this->stressEigenValue2;
};
inline double Face::getStrainEigenValue1(){
  return this->strainEigenValue1;
};
inline double Face::getStrainEigenValue2(){
  return this->strainEigenValue2;
};
inline bool Face::getDomePosition(){
  return this->domePosition;
}
inline void Face::positionToDome(){
  this->domePosition = true;
}
inline void Face::positionToCylinder(){
  this->domePosition = false;
}
inline void Face::setDivisionFactor(double newfactor){
    this->divisionFactor = newfactor;
}
inline void Face::setDivisionThreshold(double newthreshold){
    this->divisionThreshold = newthreshold;
}
inline double Face::getDivisionThreshold(){
    return this->divisionThreshold;
}
inline double Face::getFirstTerm(){
    return firstTerm;
}
inline double Face::getSecondTerm(){
    return secondTerm;
}
inline double Face::getThirdTerm(){
    return thirdTerm;
}
inline double Face::getEnergy(){
    return energy;
}
inline void Face::setAngleOfTilt(double angle)
{
    angleOfTilt = angle;
}
inline double Face::getAngleOfTilt()
{
    return angleOfTilt;
}
inline unsigned int Face::countVertices()
{
  return vertexCount;
}

/* -- inline instance methods ---------------------------------------------- */
inline void Face::setTargetArea(double temparea)
{
  targetArea = temparea;
}
inline double Face::getTargetArea()
{
  return targetArea;
}
inline Cell *Face::getCell()
{
  return cell;
}

inline unsigned int Face::getID()
{
  return id;
}

inline Edge *Face::getEdge()
{
  return edge;
}
//************added features***************************************//
/*
inline Vertex *Face::getVertex()
{
  return vertices[0];
}
inline Vertex *Face::getIthVertex(unsigned int position)
{
  return vertices[position];
}
*/
inline int Face::getVertexCount()
{
  return vertexCount;
}
inline double Face::getXCentralised()
{
  return xCentralised;
}
inline double Face::getYCentralised()
{
  return yCentralised;
}
inline double Face::getZCentralised()
{
  return zCentralised;
}
inline double Face::getMu1()
{
  return mu1;
}
inline double Face::getMu2()
{
  return mu2;
}
inline double Face::getMu3()
{
  return mu3;
}
inline double Face::getMu4()
{
  return mu4;
}
inline double Face::getTraceSquaredTargetFormMatrix()
{
  return traceSquaredTargetFormMatrix;
}
inline double Face::getTargetFormMatrixDeterminant()
{
  return targetFormMatrix[0][0]*targetFormMatrix[1][1]-targetFormMatrix[0][1]*targetFormMatrix[1][0];
}
//************ end added features **********************************//
/* ----------------------------------------------------------------------------
 * FaceEdgeIterator
 * ------------------------------------------------------------------------- */

/*
 * Enumerates the bounding edges of a given face in counterclockwise order.
 */
class FaceEdgeIterator
{

  /* -- public instance methods -------------------------------------------- */

  public:

  /*
   * Initialize this edge iterator over a given face.
   * face -> the face to iterate the edges of;
   *         must be nonnull
   */
  FaceEdgeIterator(Face *face)
  {
    // pick an arbitrary edge in the face orbit

    start = face->getEdge();
    edge  = start;
  }

  /*
   * Release the storage occupied by this edge iterator.
   */
  ~FaceEdgeIterator()
  {
  }

  /*
   * Return the next edge of this edge iterator, if any.
   * <- the next edge of this edge iterator;
   *    null if none
   */
  Edge *next()
  {
    // check for degeneracy or exhausted iteration

    Edge *current = edge;

    if (current==0)
	return 0;

    // get the next edge in the left orbit of the face, but return the current
    // edge
    // reset to null if we've come back to the start

    Edge *next = current->Lnext();

    edge = next!=start ? next : 0;

    return current;
  }

  /* -- private instance variables ----------------------------------------- */

  private:

  /*
   * The first edge to be iterated.
   * Nonnull.
   */
  Edge *start;

  /*
   * The next edge to be iterated.
   * Null if exhausted.
   */
  Edge *edge;

};

#endif /* #ifndef faceINCLUDED */

