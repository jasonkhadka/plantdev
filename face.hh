
/* ============================================================================
 * p2/cell/face.hh
 * ========================================================================= */

#ifndef faceINCLUDED
#define faceINCLUDED

#include "edge.hh"
#include "vertex.hh"

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
     
    //***************end added features****************************************//
  /* ----------Public instance methods ---------*/
 /* -- public instance methods ----------------------------------------- */
    public:
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
    * targetFormMatrix of this face (or cell in biological terms)
    * **for now this is jsut the form matrix of current shape
    */
    double targetFormMatrix[2][2] = {{0.,0.},{0.,0.}};
    /**
     * target Area: the ideal area of this cell (area from target Form)
     */
    double targetArea = 0.;
    /**
     * current form Matrix : the matrix that represents the 
     *  current form of the cell
     */
    double currentFormMatrix[2][2]= {{0.,0.},{0.,0.}};
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
  void setTempTargetFormMatrixCurrent();

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
   * get x centralised coordiante of this face
   */
   double getXCentralised();
   /**
   * set y centralised coordiante of this face
   */
   double getYCentralised();
    /**
   * set z centralised coordiante of this face
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

  friend class centralisedDerivative;
};
/* ----- inline instance methods --------------------------------*/
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

