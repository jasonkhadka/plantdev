
/* ============================================================================
 * p2/cell/vertex.hh
 * ========================================================================= */

#ifndef vertexINCLUDED
#define vertexINCLUDED
#include <map>
//#include "svl-1.5/include/svl/Vec3.h"
#include "svl-1.5/include/svl/SVL.h"
#include "edge.hh"

class Cell;

/* ----------------------------------------------------------------------------
 * Vertex
 * ------------------------------------------------------------------------- */

/*
 * A vertex of a cell, with an outgoing set of directed edges.
 *
 * pos    the three-dimensional position of the vertex
 * Cell   the cell that the vertex belongs to;
 *        nonnull
 * ID     the ID number assigned to the vertex by its cell (or the client);
 *        positive
 * data   generic data attached to the vertex by the client
 * Edges  the edges in the orbit of the vertex;
 *        all are nonnull
 */
class Vertex
{

  /* -- private class member ----------------------------------------------- */
  /*coordinates of vertex
   *all initialised to origin (0,0,0)
   */
  double Xcoordinate;
  double Ycoordinate;
  double Zcoordinate;
    /*
    * Associative array for storing the projected X coordinate
    * Key : Face ID -- unsigned int
    * Value : the projected X coordinate -- double
    * ProjectedXcoordinate -> the projected X coordinate of this vertex
    */
   std::map<unsigned int, double> ProjectedXcoordinate;
   /*
    * Associative array for storing the projected Y coordinate
    * Key : Face ID -- unsigned int
    * Value : the projected Y coordinate -- double
    * ProjectedYcoordinate -> the projected Y coordinate of this vertex
    */
   std::map<unsigned int, double> ProjectedYcoordinate;
   //========================Derivatives===============================//
   /**
    * X derivative of area of face with respect to this vertex
    */
    std::map<unsigned int, double> areaXDerivative;
    /**
    * Y derivative of area of face with respect to this vertex
    */
    std::map<unsigned int, double> areaYDerivative;
    /**
    * X derivative of Ak of face with respect to this vertex
    */
    std::map<unsigned int, double> AkXDerivative;
    /**
    * Y derivative of Ak of face with respect to this vertex
    */
    std::map<unsigned int, double> AkYDerivative;
    // %%%%%%%%%%%%%%%%%%% X derivative %%%%%%%%%%%%%%%%%%%% //
    /**
     * x derivative of mu1
     */
    std::map<unsigned int, double>mu1XDerivative;
    /**
     * y derivative of mu4
     */
    std::map<unsigned int, double>mu4XDerivative;
    
    /**
     * X derivative of mu1_squared
     */
     std::map<unsigned int, double>mu1SquaredXDerivative;
     /**
     * X derivative of mu2_squared
     */
     std::map<unsigned int, double>mu2SquaredXDerivative;
     /**
     * X derivative of mu3_squared
     */
     std::map<unsigned int, double>mu3SquaredXDerivative;
          /**
     * X derivative of mu4_squared
     */
     std::map<unsigned int, double>mu4SquaredXDerivative;
     /**
      * First term X derivative
      */
     double firstTermXDerivative;
     /**
      * Second term X derivative
      */
     double secondTermXDerivative;
     /**
     * third term X derivative
     */
     double thirdTermXDerivative;
     // %%%%%%%%%%%%%%%%%%% Y derivative %%%%%%%%%%%%%%%%%%% //
     /**
     * x derivative of mu1
     */
    std::map<unsigned int, double>mu1YDerivative;
    /**
     * y derivative of mu4
     */
    std::map<unsigned int, double>mu4YDerivative;
     /**
     * Y derivative of mu1_squared
     */
     std::map<unsigned int, double>mu1SquaredYDerivative;
     /**
     * Y derivative of mu2_squared
     */
     std::map<unsigned int, double>mu2SquaredYDerivative;
     /**
     * Y derivative of mu3_squared
     */
     std::map<unsigned int, double>mu3SquaredYDerivative;
     /**
     * Y derivative of mu4_squared
     */
     std::map<unsigned int, double>mu4SquaredYDerivative;
     /**
      * First term derivative
      */
     double firstTermYDerivative;
     /**
      * Second term Y derivative
      */
     double secondTermYDerivative;
     /**
      * third term Y derivative
      */
     double thirdTermYDerivative;
    //========================AK and Functions===============================//
/*
    * Associative array for storing the value of function1 [f1(k)]
    * Key : Face ID -- unsigned int
    * Value : function1 valued at origin and destination of this edge -- double
    * Function1 [f1(k)] -> y_k^2 + y_k y_{k+1} + y_{k+1}^2 
    */
   std::map<unsigned int, double> Function1;
   /*
    * Associative array for storing the value of function2 [f2(k)]
    * Key : Face ID -- unsigned int
    * Value : function2 valued at origin and destination of this edge -- double
    * Function2 [f2(k)] ->  x_k*y_{k+1}+2x_k*y_k+2x_{k+1}*y_{k+1}+x_{k+1}y_k
    */
   std::map<unsigned int, double> Function2;
   /*
    * Associative array for storing the value of function3 [f3(k)]
    * Key : Face ID -- unsigned int
    * Value : function3 valued at origin and destination of this edge -- double
    * Function3 [f3(k)] -> x_k^2 + x_k x_{k+1} + x_{k+1}^2 
    */
   std::map<unsigned int, double> Function3;
   // *************************************************************//
   /*
    * Associative array for storing the value of Ak [Ak(k)]
    * Key : Face ID -- unsigned int
    * Value : Ak valued at origin and destination of this edge -- double
    * Ak[ak(k)] -> x_k*y_{k+1} - x_{k+1}*y_k 
    */
   std::map<unsigned int, double> Ak;

  /* -- public class methods ----------------------------------------------- */
  public:
  /*
   * Return a new vertex at the origin with no outgoing edges.
   * cell -> the cell that the vertex belongs to;
   *         must be nonnull
   * <- the new vertex;
   *    will be nonnull
   */
  static Vertex *make(Cell *cell);

  /*
   * Release the storage occupied by a given vertex.
   * vertex -> the vertex to kill;
   *           must be nonnull
   */
  static void kill(Vertex *vertex);

  /* -- public instance variables ------------------------------------------ */

  /*
   * The three-dimensional position of this vertex.
   */
  Vec3 pos;

  /*
   * The data associated with this edge.
   */
  const void *data;

  /* -- public instance methods -------------------------------------------- */

  public:

  /*
   * Return the cell for this vertex.
   * <- the cell that the vertex belongs to;
   *    will be nonnull
   */
  Cell *getCell();
  /*
   * Return the x coordinate of this vertex.
   * xcoordinate -> the x coordinate of this vertex;
   *    will be double
   */
   double getXcoordinate();
   /*
   * Return the y coordinate of this vertex.
   * ycoordinate -> the y coordinate of this vertex;
   *    will be double
   */
   double getYcoordinate();
   /*
   * Return the z coordinate of this vertex.
   * zcoordinate -> the z coordinate of this vertex;
   *    will be double
   */
   double getZcoordinate();
   /*
   * assign the x coordinate of this vertex.
   * xcoordinate -> the x coordinate of this vertex;
   *    will be double
   */
   void setXcoordinate(double xcoord);
   /*
   * assign the y coordinate of this vertex.
   * ycoordinate -> the y coordinate of this vertex;
   *    will be double
   */
   void setYcoordinate(double ycoord);
   /*
   * assign the z coordinate of this vertex.
   * zcoordinate -> the z coordinate of this vertex;
   *    will be double
   */
   void setZcoordinate(double ycoord);
  /*
   * Return the ID of this vertex.
   * <- the ID of this vertex;
   *    will be positive
   */
  unsigned int getID();
  /*
   * Change the ID of this vertex.
   * id -> the new id for this vertex;
   *       must be positive
   */
  void setID(unsigned int id);
  /*
   * Function to insert new elements
   * in associated array of ProjectedXcoordinate
   * faceid - unsigned int 
   * xcood - double
   */
   void insertProjectedXcoordinate(unsigned int faceid, double xcood);
   /*
   * Function to insert new elements
   * in associated array of ProjectedYcoordinate
   * faceid - unsigned int 
   * ycood - double
   */
   void insertProjectedYcoordinate(unsigned int faceid, double ycood);
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
  /*
   * Function to get projected xcoordinate from the faceid
   * in associated array of ProjectedYcoordinate
   * faceid - unsigned int 
   * return - double
   */
   double getProjectedXcoordinate(unsigned int faceid);
   /*
   * Function to get projected Ycoordinate from the faceid
   * in associated array of ProjectedYcoordinate
   * faceid - unsigned int 
   * return - double
   */
   double getProjectedYcoordinate(unsigned int faceid);
   //====================Derivatives get/set functions===========================//
   /**
    * Function to set the areaDerivative with respect to x and y given vertex
    */
    void setAreaDerivative();
    /**
    * Function to get the area Derivative with respect to X for the given vertex
    * faceid : the id of the face on which the derivative is to be returned
    */
    double getAreaXDerivative(unsigned int faceid);
    /**
    * Function to get the area Derivative with respect to Y for the given vertex
    * faceid : the id of the face on which the derivative is to be returned
    */
    double getAreaYDerivative(unsigned int faceid);
    /*
    * Function to set the areaXDerivative with respect to x and y for the given vertex
    */
    void setAkDerivative();
    /**
    * Function to get the AkXDerivative with respect to X for the given vertex
    * faceid : the id of the face on which the derivative is to be returned
    */
    double getAkXDerivative(unsigned int faceid);
    /**
    * Function to get the AkXDerivative with respect to Y for the given vertex
    * faceid : the id of the face on which the derivative is to be returned
    */
    double getAkYDerivative(unsigned int faceid);
    // %%%%%%%%%%%%%%%%%%%%%% X derivative %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
    // ***************** SET functions ******************** //
    /**
     * function to set mu deriavtive of this vertex with respect to X
     */
    void setMuXDerivative();
    /**
     * Function to set the mu squared derivative of this vertex with respect to X
     */
    void setMuSquaredXDerivative();
     /**
      * first Term X Derivative
      */
    void setFirstTermXDerivative();
    /**
     * set second term X derivative
     */
    void setSecondTermXDerivative();
    /**
     * set third term X derivative
     */
    void setThirdTermXDerivative();
    // ***************** GET functions ******************** //
    /**
     * function to get mu1 deriavtive of this vertex with the key face id
     */
    double getMu1XDerivative(unsigned int);
    /**
     * function to get mu4 deriavtive of this vertex with the key face id
     */
    double getMu4XDerivative(unsigned int); 
    /**
     * get frist term x derivaitve
     */
    double getFirstTermXDerivative();
     /**
      * get second term x derivative
      */
    double getSecondTermXDerivative();
    /**
     * get third term x derivative();
     */
    double getThirdTermXDerivative();
    /**
     * get mu1sqderivative 
     */
    double getMu1SquaredXDerivative(unsigned int facid);
     /**
     * get mu2sqderivative 
     */
    double getMu2SquaredXDerivative(unsigned int faceid);
     /**
     * get mu3sqderivative 
     */
    double getMu3SquaredXDerivative(unsigned int faceid);
     /**
     * get mu4sqderivative 
     */
    double getMu4SquaredXDerivative(unsigned int faceid);
    // %%%%%%%%%%%%%%%%%%%%%% Y derivative %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
    // ********************** SET function ********************* //
    /**
     * function to set mu deriavtive of this vertex with respect to X
     */
    void setMuYDerivative();
    /**
     * Function to set mu derivative of this vertex with respect to Y
     */
    void setMuSquaredYDerivative();
     /**
      * first Term Y Derivative
      */
    void setFirstTermYDerivative();
    /**
     * set second term Y derivative
     */
    void setSecondTermYDerivative();
    /**
     * set third term Y derivative
     */
     void setThirdTermYDerivative();
     // ******************* GET function ******************** //
    /**
     * get frist term y derivaitve
     */
    double getFirstTermYDerivative();
     /**
      * get second term y derivative
      */
    double getSecondTermYDerivative();
    /**
     * get third term y derivative
     */
    double getThirdTermYDerivative();
    /**
     * get mu1 y derivative
     */
    double getMu1YDerivative(unsigned int);
    /**
     * get mu2 y derivative();
     */
    double getMu4YDerivative(unsigned int);
     /**
     * get mu1sqderivative 
     */
    double getMu1SquaredYDerivative(unsigned int facid);
     /**
     * get mu2sqderivative 
     */
    double getMu2SquaredYDerivative(unsigned int faceid);
     /**
     * get mu3sqderivative 
     */
    double getMu3SquaredYDerivative(unsigned int faceid);
     /**
     * get mu4sqderivative 
     */
    double getMu4SquaredYDerivative(unsigned int faceid);
     
     // *************************************************************//
   /*
    * Set function to set the value of function1
    * faceid : id of the face to set the key
    * edge* : the pointer to the edge, of which origin and destination is used to calculate function1
    * 
    */
    void setFunction1();
    /*
    * Set function to set the value of function2
    * faceid : id of the face to set the key
    * edge* : the pointer to the edge, of which origin and destination is used to calculate function2
    * 
    */
    void setFunction2();
    /*
    * Set function to set the value of function3
    * faceid : id of the face to set the key
    * edge* : the pointer to the edge, of which origin and destination is used to calculate function3
    * 
    */
    void setFunction3();
   /*
    * Set function to set the value of Ak
    * will set the value of Ak for both Right and Left face of this current edge
    */
    void setAk();
    //*****************************************************************//
    /**
     * function to get the value of Ak
     * faceid: face id of the face on which the value of Ak is to be returned
     */
    double getAk(unsigned int faceid);
  /**
     * function to get the value of function1
     * faceid: face id of the face on which the value of function1 is to be returned
     */
    double getFunction1(unsigned int faceid);
  /**
     * function to get the value of function2
     * faceid: face id of the face on which the value of function2 is to be returned
     */
    double getFunction2(unsigned int faceid);
  /**
     * function to get the value of function3
     * faceid: face id of the face on which the value of function3 is to be returned
     */
    double getFunction3(unsigned int faceid);
///////////////////////////////////////////////////////////////////////
  /*
   * Return an arbitrary outgoing edge from this vertex.
   * <- an edge whose origin is this vertex;
   *    null if isolated
   */
  Edge *getEdge();

  /*
   * Add a given outgoing edge to this vertex.
   * edge -> an edge whose origin is this vertex;
   *         must be nonnull
   */
  void addEdge(Edge *edge);
    
  /*
   * Remove a given outgoing edge from this vertex.
   * edge -> an edge whose origin is no longer at this vertex;
   *         must be nonnull
   */
  void removeEdge(Edge *edge);
  /* -- protected instance methods ----------------------------------------- */

  protected:
  /*
   * Initialize this vertex at the origin with no outgoing edges.
   * cell -> the cell that this vertex belongs to;
   *         must be nonnull
   */
  Vertex(Cell *cell);

  /*
   * Release the storage occupied by the contents of this vertex.
   */
  ~Vertex();

  /* -- private instance variables ----------------------------------------- */

  private:

  /*
   * The cell that this vertex belongs to.
   * Nonnull.
   */
  Cell *cell;

  /*
   * The ID of this vertex.
   * Positive.
   */
  unsigned int id;

  /*
   * An arbitrary outgoing edge of this vertex.
   * Null if isolated.
   */
  Edge *edge;

};

/* -- inline instance methods ---------------------------------------------- */

inline Cell *Vertex::getCell()
{
  return cell;
}

inline unsigned int Vertex::getID()
{
  return id;
}

inline Edge *Vertex::getEdge()
{
  return edge;
}
inline double Vertex::getAk(unsigned int faceid)
{
  return this->Ak[faceid];
}
inline double Vertex::getFunction1(unsigned int faceid)
{
  return this->Function1[faceid];
}
inline double Vertex::getFunction2(unsigned int faceid)
{
  return this->Function2[faceid];
}
inline double Vertex::getFunction3(unsigned int faceid)
{
  return this->Function3[faceid];
}
inline double Vertex::getMu1SquaredXDerivative(unsigned int faceid)
{
  return this->mu1SquaredXDerivative[faceid];
}
inline double Vertex::getMu2SquaredXDerivative(unsigned int faceid)
{
  return this->mu2SquaredXDerivative[faceid];
}
inline double Vertex::getMu3SquaredXDerivative(unsigned int faceid)
{
  return this->mu3SquaredXDerivative[faceid];
}
inline double Vertex::getMu4SquaredXDerivative(unsigned int faceid)
{
  return this->mu4SquaredXDerivative[faceid];
}
inline double Vertex::getMu1XDerivative(unsigned int faceid)
{
  return this->mu1XDerivative[faceid];
}
inline double Vertex::getMu4XDerivative(unsigned int faceid)
{
  return this->mu4XDerivative[faceid];
}
inline double Vertex::getMu1YDerivative(unsigned int faceid)
{
  return this->mu1YDerivative[faceid];
}
inline double Vertex::getMu4YDerivative(unsigned int faceid)
{
  return this->mu4YDerivative[faceid];
}
inline double Vertex::getFirstTermXDerivative()
{
  return this->firstTermXDerivative;
}
inline double Vertex::getFirstTermYDerivative()
{
  return this->firstTermYDerivative;
}
inline double Vertex::getSecondTermXDerivative()
{
  return this->secondTermXDerivative;
}
inline double Vertex::getSecondTermYDerivative()
{
  return this->secondTermYDerivative;
}
inline double Vertex::getThirdTermXDerivative()
{
  return this->thirdTermXDerivative;
}
inline double Vertex::getThirdTermYDerivative()
{
  return this->thirdTermYDerivative;
}
/* ----------------------------------------------------------------------------
 * VertexEdgeIterator
 * ------------------------------------------------------------------------- */

/*
 * Enumerates the outgoing edges of a given vertex in counterclockwise order.
 */
class VertexEdgeIterator
{

  /* -- public instance methods -------------------------------------------- */

  public:

  /*
   * Initialize this edge iterator over a given vertex.
   * vertex -> the vertex to iterate the edges of;
   *           must be nonnull
   */
  VertexEdgeIterator(Vertex *vertex)
  {
    // pick an arbitrary edge in the vertex orbit

    start = vertex->getEdge();
    edge  = start;
  }

  /*
   * Release the storage occupied by this edge iterator.
   */
  ~VertexEdgeIterator()
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

    // get the next edge in the counterclockwise orbit of the vertex, but
    // return the current edge
    // reset to null if we've come back to the start

    Edge *next = current->Onext();

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

#endif /* #ifndef vertexINCLUDED */

