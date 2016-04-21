
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
  double Xcoordinate = 0;
  double Ycoordinate = 0;
  double Zcoordinate = 0;
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

