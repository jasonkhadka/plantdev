
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
  Vertex *getVertex();
  /**
   * Return the ith vertex of this face. 
   * <- ith vertex on this face
   */
   Vertex *getIthVertex(unsigned int position);
  /**
   * Add a given vertex to this face
   * vertex -> a vertex to be added to this face
   *           is one of the vertex forming the face
   *           must be nonnull
   */
   void addVertex(Vertex *vertex);

   /**
   * Remove a given vertex from this face
   * vertex -> a vertex to be removed to this face
   *           is one of the vertex forming the face
   *           must be nonnull
   */
   
   void removeVertex(Vertex *vertex);
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
  /*
   * Area of face
   * Positive
   * stored as Area = 1/2* \sum_{edges} A_k
   */
   double areaOfFace;
 /* -- public instance methods ----------------------------------------- */
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
  //***************added features*******************************************//
   /**
    * The vertices in this face
    * nonnull
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
  //***************end added features****************************************//
  /* ----------Public instance methods ---------*/
  public:

  /*
   * Return the number of vertices in this cell.
   * <- the number of vertices
   */
  unsigned int countVertices();

};
/* ----- inline instance methods --------------------------------*/
inline unsigned int Face::countVertices()
{
  return vertexCount;
}

/* -- inline instance methods ---------------------------------------------- */

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
inline Vertex *Face::getVertex()
{
  return vertices[0];
}
inline Vertex *Face::getIthVertex(unsigned int position)
{
  return vertices[position];
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

