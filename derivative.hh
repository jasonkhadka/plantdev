#ifndef derivativeINCLUDED
#define derivativeINCLUDED

//forward declaring the classes
class Cell;
class Edge;
class Vector;
class Face;

class CentralisedDerivative
{


	 public:
	 	/**
	 	 * xtilde (x centralised coordiante) derivative with respect to 
	 	 * 			x_k coordiante. 
	 	 */
	 	double xtildeXDerivative(Vertex *first, Vertex *second, Face *face);
	 	/**
	 	 * ytilde (y centralised coordiante) derivative with respect to 
	 	 * 			x_k coordiante. 
	 	 */
	 	double ytildeXDerivative(Vertex *first, Vertex *second, Face *face);

		/* x_i derivative with respect to x_k
		 * derivative of first with respect to second
		 * face : current face on which derivatives are calculated
		 */
		double xiXDerivative(Vertex *first, Vertex *second, Face *face);
		/* ex1 derivative with respect to x_k
		 * derivative of first with respect to x coordinate of second
		 */
		double ex1XDerivative(Vertex *first, Vertex *second, Face* face);
		/**
		 * ex2 derivative with respect to x_k
		 * derivative of first with respect to x coordinate of second on Face face
		 */
		 double ex2XDerivative(Vertex* first, Vertex *second, Face * face);
		 /**
		  * ex3 derivative with respect to x_k
		  * derivative of first with respect to x coordinate of second on Face face
		  */
		  double ex3XDerivative(Vertex* first, Vertex * second, Face* face);
		 /** unit_y with respect to x_k */
		 /* ey1 derivative with respect to x_k
		 * derivative of first with respect to x coordinate of second
		 */
		double ey1XDerivative(Vertex *first, Vertex *second, Face* face);
		/**
		 * ey2 derivative with respect to x_k
		 * derivative of first with respect to x coordinate of second on Face face
		 */
		 double ey2XDerivative(Vertex* first, Vertex *second, Face * face);
		 /**
		  * ey3 derivative with respect to x_k
		  * derivative of first with respect to x coordinate of second on Face face
		  */
		  double ey3XDerivative(Vertex* first, Vertex * second, Face* face);

	// -----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%//

	// here the functions that are not directly called but are used to build bigger 
	// parts of derivatives should be put.
		/**
		 * delta function for checking the ids
		 * takes the unsigned ints or the ids of two objects
		 * and returns : 1.0 if same 
		 *               0.0 if different
		 */
		 double deltafunction(unsigned int, unsigned int);
		/**
		 * alpha_j derivative with respect to x_k
		 * derivative of first with respect to second
		 * face : current face on which derivatives are calculated
		 */
		double alphaXDerivative(Vertex*, Vertex*, Face*);
		/**
		 * beta_j derivative with respect to x_k
		 * derivative of first with respect to second
		 * face : current face on which derivatives are calculated
		 */
		double betaXDerivative(Vertex*, Vertex*, Face*);
		/**
		 * gamma_j derivative with respect to x_k
		 * derivative of first with respect to second
		 * face : current face on which derivatives are calculated
		 */
		double gammaXDerivative(Vertex*, Vertex*, Face*);
		/**
		 * CPi vector norm derivative with respect to x_k
		 * derivative of CPi vector norm of first with respect to x coordinate of second
		 * face : current face on which derivatvie is calculated
		 */
		 double cpinormXDerivative(Vertex *first, Vertex *second, Face* face);
		 /**
		  * function to calculate x derivative of x component of pi vector of face
		  */
		 double pixXDerivative(Vertex*, Vertex*, Face*);
		 /**
		  * function to calculate x derivative of y component of pi vector of face
		  */
		 double piyXDerivative(Vertex*, Vertex*, Face*);
		 /**
		  * function to calculate x derivative of z component of pi vector of face
		  */
		 double pizXDerivative(Vertex*, Vertex*, Face*);
		 /**
		  * derivative of x component of normal of the face with respect to x_k of second
		  */
		 double ncxXDerivative(Vertex*, Vertex*, Face*);
		 /**
		  * derivative of y component of normal of the face with respect to x_k of second
		  */
		 double ncyXDerivative(Vertex* , Vertex*, Face* );
		 /**
		  * derivative of z component of normal of the face with respect to x_k of second
		  */
		 double nczXDerivative(Vertex* , Vertex*, Face* );
		 
		  /**
		  * derivative of areatotal with respect to x_k of second vertex
		  */
		  double areatotalXDerivative(Vertex*, Vertex*, Face*);
		  /**
		  * function to calculate x derivative of norm of un-normalised 
		  * vector normal of the face
		  */
		 double nctildenormXDerivative(Vertex*, Vertex*, Face*);
		 /**
		  * function to calculate x derivative of X component of un-normalised 
		  * vector normal of the face
		  */
		 double ncxtildeXDerivative(Vertex*, Vertex*, Face*);
		 /**
		  * function to calculate x derivative of Y component of un-normalised 
		  * vector normal of the face
		  */
		 double ncytildeXDerivative(Vertex*, Vertex*, Face*);
		 /**
		  * function to calculate x derivative of Z component of un-normalised 
		  * vector normal of the face
		  */
		 double ncztildeXDerivative(Vertex*, Vertex*, Face*);
		 // ******************************************************************************************* //
		 // ******** Numerical Derivative of Xtilde X derivative  ********************* //
		 double numericalXtildeXDerivative(Vertex*, Vertex*, Face*, double);

		 // ******************************************************************************************* //
		 // ******** Function metioned but left to write  ********************* //

};
//********** Inline Instance methods ************************//
inline double CentralisedDerivative::deltafunction(unsigned int first, unsigned int second){
	return (first == second)*1.0;
}




#endif