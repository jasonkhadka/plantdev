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
	 	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
	 	// %%%%%%%%%%%%%%%%%% MAIN FUNCTIONS %%%%%%%%%%%%%%%%%%%% //
	 	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //

	 	//~~~~~~~~~~~ ~~~~~~~~~~~ ~~~~~~~~~~~ ~~~~~~~~~~~ ~~~//
	 	//~~~~~~~~~~~~~~~   x derivatives   ~~~~~~~~~~~~~~~ //
	 	//~~~~~~~~~~~ ~~~~~~~~~~~ ~~~~~~~~~~~ ~~~~~~~~~~~ ~~~//
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


	 	//~~~~~~~~~~~ ~~~~~~~~~~~ ~~~~~~~~~~~ ~~~~~~~~~~~ ~~~//
	 	//~~~~~~~~~~~~~~~   y derivatives   ~~~~~~~~~~~~~~~ //
	 	//~~~~~~~~~~~ ~~~~~~~~~~~ ~~~~~~~~~~~ ~~~~~~~~~~~ ~~~//
	 	
	 	


	 	// ******************************************************************* //


	 	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
		// %%%%%%%%%%%%%%%%%% Sub FUNCTIONS %%%%%%%%%%%%%%%%%%%% //
		// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //


	 	//~~~~~~~~~~~ ~~~~~~~~~~~ ~~~~~~~~~~~ ~~~~~~~~~~~ ~~~//
	 	//~~~~~~~~~~~~~~~   x derivatives   ~~~~~~~~~~~~~~~ //
	 	//~~~~~~~~~~~ ~~~~~~~~~~~ ~~~~~~~~~~~ ~~~~~~~~~~~ ~~~//


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

	     //~~~~~~~~~~~ ~~~~~~~~~~~ ~~~~~~~~~~~ ~~~~~~~~~~~ ~~~//
	 	//~~~~~~~~~~~~~~~   y derivatives   ~~~~~~~~~~~~~~~ //
	 	//~~~~~~~~~~~ ~~~~~~~~~~~ ~~~~~~~~~~~ ~~~~~~~~~~~ ~~~//
		 /**
	 	 * xtilde (x centralised coordiante) derivative with respect to 
	 	 * 			y_k coordiante. 
	 	 */
		 double xtildeYDerivative(Vertex *first, Vertex *second, Face *face);
		 /**
	 	 * ytilde (y centralised coordiante) derivative with respect to 
	 	 * 			y_k coordiante. 
	 	 */
		 double ytildeYDerivative(Vertex *first, Vertex *second, Face *face);

	 	/* y_i derivative with respect to y_k
		 * derivative of first with respect to second
		 * face : current face on which derivatives are calculated
		 */
		 double yiYDerivative(Vertex*, Vertex*, Face*);
		 /* ex1 derivative with respect to y_k
		 * derivative of first with respect to y coordinate of second
		 */
		 double ex1YDerivative(Vertex *first, Vertex *second, Face* face);

		 /**
		 * ex2 derivative with respect to y_k
		 * derivative of first with respect to y coordinate of second on Face face
		 */
		 double ex2YDerivative(Vertex* first, Vertex *second, Face * face);
		 /**
		  * ex3 derivative with respect to y_k
		  * derivative of first with respect to y coordinate of second on Face face
		  */
		 double ex3YDerivative(Vertex* first, Vertex * second, Face* face);
		 /* ey1 derivative with respect to y_k
		 * derivative of first with respect to y coordinate of second
		 */
		 double ey1YDerivative(Vertex *first, Vertex *second, Face* face);

		 /**
		 * ey2 derivative with respect to y_k
		 * derivative of first with respect to y coordinate of second on Face face
		 */
		 double ey2YDerivative(Vertex* first, Vertex *second, Face * face);
		 /**
		  * ey3 derivative with respect to y_k
		  * derivative of first with respect to y coordinate of second on Face face
		  */
		 double ey3YDerivative(Vertex* first, Vertex * second, Face* face);
		 /**
		  * derivative of x component of normal of the face with respect to y_k of second
		  */
		 double ncxYDerivative(Vertex*, Vertex*, Face*);
		 /**
		  * derivative of y component of normal of the face with respect to y_k of second
		  */
		 double ncyYDerivative(Vertex* , Vertex*, Face* );
		 /**
		  * derivative of z component of normal of the face with respect to y_k of second
		  */
		 double nczYDerivative(Vertex* , Vertex*, Face* );
		 /**
		 * CPi vector norm derivative with respect to y_k
		 * derivative of CPi vector norm of first with respect to y coordinate of second
		 * face : current face on which derivatvie is calculated
		 */
		 double cpinormYDerivative(Vertex *first, Vertex *second, Face* face);
		 /**
		  * function to calculate y derivative of x component of pi vector of face
		  */
		 double pixYDerivative(Vertex*, Vertex*, Face*);
		 /**
		  * function to calculate y derivative of y component of pi vector of face
		  */
		 double piyYDerivative(Vertex*, Vertex*, Face*);
		 /**
		  * function to calculate y derivative of z component of pi vector of face
		  */
		 double pizYDerivative(Vertex*, Vertex*, Face*);
		/**
		  * function to calculate y derivative of norm of un-normalised 
		  * vector normal of the face
		  */
		 double nctildenormYDerivative(Vertex*, Vertex*, Face*);
		 /**
		  * function to calculate y derivative of X component of un-normalised 
		  * vector normal of the face
		  */
		 double ncxtildeYDerivative(Vertex*, Vertex*, Face*);
		 /**
		  * function to calculate y derivative of Y component of un-normalised 
		  * vector normal of the face
		  */
		 double ncytildeYDerivative(Vertex*, Vertex*, Face*);
		 /**
		  * function to calculate y derivative of Z component of un-normalised 
		  * vector normal of the face
		  */
		 double ncztildeYDerivative(Vertex*, Vertex*, Face*);
		 /**
		  * derivative of areatotal with respect to y_k of second vertex
		  */
		  double areatotalYDerivative(Vertex*, Vertex*, Face*);
		  /**
		 * alpha_j derivative with respect to y_k
		 * derivative of first with respect to second
		 * face : current face on which derivatives are calculated
		 */
		double alphaYDerivative(Vertex*, Vertex*, Face*);
		/**
		 * beta_j derivative with respect to y_k
		 * derivative of first with respect to second
		 * face : current face on which derivatives are calculated
		 */
		double betaYDerivative(Vertex*, Vertex*, Face*);
		/**
		 * gamma_j derivative with respect to Y_k
		 * derivative of first with respect to second
		 * face : current face on which derivatives are calculated
		 */
		double gammaYDerivative(Vertex*, Vertex*, Face*);


	     //~~~~~~~~~~~ ~~~~~~~~~~~ ~~~~~~~~~~~ ~~~~~~~~~~~ ~~~//
	 	//~~~~~~~~~~~~~~~   z derivatives   ~~~~~~~~~~~~~~~ //
	 	//~~~~~~~~~~~ ~~~~~~~~~~~ ~~~~~~~~~~~ ~~~~~~~~~~~ ~~~//
		 /**
	 	 * xtilde (x centralised coordiante) derivative with respect to 
	 	 * 			z_k coordiante. 
	 	 */
		 double xtildeZDerivative(Vertex *first, Vertex *second, Face *face);
		 /**
	 	 * ytilde (y centralised coordiante) derivative with respect to 
	 	 * 			z_k coordiante. 
	 	 */
		 double ytildeZDerivative(Vertex *first, Vertex *second, Face *face);

	 	/* y_i derivative with respect to z_k
		 * derivative of first with respect to second
		 * face : current face on which derivatives are calculated
		 */
		 double ziZDerivative(Vertex*, Vertex*, Face*);
		 /* ex1 derivative with respect to z_k
		 * derivative of first with respect to z coordinate of second
		 */
		 double ex1ZDerivative(Vertex *first, Vertex *second, Face* face);

		 /**
		 * ex2 derivative with respect to z_k
		 * derivative of first with respect to z coordinate of second on Face face
		 */
		 double ex2ZDerivative(Vertex* first, Vertex *second, Face * face);
		 /**
		  * ex3 derivative with respect to z_k
		  * derivative of first with respect to z coordinate of second on Face face
		  */
		 double ex3ZDerivative(Vertex* first, Vertex * second, Face* face);
		 /* ey1 derivative with respect to z_k
		 * derivative of first with respect to z coordinate of second
		 */
		 double ey1ZDerivative(Vertex *first, Vertex *second, Face* face);

		 /**
		 * ey2 derivative with respect to z_k
		 * derivative of first with respect to z coordinate of second on Face face
		 */
		 double ey2ZDerivative(Vertex* first, Vertex *second, Face * face);
		 /**
		  * ey3 derivative with respect to z_k
		  * derivative of first with respect to z coordinate of second on Face face
		  */
		 double ey3ZDerivative(Vertex* first, Vertex * second, Face* face);
		 /**
		  * derivative of x component of normal of the face with respect to z_k of second
		  */
		 double ncxZDerivative(Vertex*, Vertex*, Face*);
		 /**
		  * derivative of y component of normal of the face with respect to z_k of second
		  */
		 double ncyZDerivative(Vertex* , Vertex*, Face* );
		 /**
		  * derivative of z component of normal of the face with respect to z_k of second
		  */
		 double nczZDerivative(Vertex* , Vertex*, Face* );
		 /**
		 * CPi vector norm derivative with respect to z_k
		 * derivative of CPi vector norm of first with respect to z coordinate of second
		 * face : current face on which derivatvie is calculated
		 */
		 double cpinormZDerivative(Vertex *first, Vertex *second, Face* face);
		 /**
		  * function to calculate z derivative of x component of pi vector of face
		  */
		 double pixZDerivative(Vertex*, Vertex*, Face*);
		 /**
		  * function to calculate Z derivative of y component of pi vector of face
		  */
		 double piyZDerivative(Vertex*, Vertex*, Face*);
		 /**
		  * function to calculate z derivative of z component of pi vector of face
		  */
		 double pizZDerivative(Vertex*, Vertex*, Face*);
		/**
		  * function to calculate Z derivative of norm of un-normalised 
		  * vector normal of the face
		  */
		 double nctildenormZDerivative(Vertex*, Vertex*, Face*);
		 /**
		  * function to calculate Z derivative of X component of un-normalised 
		  * vector normal of the face
		  */
		 double ncxtildeZDerivative(Vertex*, Vertex*, Face*);
		 /**
		  * function to calculate Z derivative of Y component of un-normalised 
		  * vector normal of the face
		  */
		 double ncytildeZDerivative(Vertex*, Vertex*, Face*);
		 /**
		  * function to calculate Z derivative of Z component of un-normalised 
		  * vector normal of the face
		  */
		 double ncztildeZDerivative(Vertex*, Vertex*, Face*);
		 /**
		  * derivative of areatotal with respect to Z_k of second vertex
		  */
		  double areatotalZDerivative(Vertex*, Vertex*, Face*);
		  /**
		 * alpha_j derivative with respect to Z_k
		 * derivative of first with respect to second
		 * face : current face on which derivatives are calculated
		 */
		double alphaZDerivative(Vertex*, Vertex*, Face*);
		/**
		 * beta_j derivative with respect to z_k
		 * derivative of first with respect to second
		 * face : current face on which derivatives are calculated
		 */
		double betaZDerivative(Vertex*, Vertex*, Face*);
		/**
		 * gamma_j derivative with respect to z_k
		 * derivative of first with respect to second
		 * face : current face on which derivatives are calculated
		 */
		double gammaZDerivative(Vertex*, Vertex*, Face*);




		 //&&*************************************************************************** //
		 //&********************  Numerical Derivatives   ****************************** //
		 //&*************************************************************************** //

		 // ******** Numerical Derivative of Xtilde X derivative  ********************* //
		 double numericalXtildeXDerivative(Vertex*, Vertex*, Face*, double);
		// ******** Numerical Derivative of Ytilde X derivative  ********************* //
		 double numericalYtildeXDerivative(Vertex*, Vertex*, Face*, double);
		 // ******* Numerical Derivative of Xtilde with Y_k derivative *************** //
		 double numericalXtildeYDerivative(Vertex* , Vertex*, Face*, double);
		// ******* Numerical Derivative of Xtilde with Y_k derivative *************** //
		 double numericalYtildeYDerivative(Vertex* , Vertex*, Face*, double);
		// ******* Numerical Derivative of Xtilde with Z_k derivative *************** //
		 double numericalXtildeZDerivative(Vertex* , Vertex*, Face*, double);
		// ******* Numerical Derivative of Xtilde with Z_k derivative *************** //
		 double numericalYtildeZDerivative(Vertex* , Vertex*, Face*, double);
		    
		
		// *************************************************************************** //
		 // ********************      Miscellaneous      ****************************** //
		 // *************************************************************************** //
		 /**
		 * delta function for checking the ids
		 * takes the unsigned ints or the ids of two objects
		 * and returns : 1.0 if same 
		 *               0.0 if different
		 */
		 double deltafunction(unsigned int, unsigned int);

		 // ******************************************************************************************* //
		 // ******** Function metioned but left to write  ********************* //
		 
		 
};
//********** Inline Instance methods ************************//
inline double CentralisedDerivative::deltafunction(unsigned int first, unsigned int second){
	return (first == second)*1.0;
}




#endif