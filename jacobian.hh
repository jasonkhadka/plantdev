#include "cell.hh"
#include "edge.hh"
#include "face.hh"
#include "vertex.hh"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>//for math operations sqrt and others

/**----------------------------------------------------------------------------------------------
 * Jacobian library : 
 * This calculates Energy and Jacobian of a given cell (tissue)
 * ----------------------------------------------------------------------------------------------
 */
 
void Jacobian(Cell *cell);