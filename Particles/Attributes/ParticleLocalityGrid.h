//Declaration of ParticleLocalityGrid and ParticleLocalityCell.
//@file ParticleLocalityGrid.h
//@author Ed Bachta, Wen Su

#ifndef __PARTICLELOCALITYGRID_H__
#define __PARTICLELOCALITYGRID_H__

#include "ParticleLocality.h"
#include <vector>
#include <algorithm>

class ParticlePosition;
class ParticleBoundingBox;
/**
 * A ParticleLocalityGrid contains a 3D grid of ParticleLocalityCells.
 * This structure is useful for finding the neighbors of a particle
 * without interrogating all particles in the system.
 */
class ParticleLocalityGrid : public ParticleLocality
{
public:
	MAKE_PARTICLESTUFF_NAME();

	int grid_x;		///< X dimension of grid.
	int grid_y;		///< Y dimension of grid.
	int grid_z;		///< Z dimension of grid.

	int new_grid_x;	///< Requested X dimension
	int new_grid_y;	///< Requested Y dimension
	int new_grid_z;	///< Requested Z dimension

	int cell_x;		///< Holds cell X index after call to findCell.
	int cell_y;		///< Holds cell Y index after call to findCell.
	int cell_z;		///< Holds cell Z index after call to findCell.

	gmVector3 cell_dim;    ///< world space dimensions of one cell.

	int *grid; ///< The grid of cells.

	int *next;

	int next_size;

	/** bbox points to the ParticleBoundingBox used to indicate the spatial
	 * location of the grid.
	 */
	ParticleBoundingBox *bbox;

	/** Computes the cell corresponding to a particle's position.
	 */
	unsigned int findCell(unsigned int);

	int at(int i, int j, int k) { return k*grid_x*grid_y + j*grid_x + i; }

	/// Constructor.
	ParticleLocalityGrid(Particles *ps=NULL,const std::string& name=std::string("ParticleLocalityGrid"));

	~ParticleLocalityGrid();

	/** Updates the grid data structure to reflect the new particle positions.
	 */
	virtual void update();

	/// Fills a list with neighbors within queryRadius of particle i.
	virtual void getNeighbors(const unsigned int i, const double queryRadius, std::list<unsigned int> &neighbors);

	void cellBox(unsigned int i, gmVector3 &a, gmVector3 &b);
};

#endif
