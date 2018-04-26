
#include "ParticleLocalityGrid.h"
#include "ParticlePosition.h"
#include "ParticleBoundingBox.h"
#include "pstools.h"

REGISTER_PARTICLESTUFF(ParticleLocalityGrid,"Attribute:ParticleLocalityGrid");

ParticleLocalityGrid::ParticleLocalityGrid(Particles *ps,const std::string& name)
		:ParticleLocality(ps,name)
{
	new PSParamInt(this,&new_grid_x,1,"xres","X Resolution",
		"Number of grid cells in X direction.");
	new PSParamInt(this,&new_grid_y,1,"yres","Y Resolution",
		"Number of grid cells in Y direction.");
	new PSParamInt(this,&new_grid_z,1,"zres","Z Resolution",
		"Number of grid cells in Z direction.");
	new Attached<ParticleBoundingBox>(this,&bbox);

	// initialize to one big cell
	grid_x = grid_y = grid_z = 1;
	grid = new int[1];
	grid[0] = -1;
	next = new int[1];
	next[0] = -1;

	next_size = 1;
}

ParticleLocalityGrid::~ParticleLocalityGrid()
{
	delete [] grid;
	delete [] next;
}


/**
 * Fills a list with neighbor candidates.
 * @param p				The particle in question.
 * @param neighbors		List of indices of neighbors.
 */
void ParticleLocalityGrid::getNeighbors(const unsigned int p, const double queryRadius, std::list<unsigned int> &neighbors)
{
	// Find the cell containing particle p
	findCell(p);

	//We need to know how many cells are touched by the 
	//readius, therefore we need the ceil value
	double reach_x = ceil(queryRadius/cell_dim[0]);
	double reach_y = ceil(queryRadius/cell_dim[1]);
	double reach_z = ceil(queryRadius/cell_dim[2]);

	// Determine from which cells to get neighbors
	int min_x = (int)gmMax(0, cell_x - reach_x);
	int max_x = (int)gmMin(grid_x, cell_x + reach_x + 1.0);
	int min_y = (int)gmMax(0, cell_y - reach_y);
	int max_y = (int)gmMin(grid_y, cell_y + reach_y + 1.0);
	int min_z = (int)gmMax(0, cell_z - reach_z);
	int max_z = (int)gmMin(grid_z, cell_z + reach_z + 1.0);

	gmVector3 px = position->getPosition(p);
	gmVector3 d;
	double qr2 = queryRadius * queryRadius;

//I broke up the loop to save several times the if tests for identity
//also particles lying slightly outside the radius are no longer excluded.
//We would gain a small amount of computation for each excluded particle 
//during the energy evaluation, but the probability of having particles 
//in a voxel close to the current position that are outside this radius is small
#if 0
	for ( int x = min_x; x < cell_x; x++) {
		for ( int y = min_y; y < max_y; y++) {
			for ( int z = min_z; z < max_z; z++) {
				for (int i = grid[at(x,y,z)]; i != -1; i = next[i]) {
//					if (i == p) continue;
//					d = position->getPosition(i) - px;
//					if (d.lengthSquared() <= qr2) {
						neighbors.push_back(i);
//					}
				}
			}
		}
	}
	for ( int x = cell_x+1; x < max_x; x++) {
		for ( int y = min_y; y < max_y; y++) {
			for ( int z = min_z; z < max_z; z++) {
				for (int i = grid[at(x,y,z)]; i != -1; i = next[i]) {
//					if (i == p) continue;
//					d = position->getPosition(i) - px;
//					if (d.lengthSquared() <= qr2) {
						neighbors.push_back(i);
//					}
				}
			}
		}
	}

	for ( int y = min_y; y < cell_y; y++) {
		for ( int z = min_z; z < max_z; z++) {
			for (int i = grid[at(cell_x,y,z)]; i != -1; i = next[i]) {
//				if (i == p) continue;
//				d = position->getPosition(i) - px;
//				if (d.lengthSquared() <= qr2) {
					neighbors.push_back(i);
//				}
			}
		}
	}

	for ( int y = cell_y+1; y < max_y; y++) {
		for ( int z = min_z; z < max_z; z++) {
			for (int i = grid[at(cell_x,y,z)]; i != -1; i = next[i]) {
//				if (i == p) continue;
//				d = position->getPosition(i) - px;
//				if (d.lengthSquared() <= qr2) {
					neighbors.push_back(i);
//				}
			}
		}
	}

	for ( int z = min_z; z <cell_z; z++) {
		for (int i = grid[at(cell_x,cell_y,z)]; i != -1; i = next[i]) {
//			if (i == p) continue;
//			d = position->getPosition(i) - px;
//			if (d.lengthSquared() <= qr2) {
				neighbors.push_back(i);
//			}
		}
	}

	for ( int z = cell_z+1; z <max_z; z++ ) {
		for (int i = grid[at(cell_x,cell_y,z)]; i != -1; i = next[i]) {
//			if (i == p) continue;
//			d = position->getPosition(i) - px;
//			if (d.lengthSquared() <= qr2) {
				neighbors.push_back(i);
//			}
		}
	}

	for (int i = grid[at(cell_x,cell_y,cell_z)]; i != -1; i = next[i]) {
		//for the center cell we need this test!
		if (i == p) continue;
//		d = position->getPosition(i) - px;
//		if (d.lengthSquared() <= qr2) {
			neighbors.push_back(i);
//		}
	}
#else
	
	for ( int x = min_x; x < max_x; x++) {
		for ( int y = min_y; y < max_y; y++) {
			for ( int z = min_z; z < max_z; z++) {
				for (int i = grid[at(x,y,z)]; i != -1; i = next[i]) {
					if (i == p) continue;
					d = position->getPosition(i) - px;
					if (d.lengthSquared() <= qr2) {
						neighbors.push_back(i);
					}
				}
			}
		}
	}
#endif
}

/** This update procedure fills the grid with a linear pass through the particles.
 * It first clears the grid with -1 entries. Then for each particle, enters the index of
 * the particle into its corresponding grid cell and the contents of that grid cell into
 * the particle's corresponding gridlist entry. Thus the grid contains the index of a
 * particle in it, and that particle's gridlist entry points to another particle in it.
 * A -1 entry terminates the list.
 */

void ParticleLocalityGrid::update()
{
	if (new_grid_x != grid_x || new_grid_y != grid_y || new_grid_z != grid_z) {
		delete [] grid;
		grid_x = new_grid_x; grid_y = new_grid_y; grid_z = new_grid_z;
		grid = new int[grid_x*grid_y*grid_z];
		assert(grid);
	}
	for ( int i = 0; i < grid_x*grid_y*grid_z; i++)
		grid[i] = -1;

	if ((int) ps->size() > next_size) {
		delete [] next;
		next_size = (int) 1.5*ps->size();
		next = new int[next_size];
		assert(next);
	}

	for (unsigned int i = 0; i < ps->size(); i++)
		next[i] = -1;	

	// calculate cell dimensions may change when bounding box is changing
	gmVector3 gs(bbox->max - bbox->min);
	cell_dim = gmVector3(gs[0]/grid_x, gs[1]/grid_y, gs[2]/grid_z);

	for(unsigned int i = 0; i < ps->size(); i++) {
		unsigned int xyz = findCell(i);
		next[i] = grid[xyz];
		grid[xyz] = i;
	}
}

/**
 * Finds the appropriate cell for a particle. The result is stored
 * in private member variables cell_x, cell_y, cell_z.
 */
unsigned int ParticleLocalityGrid::findCell(unsigned int i)
{
	const gmVector3 &grid_origin=bbox->min;
	cell_x = (int)((position->getPosition(i)[0] - grid_origin[0]) / cell_dim[0]);
	cell_y = (int)((position->getPosition(i)[1] - grid_origin[1]) / cell_dim[1]);
	cell_z = (int)((position->getPosition(i)[2] - grid_origin[2]) / cell_dim[2]);

	clampInt(cell_x, grid_x); // [0, grid_x)
	clampInt(cell_y, grid_y); // [0, grid_y)
	clampInt(cell_z, grid_z); // [0, grid_z)

	return at(cell_x,cell_y,cell_z);
}

/**
 * Finds the vertices of the appropriate cell for a particle.
 */
void ParticleLocalityGrid::cellBox(unsigned int i, gmVector3 &a, gmVector3 &b)
{
	findCell(i);
	a = bbox->min + gmVector3(cell_x * cell_dim[0], cell_y * cell_dim[1], cell_z*cell_dim[2]);
	b = a + cell_dim;
}

