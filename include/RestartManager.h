#ifndef RestartManager_H
#define RestartManager_H

#include <string>
#include <sstream>
#include <iomanip>
#include <vector>
#include "VCtopology3D.h"
#include "Grid3DCU.h"
#include "aligned_vector.h"
#include "arrays.h"

using namespace std;

/*! RestartManager encapsulates all state and logic related to restarting a
 *  simulation from a previously saved checkpoint. It is responsible for:
 *    1. Storing the state of the previous simulation (RestartState).
 *    2. Computing, for each process, the intersection between its new
 *       subdomain and the old subdomains (file_slab / memory_slab pairs).
 *    3. Reading fields and particles from the restart files accordingly.
 *
 *  The old and the new topology are treated symmetrically: both are just
 *  partitions of the same global grid, so growing, shrinking and same-size
 *  restarts all follow from the same intersection logic.
 */
class RestartManager
{

public:

     struct RestartState {
        int nc[3]; // number of cells in each direction
        int pastLEN[3]; // number of subdomains in each direction before restart
        int currLEN[3]; // number of subdomains in each direction after restart
        int last_cycle;
        double L[3]; // physical box lengths, used for particle domain checks
    }

    /*! A 3D slab: offset and extent (in cells) along each direction */
    struct Slab {
        int offset[3];
        int count[3];
    };

    /*! Which restart file to read from, and which part of it */
    struct FileSlab {
        int file_coord[3]; // coordinates of the old subdomain -> restart_x*_y*_z*.hdf
        Slab slab;         // hyperslab inside that file (cells; fields use count+1 nodes)
    };

    /*! Where to put the data in this process's local arrays */
    struct MemorySlab {
        Slab slab;         // offset/count in the local grid (cells)
    };

    /*! 1:1 pairing of a file slab with the memory slab it maps to */
    struct RestartSlab {
        FileSlab   file;
        MemorySlab mem;
    }

private:

    /*! status of the restart: 0 = fresh start, 1 = restart */
    int restart_state;

    /*! directory where the restart files live */
    string RestartDirName;

    /*! full restart state, passed in at construction time */
    RestartState state;

    /*! Build the path of the restart file of the old subdomain coord */
    string file_name(const int coord[3]) const;

    /*! Check whether a particle at (px,py,pz) belongs to this process's
     *  new subdomain, based on its coordinates */
    bool owns_particle(double px, double py, double pz, const VCtopology3D* vct) const;

public:
    /*! Construct a RestartManager for the given restart directory and state.
     *  The RestartState must be fully populated (old/new topology, grid size,
     *  box lengths and last cycle) before construction. */
    RestartManager(const string& restartDirName, const RestartState& state);

    /*! Compute fresh, for this process, the intersection between its new
     *  subdomain (from vct) and the old subdomains: one RestartSlab per
     *  overlapping old subdomain. Slab counts are in cells. */
    vector<RestartSlab> compute_slabs(const VCtopology3D* vct) const;

    /*! Read the last cycle number from the rank-0 restart file.
     *  @return the last cycle number, or -1 on error. */
    int read_last_cycle() const;

    /*! Read electromagnetic fields (B, E) and charge densities (rho) from
     *  the restart file.  Interface matches the old Collective method. */
    void read_field_restart(
        const VCtopology3D* vct,
        const Grid* grid,
        arr3_double Bxn, arr3_double Byn, arr3_double Bzn,
        arr3_double Ex,  arr3_double Ey,  arr3_double Ez,
        array4_double* rhons_, int ns) const;

    /*! Read particle data for one species from the restart file.
     *  Only particles that belong to this process's subdomain are kept.
     *  Interface matches the old Collective method. */
    void read_particles_restart(
        const VCtopology3D* vct,
        int species_number,
        vector_double& u, vector_double& v, vector_double& w,
        vector_double& q, vector_double& x, vector_double& y,
        vector_double& z, vector_double& t) const;

    // ----------------------------------------------------------------
    // getters
    // ----------------------------------------------------------------
    int get_XLEN_prev() const { return state.pastLEN[0]; }
    int get_YLEN_prev() const { return state.pastLEN[1]; }
    int get_ZLEN_prev() const { return state.pastLEN[2]; }

    int get_restart_state() const { return restart_state; }

    string getRestartDirName() const { return RestartDirName; }

    const RestartState& getState() const { return state; }
};

#endif // RestartManager_H
