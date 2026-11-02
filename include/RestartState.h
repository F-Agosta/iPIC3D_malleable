#ifndef RestartState_H
#define RestartState_H

#include <string>
#include <sstream>
#include <iomanip>
#include "VCtopology3D.h"
#include "Grid3DCU.h"
#include "aligned_vector.h"
#include "arrays.h"

using namespace std;

/*! RestartState encapsulates all state and logic related to restarting a
 *  simulation from a previously saved checkpoint. It is responsible for:
 *    1. Reading and storing the topology of the previous simulation.
 *    2. Validating that the new topology is compatible with the old one.
 *    3. Providing helper methods to read fields and particles from restart files.
 */
class RestartState
{
public:
    struct RestartCoordinates {
        int prev[3];
        int internal[3];
    };

private:
    /*! previous subdomain dimensions */
    int XLEN_prev;
    int YLEN_prev;
    int ZLEN_prev;

    /*! multiplication factors: new_dim / prev_dim */
    int x_f;
    int y_f;
    int z_f;

    /*! remainder along each direction (should be 0 for a valid restart) */
    int x_m;
    int y_m;
    int z_m;

    /*! status of the restart: 0 = fresh start, 1 = restart */
    int restart_state;

    /*! directory where the restart files live */
    string RestartDirName;

    /*! last cycle*/
    mutable int last_cycle;

public:
    /*! Construct a RestartState for the given restart directory.
     *  After construction, call read_old_topology() and validate_new_topology()
     *  before using the read methods. */
    RestartState(const string& restartDirName);

    /*! Read the topology (XLEN, YLEN, ZLEN) of the previous simulation
     *  from settings.hdf in RestartDirName.
     *  @return 0 on success, -1 on error. */
    int read_old_topology();

    /*! Check that the new topology (passed in) is compatible with the old one,
     *  i.e. that each new dimension is an exact multiple of the old dimension.
     *  Also computes the x_f, y_f, z_f multiplication factors.
     *  @param XLEN_new  number of MPI subdomains in X for the new run
     *  @param YLEN_new  number of MPI subdomains in Y for the new run
     *  @param ZLEN_new  number of MPI subdomains in Z for the new run
     *  @return 1 if valid, -1 if invalid. */
    int validate_new_topology(int XLEN_new, int YLEN_new, int ZLEN_new);

    /*! Calculate the per-process coordinates used to derive file names and
     *  hyperslab offsets.
     *  @param vct  the virtual Cartesian topology of the current run */
    RestartCoordinates setup_coords(const VCtopology3D* vct) const;

    /*! Read the last cycle number from the restart file for this process.
     *  @return the last cycle number, or -1 on error. */
    int read_last_cycle() const;

    /*! Generate the name of the restart HDF5 file that this process should read.
     *  @return the full path to the restart file. */
    string generate_file_name(const RestartCoordinates& coords) const;

    /*! Read electromagnetic fields (B, E) and charge densities (rho) from
     *  the restart file.  Interface matches the old Collective method. */
    void read_field_restart(
        const VCtopology3D* vct,
        const Grid* grid,
        arr3_double Bxn, arr3_double Byn, arr3_double Bzn,
        arr3_double Ex,  arr3_double Ey,  arr3_double Ez,
        array4_double* rhons_, int ns) const;

    /*! Read particle data for one species from the restart file.
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
    int get_XLEN_prev() const { return XLEN_prev; }
    int get_YLEN_prev() const { return YLEN_prev; }
    int get_ZLEN_prev() const { return ZLEN_prev; }

    int get_x_f() const { return x_f; }
    int get_y_f() const { return y_f; }
    int get_z_f() const { return z_f; }

    int get_x_m() const { return x_m; }
    int get_y_m() const { return y_m; }
    int get_z_m() const { return z_m; }

    int get_restart_state() const { return restart_state; }

    string getRestartDirName() const { return RestartDirName; }
};

#endif // RestartState_H
