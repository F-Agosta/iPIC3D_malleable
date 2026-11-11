/* iPIC3D was originally developed by Stefano Markidis and Giovanni Lapenta.
 * This release was contributed by Alec Johnson and Ivy Bo Peng.
 * Publications that use results from iPIC3D need to properly cite
 * 'S. Markidis, G. Lapenta, and Rizwan-uddin. "Multi-scale simulations of
 * plasma with iPIC3D." Mathematics and Computers in Simulation 80.7 (2010): 1509-1519.'
 *
 *        Copyright 2015 KTH Royal Institute of Technology
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "RestartState.h"
#include "ipichdf5.h"
#include "errors.h"
#include "ipicdefs.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>

using namespace std;

// ============================================================
// Constructor
// ============================================================

RestartState::RestartState(const string& restartDirName)
    : XLEN_prev(0), YLEN_prev(0), ZLEN_prev(0),
      x_f(1), y_f(1), z_f(1),
      x_m(0), y_m(0), z_m(0),
      restart_state(0),
      RestartDirName(restartDirName),
      last_cycle(0)
{
}

// ============================================================
// read_old_topology
// ============================================================

/*! Read the topology of the previous simulation from settings.hdf.
 *  Fills XLEN_prev, YLEN_prev, ZLEN_prev. */
int RestartState::read_old_topology()
{
#ifdef NO_HDF5
    eprintf("restart requires compiling with HDF5");
    return -1;
#else
    hid_t file_id;
    hid_t dataset_id;
    herr_t status;

    file_id = H5Fopen((RestartDirName + "/settings.hdf").c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0) {
        cerr << "Unable to open restart settings file: "
             << RestartDirName << "/settings.hdf" << endl;
        return -1;
    }

    // read old XLEN
    dataset_id = H5Dopen2(file_id, "/topology/XLEN", H5P_DEFAULT);
    status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &XLEN_prev);
    status = H5Dclose(dataset_id);

    // read old YLEN
    dataset_id = H5Dopen2(file_id, "/topology/YLEN", H5P_DEFAULT);
    status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &YLEN_prev);
    status = H5Dclose(dataset_id);

    // read old ZLEN
    dataset_id = H5Dopen2(file_id, "/topology/ZLEN", H5P_DEFAULT);
    status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &ZLEN_prev);
    status = H5Dclose(dataset_id);

    status = H5Fclose(file_id);

    restart_state = 1;
    return 0;
#endif
}

// ============================================================
// validate_new_topology
// ============================================================

/*! Check that the new topology is a multiple of the old one.
 *  Also computes x_f, y_f, z_f and x_m, y_m, z_m. */
int RestartState::validate_new_topology(int XLEN_new, int YLEN_new, int ZLEN_new)
{
    x_m = XLEN_new % XLEN_prev;
    y_m = YLEN_new % YLEN_prev;
    z_m = ZLEN_new % ZLEN_prev;

    if (x_m || y_m || z_m) {
        return -1;
    }

    x_f = XLEN_new / XLEN_prev;
    y_f = YLEN_new / YLEN_prev;
    z_f = ZLEN_new / ZLEN_prev;

    return 1;
}

// ============================================================
// setup_coords
// ============================================================

/*! Compute per-process file coordinates and internal offsets from the VCtopology. */
RestartState::RestartCoordinates
RestartState::setup_coords(const VCtopology3D* vct) const
{
    RestartCoordinates coords;
    coords.prev[0]     = vct->getCoordinates(0) / x_f;
    coords.prev[1]     = vct->getCoordinates(1) / y_f;
    coords.prev[2]     = vct->getCoordinates(2) / z_f;

    coords.internal[0] = vct->getCoordinates(0) % x_f;
    coords.internal[1] = vct->getCoordinates(1) % y_f;
    coords.internal[2] = vct->getCoordinates(2) % z_f;
    return coords;
}

// ============================================================
// generate_file_name
// ============================================================

/*! Return the path of the HDF5 restart file for this process. */
string RestartState::generate_file_name(
    const RestartCoordinates& coords) const
{
    stringstream ss;
    ss << RestartDirName << "/restart_x"
       << setw(3) << setfill('0') << coords.prev[0] << "_y"
       << setw(3) << setfill('0') << coords.prev[1] << "_z"
       << setw(3) << setfill('0') << coords.prev[2] << ".hdf";
    return ss.str();
}

// ============================================================
// read_last_cycle
// ============================================================

/*! Read the last cycle number stored in this process's restart file. */
int RestartState::read_last_cycle() const
{
#ifdef NO_HDF5
    eprintf("restart requires compiling with HDF5");
    return -1;
#else
    const RestartCoordinates coords = {};
    const string name_file = generate_file_name(coords);
    hid_t file_id = H5Fopen(name_file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0) {
        cout << "couldn't open restart file for last_cycle: " << name_file << endl;
        return -1;
    }

    last_cycle = 0;
    hid_t dataset_id = H5Dopen2(file_id, "/last_cycle", H5P_DEFAULT);
    if (dataset_id >= 0) {
        H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &last_cycle);
        H5Dclose(dataset_id);
    }
    H5Fclose(file_id);
    return last_cycle;
#endif
}

void RestartState::determineResizeMode(int old_proc, int new_proc) {
    if (old_proc > new_proc)
        resize_mode = SHRINKING;
    else if (old_proc > new_proc)
        resize_mode = GROWING;
    else
        resize_mode = SAME;
}

// ============================================================
// read_field_restart
// ============================================================

/*! Read electromagnetic fields and charge densities from the restart file.
 *  Corresponds to the old Collective::read_field_restart(). */
void RestartState::read_field_restart(
    const VCtopology3D* vct,
    const Grid* grid,
    arr3_double Bxn, arr3_double Byn, arr3_double Bzn,
    arr3_double Ex,  arr3_double Ey,  arr3_double Ez,
    array4_double* rhons_, int ns) const
{
#ifdef NO_HDF5
    eprintf("Require HDF5 to read from restart file.");
#else
    const int nxn = grid->getNXN();
    const int nyn = grid->getNYN();
    const int nzn = grid->getNZN();

    if (vct->getCartesian_rank() == 0) {
        printf("LOADING EM FIELD FROM RESTART FILE in %s/restart.hdf\n", RestartDirName.c_str());
    }

    const RestartCoordinates coords = setup_coords(vct);
    // 1. File & Internal Coordinates Setup
    // int file_coord[3] = {
    //     vct->getCoordinates(0) / x_f,
    //     vct->getCoordinates(1) / y_f,
    //     vct->getCoordinates(2) / z_f
    // };
    // int internal_coord_local[3] = {
    //     vct->getCoordinates(0) % x_f,
    //     vct->getCoordinates(1) % y_f,
    //     vct->getCoordinates(2) % z_f
    // };

    string name_file = generate_file_name(coords);
    // stringstream restart_name;
    // restart_name << RestartDirName << "/restart_x"
    //              << setw(3) << setfill('0') << file_coord[0] << "_y"
    //              << setw(3) << setfill('0') << file_coord[1] << "_z"
    //              << setw(3) << setfill('0') << file_coord[2] << ".hdf";
    // string name_file = restart_name.str();

    // 2. Open File and Get Last Cycle
    hid_t file_id = H5Fopen(name_file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0) {
        eprintf("Failed to open file: %s\n", name_file.c_str());
    }


    //NB: this piece is useless because last_cycle has already been read
    read_last_cycle();
    // hid_t dataset_id = H5Dopen2(file_id, "/last_cycle", H5P_DEFAULT);
    // H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &last_cycle);
    // H5Dclose(dataset_id);

    // 3. Grid Math & Shared Dataspace Initialization (Using Bx)
    stringstream ss;
    ss << "/fields/Bx/cycle_" << last_cycle;
    hid_t dataset_id = H5Dopen2(file_id, ss.str().c_str(), H5P_DEFAULT);
    hid_t shared_dataspace = H5Dget_space(dataset_id);

    hsize_t dims_out[3];
    H5Sget_simple_extent_dims(shared_dataspace, dims_out, NULL);

    hsize_t new_cells[3] = {
        (dims_out[0] - 1) / x_f,
        (dims_out[1] - 1) / y_f,
        (dims_out[2] - 1) / z_f
    };

    hsize_t curr_dims[3] = { new_cells[0] + 1, new_cells[1] + 1, new_cells[2] + 1 };
    hsize_t offset[3] = {
        (hsize_t)coords.internal[0] * new_cells[0],
        (hsize_t)coords.internal[1] * new_cells[1],
        (hsize_t)coords.internal[2] * new_cells[2]
    };

    // Allocate memory and setup hyperslab for ALL fields
    std::vector<double> temp_storage(curr_dims[0] * curr_dims[1] * curr_dims[2]);
    H5Sselect_hyperslab(shared_dataspace, H5S_SELECT_SET, offset, NULL, curr_dims, NULL);
    hid_t memspace_id = H5Screate_simple(3, curr_dims, NULL);

    // 4. Lambda helper to eliminate repetitive read/map code
    auto read_and_map = [&](const string& path, arr3_double arr) {
        hid_t d_id = H5Dopen2(file_id, path.c_str(), H5P_DEFAULT);
        H5Dread(d_id, H5T_NATIVE_DOUBLE, memspace_id, shared_dataspace, H5P_DEFAULT, temp_storage.data());
        int k = 0;
        for (int i = 1; i < nxn - 1; i++) {
            for (int j = 1; j < nyn - 1; j++) {
                for (int jj = 1; jj < nzn - 1; jj++) {
                    arr[i][j][jj] = temp_storage[k++];
                }
            }
        }
        H5Dclose(d_id);
    };

    // 5. Execute Reads
    // Bx was already opened to get dims, so read it and close it first
    H5Dread(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, shared_dataspace, H5P_DEFAULT, temp_storage.data());
    int k = 0;
    for (int i = 1; i < nxn - 1; i++)
        for (int j = 1; j < nyn - 1; j++)
            for (int jj = 1; jj < nzn - 1; jj++)
                Bxn[i][j][jj] = temp_storage[k++];
    H5Dclose(dataset_id);

    // Read remaining fields cleanly using the lambda
    ss.str(""); ss << "/fields/By/cycle_" << last_cycle; read_and_map(ss.str(), Byn);
    ss.str(""); ss << "/fields/Bz/cycle_" << last_cycle; read_and_map(ss.str(), Bzn);
    ss.str(""); ss << "/fields/Ex/cycle_" << last_cycle; read_and_map(ss.str(), Ex);
    ss.str(""); ss << "/fields/Ey/cycle_" << last_cycle; read_and_map(ss.str(), Ey);
    ss.str(""); ss << "/fields/Ez/cycle_" << last_cycle; read_and_map(ss.str(), Ez);

    // Read charge density (rho) for each species
    for (int is = 0; is < ns; is++) {
        ss.str(""); ss << "/moments/species_" << is << "/rho/cycle_" << last_cycle;
        dataset_id = H5Dopen2(file_id, ss.str().c_str(), H5P_DEFAULT);
        H5Sselect_hyperslab(shared_dataspace, H5S_SELECT_SET, offset, NULL, curr_dims, NULL);
        H5Dread(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, shared_dataspace, H5P_DEFAULT, temp_storage.data());
        array4_double& rhons = *rhons_;
        k = 0;
        for (int i = 1; i < nxn - 1; i++) {
            for (int j = 1; j < nyn - 1; j++) {
                for (int jj = 1; jj < nzn - 1; jj++) {
                    rhons[is][i][j][jj] = temp_storage[k++];
                }
            }
        }
        H5Dclose(dataset_id);
    }

    // 6. Cleanup
    H5Sclose(shared_dataspace);
    H5Sclose(memspace_id);
    H5Fclose(file_id);
#endif
}

// ============================================================
// read_particles_restart
// ============================================================

/*! Read particle data for one species from the restart file.
 *  Corresponds to the old Collective::read_particles_restart(). */
void RestartState::read_particles_restart(
    const VCtopology3D* vct,
    int species_number,
    vector_double& u, vector_double& v, vector_double& w,
    vector_double& q, vector_double& x, vector_double& y,
    vector_double& z, vector_double& t) const
{
#ifdef NO_HDF5
    eprintf("Require HDF5 to read from restart file.");
#else
    if (vct->getCartesian_rank() == 0 && species_number == 0) {
        printf("LOADING PARTICLES FROM RESTART FILE in %s/restart.hdf\n", RestartDirName.c_str());
    }

    // 1. File Coordinates & Internal Chunk ID
    const RestartCoordinates coords = setup_coords(vct);
    // int file_coord[3] = {
    //     vct->getCoordinates(0) / x_f,
    //     vct->getCoordinates(1) / y_f,
    //     vct->getCoordinates(2) / z_f
    // };
    // int internal_coord_local[3] = {
    //     vct->getCoordinates(0) % x_f,
    //     vct->getCoordinates(1) % y_f,
    //     vct->getCoordinates(2) % z_f
    // };

    string name_file = generate_file_name(coords);
    // stringstream restart_name;
    // restart_name << RestartDirName << "/restart_x"
    //              << setw(3) << setfill('0') << file_coord[0] << "_y"
    //              << setw(3) << setfill('0') << file_coord[1] << "_z"
    //              << setw(3) << setfill('0') << file_coord[2] << ".hdf";
    // string name_file = restart_name.str();

    // // 2. Open File and Get Last Cycle
    // 
    // if (file_id < 0) {
    //     eprintf("couldn't open file: %s\n\tRESTART NOT POSSIBLE", name_file.c_str());
    // }

    hid_t file_id = H5Fopen(name_file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

    //NB: useless as above
    // int last_cycle = 0;
    // hid_t dataset_id = H5Dopen2(file_id, "/last_cycle", H5P_DEFAULT);
    // H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &last_cycle);
    // H5Dclose(dataset_id);

    // 3. Setup Shared Dataspace and Dimensions (Using 'x' dataset)
    stringstream ss;
    ss << "/particles/species_" << species_number << "/x/cycle_" << last_cycle;
    hid_t dataset_id = H5Dopen2(file_id, ss.str().c_str(), H5P_DEFAULT);

    hid_t shared_dataspace = H5Dget_space(dataset_id);
    hsize_t dims_out[1];
    H5Sget_simple_extent_dims(shared_dataspace, dims_out, NULL);
    H5Dclose(dataset_id); // Close the dataset, we just needed the dataspace to get dims

    // 4. Calculate 1D chunking
    int tot_mult_fact = x_f * y_f * z_f;
    int chunk_id = coords.internal[0]
                 + (coords.internal[1] * x_f)
                 + (coords.internal[2] * x_f * y_f);

    hsize_t count[1] = { dims_out[0] / tot_mult_fact };
    if (chunk_id == tot_mult_fact - 1) {
        count[0] += dims_out[0] % tot_mult_fact; // Give remainder to the last chunk
    }

    hsize_t offset[1] = { (hsize_t)chunk_id * (dims_out[0] / tot_mult_fact) };

    // 5. Memory Allocation (Preserving original padding logic)
    const int nop        = count[0];
    const int padded_nop = roundup_to_multiple(nop, DVECWIDTH);

    u.reserve(padded_nop); v.reserve(padded_nop); w.reserve(padded_nop);
    q.reserve(padded_nop); x.reserve(padded_nop); y.reserve(padded_nop);
    z.reserve(padded_nop); t.reserve(padded_nop);

    u.resize(nop); v.resize(nop); w.resize(nop);
    q.resize(nop); x.resize(nop); y.resize(nop);
    z.resize(nop); t.resize(nop);

    // 6. Setup 1D Hyperslab
    hid_t memspace_id = H5Screate_simple(1, count, NULL);
    H5Sselect_hyperslab(shared_dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);

    // 7. Lambda helper to eliminate repetitive reads
    auto read_particle_array = [&](const string& path, vector_double& vec) {
        hid_t d_id = H5Dopen2(file_id, path.c_str(), H5P_DEFAULT);
        if (d_id >= 0) { // Safety check
            H5Dread(d_id, H5T_NATIVE_DOUBLE, memspace_id, shared_dataspace, H5P_DEFAULT, &vec[0]);
            H5Dclose(d_id);
        }
    };

    // 8. Execute Reads (Reusing the shared hyperslab for all arrays)
    ss.str(""); ss << "/particles/species_" << species_number << "/x/cycle_" << last_cycle;
    read_particle_array(ss.str(), x);

    ss.str(""); ss << "/particles/species_" << species_number << "/y/cycle_" << last_cycle;
    read_particle_array(ss.str(), y);

    ss.str(""); ss << "/particles/species_" << species_number << "/z/cycle_" << last_cycle;
    read_particle_array(ss.str(), z);

    ss.str(""); ss << "/particles/species_" << species_number << "/u/cycle_" << last_cycle;
    read_particle_array(ss.str(), u);

    ss.str(""); ss << "/particles/species_" << species_number << "/v/cycle_" << last_cycle;
    read_particle_array(ss.str(), v);

    ss.str(""); ss << "/particles/species_" << species_number << "/w/cycle_" << last_cycle;
    read_particle_array(ss.str(), w);

    ss.str(""); ss << "/particles/species_" << species_number << "/q/cycle_" << last_cycle;
    read_particle_array(ss.str(), q);

    // Following original code logic: read 'q' into 't' as fallback for particle ID
    ss.str(""); ss << "/particles/species_" << species_number << "/q/cycle_" << last_cycle;
    read_particle_array(ss.str(), t);

    // 9. Cleanup
    H5Sclose(shared_dataspace);
    H5Sclose(memspace_id);
    H5Fclose(file_id);
#endif
}
