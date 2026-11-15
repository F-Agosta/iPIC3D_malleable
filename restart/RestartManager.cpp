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

#include "RestartManager.h"
#include "ipichdf5.h"
#include "errors.h"
#include "ipicdefs.h"
#include "ipicmath.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <algorithm>

using namespace std;

// ============================================================
// Constructor
// ============================================================

/*! The RestartState is fully populated by the caller (Collective::
 *  read_from_settings); the manager is immutable after construction. */
RestartManager::RestartManager(const string& restartDirName, const RestartState& state_)
    : RestartDirName(restartDirName),
      state(state_),
      restart_state(1)
{
}

// ============================================================
// file_name
// ============================================================

string RestartManager::file_name(const int coord[3]) const
{
    stringstream ss;
    ss << RestartDirName << "/restart_x"
       << setw(3) << setfill('0') << coord[0] << "_y"
       << setw(3) << setfill('0') << coord[1] << "_z"
       << setw(3) << setfill('0') << coord[2] << ".hdf";
    return ss.str();
}

// ============================================================
// compute_slabs
// ============================================================

/*! Intersection of this process's new subdomain with the old subdomains.
 *  Per direction d, the new subdomain of coordinate a covers the global
 *  cells [a*w_new, (a+1)*w_new) with w_new = nc[d]/currLEN[d]; the old
 *  subdomain i covers [i*w_old, (i+1)*w_old) with w_old = nc[d]/pastLEN[d].
 *  For every old subdomain i overlapping the new one, the intersection
 *  [lo, hi) gives the file slab (offset lo - i*w_old inside file i) and the
 *  memory slab (offset lo - a*w_new inside the local grid), both with
 *  count hi - lo. The 3D slabs are the Cartesian product over directions.
 *  Growing, shrinking and same-size restarts are all just special cases of
 *  this intersection: per direction, the new subdomain overlaps exactly one
 *  old subdomain when currLEN >= pastLEN, and w_old/w_new of them when
 *  currLEN < pastLEN. */
vector<RestartManager::RestartSlab>
RestartManager::compute_slabs(const VCtopology3D* vct) const
{
    struct DirSlab {
        int old_idx;  // index of the old subdomain along this direction
        int file_off; // offset of the intersection inside the old subdomain (cells)
        int mem_off;  // offset of the intersection inside the new subdomain (cells)
        int count;    // extent of the intersection (cells)
    };

    vector<DirSlab> dir[3];
    for (int d = 0; d < 3; d++) {
        const int G     = state.nc[d];
        const int p     = state.pastLEN[d];
        const int n     = state.currLEN[d];
        const int w_old = G / p;
        const int w_new = G / n;
        const int a     = vct->getCoordinates(d);
        const int new_lo = a * w_new;
        const int new_hi = (a + 1) * w_new;
        const int i_min = new_lo / w_old;
        const int i_max = (new_hi - 1) / w_old;
        for (int i = i_min; i <= i_max; i++) {
            const int lo = max(new_lo, i * w_old);
            const int hi = min(new_hi, (i + 1) * w_old);
            DirSlab s;
            s.old_idx  = i;
            s.file_off = lo - i * w_old;
            s.mem_off  = lo - new_lo;
            s.count    = hi - lo;
            dir[d].push_back(s);
        }
    }

    vector<RestartSlab> slabs;
    for (const auto& sx : dir[0])
        for (const auto& sy : dir[1])
            for (const auto& sz : dir[2]) {
                RestartSlab rs;
                const DirSlab* ds[3] = {&sx, &sy, &sz};
                for (int d = 0; d < 3; d++) {
                    rs.file.file_coord[d]  = ds[d]->old_idx;
                    rs.file.slab.offset[d] = ds[d]->file_off;
                    rs.file.slab.count[d]  = ds[d]->count;
                    rs.mem.slab.offset[d]  = ds[d]->mem_off;
                    rs.mem.slab.count[d]   = ds[d]->count;
                }
                slabs.push_back(rs);
            }
    return slabs;
}

// ============================================================
// owns_particle
// ============================================================

/*! A particle belongs to this process if its global cell index falls
 *  inside this process's new subdomain along every direction. */
bool RestartManager::owns_particle(double px, double py, double pz,
                                   const VCtopology3D* vct) const
{
    const double pos[3] = {px, py, pz};
    for (int d = 0; d < 3; d++) {
        const double cell_size = state.L[d] / state.nc[d];
        const double cell = pos[d] / cell_size;
        const int a = vct->getCoordinates(d);
        const int w_new = state.nc[d] / state.currLEN[d];
        if (cell < a * w_new || cell >= (a + 1) * w_new)
            return false;
    }
    return true;
}

// ============================================================
// read_last_cycle
// ============================================================

/*! Read the last cycle number stored in the rank-0 restart file. */
int RestartManager::read_last_cycle() const
{
#ifdef NO_HDF5
    eprintf("restart requires compiling with HDF5");
    return -1;
#else
    const int coord[3] = {0, 0, 0};
    const string name_file = file_name(coord);

    hid_t file_id = H5Fopen(name_file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0) {
        cout << "couldn't open restart file for last_cycle: " << name_file << endl;
        return -1;
    }

    int last_cycle = state.last_cycle;
    hid_t dataset_id = H5Dopen2(file_id, "/last_cycle", H5P_DEFAULT);
    if (dataset_id >= 0) {
        H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &last_cycle);
        H5Dclose(dataset_id);
    }
    H5Fclose(file_id);
    return last_cycle;
#endif
}

// ============================================================
// read_field_restart
// ============================================================

/*! Read electromagnetic fields and charge densities from the restart files.
 *  For each old subdomain intersecting this process's new subdomain, read
 *  the file slab and scatter it into the corresponding memory slab.
 *  Corresponds to the old Collective::read_field_restart(). */
void RestartManager::read_field_restart(
    const VCtopology3D* vct,
    const Grid* grid,
    arr3_double Bxn, arr3_double Byn, arr3_double Bzn,
    arr3_double Ex,  arr3_double Ey,  arr3_double Ez,
    array4_double* rhons_, int ns) const
{
#ifdef NO_HDF5
    eprintf("Require HDF5 to read from restart file.");
#else
    if (vct->getCartesian_rank() == 0) {
        printf("LOADING EM FIELD FROM RESTART FILE in %s\n", RestartDirName.c_str());
    }

    const int last_cycle = state.last_cycle;
    const vector<RestartSlab> slabs = compute_slabs(vct);

    for (const auto& rs : slabs) {
        const string name_file = file_name(rs.file.file_coord);
        hid_t file_id = H5Fopen(name_file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        if (file_id < 0) {
            eprintf("Failed to open file: %s\n", name_file.c_str());
        }

        // field datasets store nodes, so a slab of `count` cells spans
        // `count + 1` nodes; the left edges coincide
        hsize_t offset[3], count[3];
        for (int d = 0; d < 3; d++) {
            offset[d] = (hsize_t)rs.file.slab.offset[d];
            count[d]  = (hsize_t)rs.file.slab.count[d] + 1;
        }

        vector<double> buf(count[0] * count[1] * count[2]);
        hid_t memspace = H5Screate_simple(3, count, NULL);

        // read one field dataset and scatter it into the memory slab
        auto read_and_map = [&](const string& path, arr3_double arr) {
            hid_t d_id = H5Dopen2(file_id, path.c_str(), H5P_DEFAULT);
            hid_t fspace = H5Dget_space(d_id);
            H5Sselect_hyperslab(fspace, H5S_SELECT_SET, offset, NULL, count, NULL);
            H5Dread(d_id, H5T_NATIVE_DOUBLE, memspace, fspace, H5P_DEFAULT, buf.data());
            int k = 0;
            for (int i = 0; i < (int)count[0]; i++)
                for (int j = 0; j < (int)count[1]; j++)
                    for (int jj = 0; jj < (int)count[2]; jj++)
                        arr[rs.mem.slab.offset[0] + i]
                           [rs.mem.slab.offset[1] + j]
                           [rs.mem.slab.offset[2] + jj] = buf[k++];
            H5Sclose(fspace);
            H5Dclose(d_id);
        };

        stringstream ss;
        ss << "/fields/Bx/cycle_" << last_cycle; read_and_map(ss.str(), Bxn);
        ss.str(""); ss << "/fields/By/cycle_" << last_cycle; read_and_map(ss.str(), Byn);
        ss.str(""); ss << "/fields/Bz/cycle_" << last_cycle; read_and_map(ss.str(), Bzn);
        ss.str(""); ss << "/fields/Ex/cycle_" << last_cycle; read_and_map(ss.str(), Ex);
        ss.str(""); ss << "/fields/Ey/cycle_" << last_cycle; read_and_map(ss.str(), Ey);
        ss.str(""); ss << "/fields/Ez/cycle_" << last_cycle; read_and_map(ss.str(), Ez);

        // charge density (rho) for each species
        array4_double& rhons = *rhons_;
        for (int is = 0; is < ns; is++) {
            ss.str(""); ss << "/moments/species_" << is << "/rho/cycle_" << last_cycle;
            hid_t d_id = H5Dopen2(file_id, ss.str().c_str(), H5P_DEFAULT);
            hid_t fspace = H5Dget_space(d_id);
            H5Sselect_hyperslab(fspace, H5S_SELECT_SET, offset, NULL, count, NULL);
            H5Dread(d_id, H5T_NATIVE_DOUBLE, memspace, fspace, H5P_DEFAULT, buf.data());
            int k = 0;
            for (int i = 0; i < (int)count[0]; i++)
                for (int j = 0; j < (int)count[1]; j++)
                    for (int jj = 0; jj < (int)count[2]; jj++)
                        rhons[is][rs.mem.slab.offset[0] + i]
                                [rs.mem.slab.offset[1] + j]
                                [rs.mem.slab.offset[2] + jj] = buf[k++];
            H5Sclose(fspace);
            H5Dclose(d_id);
        }

        H5Sclose(memspace);
        H5Fclose(file_id);
    }
#endif
}

// ============================================================
// read_particles_restart
// ============================================================

/*! Read particle data for one species from the restart files.
 *  Each intersecting old file is read as a 1D slab (the whole particle
 *  array); every particle is then checked against this process's subdomain
 *  and only the ones inside it are appended to the local vectors (the
 *  memory slab grows at the current end of the vectors).
 *  Corresponds to the old Collective::read_particles_restart(). */
void RestartManager::read_particles_restart(
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
        printf("LOADING PARTICLES FROM RESTART FILE in %s\n", RestartDirName.c_str());
    }

    const int last_cycle = state.last_cycle;
    const vector<RestartSlab> slabs = compute_slabs(vct);

    // distinct old files intersecting this process's subdomain
    vector<FileSlab> files;
    for (const auto& rs : slabs) {
        bool seen = false;
        for (const auto& f : files) {
            if (f.file_coord[0] == rs.file.file_coord[0] &&
                f.file_coord[1] == rs.file.file_coord[1] &&
                f.file_coord[2] == rs.file.file_coord[2]) {
                seen = true;
                break;
            }
        }
        if (!seen) {
            FileSlab f;
            for (int d = 0; d < 3; d++) {
                f.file_coord[d]  = rs.file.file_coord[d];
                f.slab.offset[d] = 0;
                f.slab.count[d]  = 1;
            }
            f.slab.count[0] = 0; // set to the array length once the file is open
            files.push_back(f);
        }
    }

    for (auto& fs : files) {
        const string name_file = file_name(fs.file_coord);
        hid_t file_id = H5Fopen(name_file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        if (file_id < 0) {
            eprintf("couldn't open file: %s\n\tRESTART NOT POSSIBLE", name_file.c_str());
        }

        // 1D file slab: the whole particle array of this old file
        stringstream ss;
        ss << "/particles/species_" << species_number << "/x/cycle_" << last_cycle;
        hid_t dataset_id = H5Dopen2(file_id, ss.str().c_str(), H5P_DEFAULT);
        hid_t fspace = H5Dget_space(dataset_id);
        hsize_t dims[1];
        H5Sget_simple_extent_dims(fspace, dims, NULL);
        H5Dclose(dataset_id);

        fs.slab.count[0] = (int)dims[0];
        const hsize_t count[1]  = { dims[0] };
        const hsize_t offset[1] = { 0 };
        hid_t memspace = H5Screate_simple(1, count, NULL);
        H5Sselect_hyperslab(fspace, H5S_SELECT_SET, offset, NULL, count, NULL);

        vector<double> px(dims[0]), py(dims[0]), pz(dims[0]);
        vector<double> pu(dims[0]), pv(dims[0]), pw(dims[0]), pq(dims[0]);

        auto read_arr = [&](const string& path, vector<double>& vec) {
            hid_t d_id = H5Dopen2(file_id, path.c_str(), H5P_DEFAULT);
            if (d_id >= 0) {
                H5Dread(d_id, H5T_NATIVE_DOUBLE, memspace, fspace, H5P_DEFAULT, vec.data());
                H5Dclose(d_id);
            }
        };

        ss.str(""); ss << "/particles/species_" << species_number << "/x/cycle_" << last_cycle; read_arr(ss.str(), px);
        ss.str(""); ss << "/particles/species_" << species_number << "/y/cycle_" << last_cycle; read_arr(ss.str(), py);
        ss.str(""); ss << "/particles/species_" << species_number << "/z/cycle_" << last_cycle; read_arr(ss.str(), pz);
        ss.str(""); ss << "/particles/species_" << species_number << "/u/cycle_" << last_cycle; read_arr(ss.str(), pu);
        ss.str(""); ss << "/particles/species_" << species_number << "/v/cycle_" << last_cycle; read_arr(ss.str(), pv);
        ss.str(""); ss << "/particles/species_" << species_number << "/w/cycle_" << last_cycle; read_arr(ss.str(), pw);
        ss.str(""); ss << "/particles/species_" << species_number << "/q/cycle_" << last_cycle; read_arr(ss.str(), pq);

        // memory slab: append at the current end of the local vectors,
        // keeping only the particles that belong to this subdomain
        MemorySlab ms;
        ms.slab.offset[0] = x.size();
        ms.slab.offset[1] = 0;
        ms.slab.offset[2] = 0;
        for (int p = 0; p < (int)dims[0]; p++) {
            if (owns_particle(px[p], py[p], pz[p], vct)) {
                x.push_back(px[p]); y.push_back(py[p]); z.push_back(pz[p]);
                u.push_back(pu[p]); v.push_back(pv[p]); w.push_back(pw[p]);
                q.push_back(pq[p]);
                t.push_back(pq[p]); // 'q' as fallback for particle ID, as before
            }
        }
        ms.slab.count[0] = x.size() - ms.slab.offset[0];
        ms.slab.count[1] = 1;
        ms.slab.count[2] = 1;

        H5Sclose(fspace);
        H5Sclose(memspace);
        H5Fclose(file_id);
    }

    // pad the local vectors to the vector width
    const int nop = x.size();
    const int padded_nop = roundup_to_multiple(nop, DVECWIDTH);
    u.reserve(padded_nop); v.reserve(padded_nop); w.reserve(padded_nop);
    q.reserve(padded_nop); x.reserve(padded_nop); y.reserve(padded_nop);
    z.reserve(padded_nop); t.reserve(padded_nop);
#endif
}
