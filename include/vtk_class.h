/*
 *  Copyright (C) 2024-2025 Vincenzo Di Florio
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

/**
 * @file vtk_class.h
 * @brief VTK UnstructuredGrid (VTU) writer for distributed octree fields.
 *
 * Provides VTKWriter, which serialises a BIM++ distributed_vector field on a
 * tmesh_3d octree to a binary VTU file per MPI rank and an accompanying
 * parallel PVTU descriptor.  Output is compatible with ParaView and VisIt.
 */
#ifndef VTKCLASS_H
#define VTKCLASS_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdint>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <bim_distributed_vector.h>
#include <quad_operators_3d.h>
#include <tmesh_3d.h>

typedef std::vector<uint8_t> ByteArray; ///< Raw byte buffer for binary VTK data

/**
 * @brief Writes a scalar field on an adaptive octree to binary VTU/PVTU files.
 *
 * Each MPI rank writes its own `<basename>_XXXX.vtu` file.  After all ranks
 * have written, rank 0 creates a `<basename>.pvtu` that links them together.
 *
 * Typical usage:
 * @code
 *   VTKWriter writer("phi", rank);
 *   writer.writeFieldVtuBinary(tmsh, *phi, fieldname);
 *   writer.createPvtuFile(fieldNames, baseNames, size);
 * @endcode
 */
class VTKWriter
{
public:

  /// @brief Construct a VTKWriter for the given base filename and MPI rank.
  /// @param basename_ Base name for output files (without extension).
  /// @param rank_     MPI rank of this process.
  VTKWriter (std::string basename_ = "default", int rank_ = 0)
    : basename (basename_.empty() ? "default" : basename_),
      rank (rank_), offsetCounter (0)
  {
    snprintf (filename, sizeof (filename), "%s_%4.4d.vtu", basename.c_str(), rank);
  }

  /// @brief Change the output base name and rank after construction.
  void setBaseName (std::string basename_, int rank_);

  /// @brief Write a distributed scalar field on @p tmsh to a binary VTU file.
  /// @param tmsh       The adaptive octree mesh.
  /// @param field      The distributed scalar field to write.
  /// @param fieldname_ Name tag written inside the VTU XML.
  void writeFieldVtuBinary (tmesh_3d& tmsh,
                            const distributed_vector& field,
                            std::string &fieldname_);

  /// @brief Write only the sub-region [@p l_cr, @p r_cr] of the field to a VTU file.
  /// @param tmsh       The adaptive octree mesh.
  /// @param field      The distributed scalar field to write.
  /// @param fieldname_ Name tag written inside the VTU XML.
  /// @param l_cr       Lower corner of the sub-region [Å].
  /// @param r_cr       Upper corner of the sub-region [Å].
  void writeFieldVtuBinary_local (tmesh_3d& tmsh,
                                  const distributed_vector& field,
                                  std::string &fieldname_, double* l_cr, double* r_cr);

  /// @brief Create the parallel PVTU metadata file that links all per-rank VTU files.
  /// @param fieldNames List of field names included in this dataset.
  /// @param baseNames  Corresponding VTU base names.
  /// @param numRanks   Total number of MPI ranks.
  void createPvtuFile (std::vector<std::string> &fieldNames, std::vector<std::string> &baseNames, int numRanks);

private:
  std::string basename;   ///< Base filename (no rank suffix, no extension)
  std::string fieldname;  ///< Current field name tag
  char filename[255];     ///< Full per-rank filename (e.g. "phi_0001.vtu")
  std::ofstream file;     ///< Output file stream
  ByteArray dataContainer;///< Binary data buffer for appended VTK section
  int nnodes;             ///< Number of nodes in the local piece
  int nelems;             ///< Number of elements in the local piece
  std::vector<double> p;  ///< Node coordinates (x0,y0,z0, x1,y1,z1, …)
  std::vector<double> f_loc; ///< Local field values at nodes
  std::vector<int> t;     ///< Element connectivity
  int offsetCounter;      ///< Running byte offset for the appended data section
  int rank;               ///< MPI rank

  /// @brief Extract and sort nodal coordinates and field values from the mesh.
  void preparingData (const distributed_vector& field, tmesh_3d& tmsh);

  /// @brief Extract nodal data restricted to the sub-region [@p l_cr, @p r_cr].
  void preparingData_local (const distributed_vector& field, tmesh_3d& tmsh, double* l_cr, double* r_cr);

  void setFieldName (std::string &fieldname_); ///< Store field name for VTK header
  void initializeFile();      ///< Write VTU XML header and open appended data section
  void writePieceHeader();    ///< Write VTK Piece element with point/cell counts
  void writeGrid();           ///< Write Points + Cells sections
  void writeDataPoints();     ///< Write node coordinate DataArray
  /// @brief Write a named DataArray with binary offset reference.
  void writeDataArray (const std::string& type, const std::string& name, int numComponents,
                       const std::vector<double>& data, int numEntries);
  void writeConnectivity();   ///< Write element connectivity DataArray
  void writeOffsets();        ///< Write element offset DataArray
  void writeCellTypes();      ///< Write element type DataArray (all VTK_HEXAHEDRON)
  void appendData (const ByteArray& rawData); ///< Append binary data block
  void finalizeFile();        ///< Close VTU XML and flush file
  bool fileExists (std::string& name); ///< Return true if @p name already exists on disk
};

#endif //VTKCLASS_H

