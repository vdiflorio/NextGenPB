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

#include "vtk_class.h"

void
VTKWriter:: writeFieldVtuBinary (tmesh_3d &tmsh,
                                 const distributed_vector &field,
                                 std::string &fieldname_)
{
  preparingData (field, tmsh);
  setFieldName (fieldname_);
  initializeFile();
  writePieceHeader();
  writeGrid();
  writeDataPoints();
  finalizeFile();
}

void
VTKWriter:: writeFieldVtuBinary_local (tmesh_3d &tmsh,
                                       const distributed_vector &field,
                                       std::string &fieldname_, double *l_cr, double *r_cr)
{
  preparingData_local (field, tmsh, l_cr, r_cr);
  setFieldName (fieldname_);
  initializeFile();
  writePieceHeader();
  writeGrid();
  writeDataPoints();
  finalizeFile();
}

void
VTKWriter::createPvtuFile (std::vector<std::string> &fieldNames, std::vector<std::string> &baseNames, int numRanks)
{
  // Build the .pvtu filename
  if (fieldNames.size() != baseNames.size()) {
    throw std::runtime_error ("Error: different number of field and names! ");
  }

  for (size_t i = 0; i < baseNames.size(); ++i) {

    std::string pvtuFilename = "total_" + baseNames[i] + ".pvtu";

    // Open the .pvtu file for writing
    std::ofstream pvtuFile (pvtuFilename);

    if (!pvtuFile.is_open()) {
      throw std::runtime_error ("Error: impossible to create the file " + pvtuFilename);
    }

    // Write the XML header
    pvtuFile << "<?xml version=\"1.0\"?>\n";
    pvtuFile << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    pvtuFile << "  <PUnstructuredGrid>\n";

    // Write PPoints, PCells and PPointData sections
    pvtuFile << "    <PPoints>\n";
    pvtuFile << "      <PDataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"Array\"/>\n";
    pvtuFile << "    </PPoints>\n";

    pvtuFile << "    <PCells>\n";
    pvtuFile << "      <PDataArray type=\"Int32\" Name=\"connectivity\"/>\n";
    pvtuFile << "      <PDataArray type=\"Int32\" Name=\"offsets\"/>\n";
    pvtuFile << "      <PDataArray type=\"Int32\" Name=\"types\"/>\n";
    pvtuFile << "    </PCells>\n";

    // Add point data for each field
    pvtuFile << "    <PPointData>\n";
    pvtuFile << "      <PDataArray type=\"Float64\" Name=\"" << fieldNames[i] << "\"/>\n";
    pvtuFile << "    </PPointData>\n";

    // Add per-rank .vtu piece references
    for (int rank = 0; rank < numRanks; ++rank)
      pvtuFile << "    <Piece Source=\""
               << baseNames[i] << "_" << std::setfill ('0') << std::setw (4) << rank << ".vtu\"/>\n";

    // Close the root sections
    pvtuFile << "  </PUnstructuredGrid>\n";
    pvtuFile << "</VTKFile>\n";

    pvtuFile.close();
    std::cout << ".pvtu file created: " << pvtuFilename << std::endl;
  }
}

void
VTKWriter::preparingData (const distributed_vector &field, tmesh_3d &tmsh)
{
  int real_node = tmsh.num_owned_nodes();
  nelems = tmsh.num_local_quadrants ();

  p.assign (3 * real_node, 0);
  f_loc.assign (real_node, 0);
  t.assign (8 * nelems, 0);

  int ij = 0;
  int triang = 0;

  for (auto quadrant = tmsh.begin_quadrant_sweep ();
       quadrant != tmsh.end_quadrant_sweep ();
       ++quadrant) {
    ij = 0;

    for (int ii = 0; ii < 8; ++ii) {
      triang = quadrant->t (ii);

      if (! quadrant->is_hanging (ii))
        if (triang < real_node) {
          for (int jj = 0; jj < 3; ++jj)
            p[3 * triang + jj] = quadrant->p (jj, ii);

          f_loc[triang] = field[quadrant->gt (ii)];
          t[8 * quadrant->get_forest_quad_idx () + (ij++)] = triang;
        } else {
          for (int jj = 0; jj < 3; ++jj)
            p.push_back (quadrant->p (jj, ii));

          f_loc.push_back (field[quadrant->gt (ii)]);
          t[8 * quadrant->get_forest_quad_idx () + (ij++)] =
            f_loc.size () - 1;
        } else {
        for (int jj = 0; jj < 3; ++jj)
          p.push_back (quadrant->p (jj, ii));

        int pp;
        double fbuff = 0;

        for (pp = 0; pp < quadrant->num_parents (ii); ++pp)
          fbuff += field [quadrant->gparent (pp, ii)];

        f_loc.push_back (fbuff / pp);
        t[8 * quadrant->get_forest_quad_idx () + (ij++)] =
          f_loc.size () - 1;
      }
    }
  }

  nnodes = f_loc.size ();
}

void VTKWriter::preparingData_local (const distributed_vector &field,
                                     tmesh_3d &tmsh,
                                     double *l_cr,
                                     double *r_cr)
{
  nelems = 0;
  nnodes = 0;

  // Index map from global node index to local VTK node index; -1 = not yet added
  std::vector<int> global_to_local (tmsh.num_owned_nodes(), -1);

  for (auto quadrant = tmsh.begin_quadrant_sweep();
       quadrant != tmsh.end_quadrant_sweep();
       ++quadrant) {

    // Corner coordinates of the current quadrant
    double x1 = quadrant->p (0, 0), y1 = quadrant->p (1, 0), z1 = quadrant->p (2, 0);
    double x2 = quadrant->p (0, 7), y2 = quadrant->p (1, 7), z2 = quadrant->p (2, 7);

    // Skip quadrants that lie outside the requested region
    if ((x1 > l_cr[0]) && (x2 < r_cr[0]) &&
        (y1 > l_cr[1]) && (y2 < r_cr[1]) &&
        (z1 > l_cr[2]) && (z2 < r_cr[2])) {
      // Iterate over the 8 corner nodes of the quadrant
      for (int ii = 0; ii < 8; ++ii) {
        int global_index = quadrant->gt (ii);

        // Check whether this node has already been added to the local list
        if (global_to_local[global_index] == -1) {
          if (rank == 1) {
            std::cout << f_loc.size() << " " << global_to_local[global_index] << "  " << global_index << std::endl;
          }

          // New node: append its coordinates and update the global→local map
          for (int jj = 0; jj < 3; ++jj) {
            p.push_back (quadrant->p (jj, ii));
          }

          // Compute the field value at this node
          double field_value = 0.0;

          if (!quadrant->is_hanging (ii)) {
            field_value = field[global_index];
          } else {
            // Hanging node: average over parent node values
            int num_parents = quadrant->num_parents (ii);

            for (int pp = 0; pp < num_parents; ++pp) {
              field_value += field[quadrant->gparent (pp, ii)];
            }

            field_value /= num_parents;
          }

          f_loc.push_back (field_value);

          // Record the new local index for this node
          global_to_local[global_index] = nnodes++;
        }

        // Append the local node index to the connectivity vector
        t.push_back (global_to_local[global_index]);
      }

      // Incrementa il numero di elementi
      ++nelems;
    }
  }

  // Update the total node count from the field value vector
  nnodes = f_loc.size();
}

void
VTKWriter::setBaseName (std::string basename_, int rank_)
{
  // Fall back to "default" when no basename is provided
  basename = basename_.empty() ? "default" : basename_;
  rank = rank_;

  // Build the per-rank .vtu filename
  snprintf (filename, sizeof (filename), "%s_%4.4d.vtu", basename.c_str(), rank);

  // Reset the byte offset counter whenever the filename changes
  offsetCounter = 0;
  // Reset the binary data container and associated buffers
  ByteArray().swap (dataContainer);
  int nnodes;
  int nelems;
  std::vector<double>().swap (p);
  std::vector<double>().swap (f_loc);
  std::vector<int>().swap (t);
}

void
VTKWriter::setFieldName (std::string &fieldname_)
{
  if (fieldname_.empty()) {
    fieldname = "field";
    std::cerr << "Warning: fieldname is empty. Setting basename to 'fied'.\n";
  } else
    fieldname = fieldname_;
}


void
VTKWriter::initializeFile()
{
  // Open the file in binary write mode, truncating any existing content
  file.open (filename, std::ios::binary | std::ios::trunc);

  if (!file.is_open()) {
    throw std::runtime_error ("Error: cannot open file " + std::string (filename));
  }

  // Write the VTK XML file header
  file << "<?xml version=\"1.0\"?>\n";
  file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
  file << "<UnstructuredGrid>\n";
}


void
VTKWriter::writePieceHeader()
{
  file << "<Piece NumberOfPoints=\"" << nnodes << "\" NumberOfCells=\"" << nelems << "\">\n";
}


template <typename T>
ByteArray serializeData (const std::vector<T> &data)
{
  ByteArray rawData;
  rawData.reserve (data.size() * sizeof (T));

  for (const T& value : data) {
    const char* rawBytes = reinterpret_cast<const char *> (&value);
    rawData.insert (rawData.end(), rawBytes, rawBytes + sizeof (T));
  }

  return rawData;
}

template <typename T>
ByteArray serializeDataWithHeader (const std::vector<T> &data)
{
  ByteArray rawData;
  // Prepend the block size (in bytes) as a 4-byte header, as required by VTK appended format
  int32_t blockSize = static_cast<int32_t> (data.size() * sizeof (T));
  const char* headerBytes = reinterpret_cast<const char *> (&blockSize);
  rawData.insert (rawData.end(), headerBytes, headerBytes + sizeof (int32_t));

  // Append the actual data bytes
  const char* dataBytes = reinterpret_cast<const char *> (data.data());
  rawData.insert (rawData.end(), dataBytes, dataBytes + blockSize);

  return rawData;
}

void
VTKWriter::appendData (const ByteArray &rawData)
{
  // Append binary data to the main container
  dataContainer.insert (dataContainer.end(), rawData.begin(), rawData.end());

  // Update the running byte offset (total bytes written so far)
  offsetCounter += rawData.size();
  // size_t alignTo = 8; // Multiplo di 8 per l'allineamento
  // offsetCounter = (offsetCounter + alignTo - 1) & ~(alignTo - 1);
}

void
VTKWriter::writeDataArray (const std::string &type, const std::string &name, int numComponents,
                           const std::vector<double> &data, int numEntries)
{
  file << "<DataArray type=\"" << type << "\" Name=\"" << name
       << "\" NumberOfComponents=\"" << numComponents
       << "\" format=\"appended\" offset=\"" << offsetCounter << "\">\n";

  // ByteArray rawData = serializeData(data);
  ByteArray rawData = serializeDataWithHeader (data);
  appendData (rawData);

  file << "</DataArray>\n";
}

void VTKWriter::writeConnectivity()
{
  file << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\"" << offsetCounter << "\">\n";
  file << "</DataArray>\n";

  // ByteArray rawData = serializeData(t);
  ByteArray rawData = serializeDataWithHeader (t); // Int32 connectivity with VTK block header
  appendData (rawData);
}

void VTKWriter::writeOffsets()
{
  std::vector<int> offsets (nelems);
  int nodesPerCell = 8; // 3D hexahedral cells

  for (int i = 0; i < nelems; ++i) {
    offsets[i] = (i + 1) * nodesPerCell;
  }

  file << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\"" << offsetCounter << "\">\n";
  file << "</DataArray>\n";

  // ByteArray rawData = serializeData(offsets);
  ByteArray rawData = serializeDataWithHeader (offsets);
  appendData (rawData);
}

void VTKWriter::writeCellTypes()
{
  int cellType = 11; // VTK_VOXEL (3D hexahedral cell type)
  std::vector<int> types (nelems, cellType);

  file << "<DataArray type=\"Int32\" Name=\"types\" format=\"appended\" offset=\"" << offsetCounter << "\">\n";
  file << "</DataArray>\n";

  // ByteArray rawData = serializeData(types);
  ByteArray rawData = serializeDataWithHeader (types);
  appendData (rawData);
}

void
VTKWriter::writeGrid()
{
  file << "<Points>\n";
  writeDataArray ("Float64", "Array", 3, p, nnodes); // Node coordinates
  file << "</Points>\n";

  file << "<Cells>\n";
  writeConnectivity(); // Cell-to-node connectivity
  writeOffsets(); // Per-cell node count offsets
  writeCellTypes(); // VTK cell type codes
  file << "</Cells>\n";
}

void
VTKWriter::writeDataPoints()
{

  file << "<PointData>\n";
  writeDataArray ("Float64", fieldname, 1, f_loc, nnodes);

  file << "</PointData>\n";

}



void
VTKWriter::finalizeFile()
{
  file << "</Piece>\n";
  file << "</UnstructuredGrid>\n";
  file << " <AppendedData encoding=\"raw\">\n";
  file << "_";
  file.write (reinterpret_cast<const char *> (dataContainer.data()), dataContainer.size());
  file << " </AppendedData>\n";
  file << "</VTKFile>\n";
  file.close();
}

bool
VTKWriter::fileExists (std::string &name)
{
  std::ifstream f (name.c_str());
  return f.good();
}
