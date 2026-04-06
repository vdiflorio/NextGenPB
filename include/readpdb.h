/**
 * @file readpdb.h
 * @brief PDB/PQR file I/O utilities for reading atomic structures.
 *
 * Provides functions to read PDB and PQR files, load van der Waals radii and
 * partial charges from parameter files, and write PQR output files.
 */
#ifndef READPDB_H
#define READPDB_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <array>
#include <cctype>
#include <cstdlib>

#include <nanoshaper.h>

/// @brief Read atomic coordinates (ATOM/HETATM records) from a PDB file.
/// @param filename Path to the PDB file.
/// @param atoms    Output atom list (name, residue, chain, position populated).
void read_pdb (const std::string& filename, std::vector<NS::Atom>& atoms);

/// @brief Load van der Waals radii from a parameter file and assign them to @p atoms.
/// @param filename Path to the radius parameter file.
/// @param atoms    Atom list to update (matched by atom name / residue / chain).
void load_radii (const std::string& filename, std::vector<NS::Atom>& atoms);

/// @brief Load partial charges from a parameter file and assign them to @p atoms.
/// @param filename Path to the charge parameter file.
/// @param atoms    Atom list to update.
void load_charges (const std::string& filename, std::vector<NS::Atom>& atoms);

/// @brief Write atom data (coordinates, charges, radii) to a PQR file.
/// @param filename Path to the output PQR file.
/// @param atoms    Atom list to write.
void write_pqr (const std::string& filename, const std::vector<NS::Atom>& atoms);

/// @brief Parse a single line from a radius parameter file and update @p atom.
/// @param line Input line string.
/// @param atom Atom to update if the line matches.
/// @return Number of fields successfully parsed (0 = no match).
int parse_line_to_atom_radius (const std::string& line, NS::Atom& atom);

/// @brief Parse a single line from a charge parameter file and update @p atom.
/// @param line Input line string.
/// @param atom Atom to update if the line matches.
/// @return Number of fields successfully parsed (0 = no match).
int parse_line_to_atom_charge (const std::string& line, NS::Atom& atom);

#endif
