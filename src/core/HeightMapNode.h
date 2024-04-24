/*
 *  OpenSCAD (www.openscad.org)
 *  Copyright (C) 2009-2011 Clifford Wolf <clifford@clifford.at> and
 *                          Marius Kintel <marius@kintel.net>
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  As a special exception, you have permission to link this program
 *  with the CGAL library and distribute executables, as long as you
 *  follow the requirements of the GNU GPL in regard to all of the
 *  software in the executable aside from CGAL.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

#include "FunctionType.h"
#include "Value.h"
#include "node.h"
#include "linalg.h"
#include <memory>

struct map_data_t
{
public:
  using storage_type = double; // float could be enough here

  map_data_t() { height = width = 0; }

  void clear() { height = width = 0; storage.clear(); }

  void reserve(size_t x) { storage.reserve(x); }

  void resize(size_t x) { storage.resize(x); }

  storage_type& operator[](int x) { return storage[x]; }
  const storage_type& operator[](int x) const { return storage[x]; }

public:
  unsigned int height; // rows
  unsigned int width; // columns
  std::vector<storage_type> storage;

};


class HeightMapNode : public LeafNode
{
public:
  VISITABLE();
  HeightMapNode(const ModuleInstantiation *mi) : LeafNode(mi) { }
  std::string toString() const override;
  std::string name() const override { return "heightmap"; }

  Filename filename;
  Vector3d size{1,1,1};
  bool center{false};
  bool doubleSided{false};
  float thickness{1};
  int convexity{1};
  int pixelStep{1};
  std::array<double, 4> fadeTo;
  std::array<double, 4> fadeWidth;
  std::string shapeExpr;
  std::vector<Vector3d> shapeGridPoints;
  std::vector<Vector3d> shapeGridNormals;

  std::pair<unsigned int, unsigned int> getDataSize(std::string filename) const;
  std::unique_ptr<const Geometry> createGeometry() const override;
private:
  void convert_image(map_data_t& data, std::vector<uint8_t>& img, unsigned int width, unsigned int height) const;
  bool is_png(std::vector<uint8_t>& img) const;
  map_data_t read_dat(std::string filename) const;
  map_data_t read_png_or_dat(std::string filename) const;
};

