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
#include "module.h"
#include "ModuleInstantiation.h"
#include "core/node.h"
#include "PolySet.h"
#include "PolySetBuilder.h"
#include "Builtins.h"
#include "Children.h"
#include "Parameters.h"
#include "printutils.h"
#include "io/fileutils.h"
#include "handle_dep.h"
#include "ext/lodepng/lodepng.h"
#include "HeightMapNode.h"
#include "Expression.h"

#include <CGAL/Algebraic_kernel_for_spheres_2_3.h>
#include <CGAL/Distance_2/Line_2_Line_2.h>
#include <CGAL/Kernel/global_functions_3.h>
#include <Eigen/Sparse>
#include <algorithm>
#include <boost/spirit/home/support/common_terminals.hpp>
#include <cstddef>
#include <cstdint>
#include <cmath>
#include <memory>
#include <sstream>
#include <fstream>
#include <unordered_map>
#include <boost/functional/hash.hpp>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/assign/std/vector.hpp>
#include <utility>
using namespace boost::assign; // bring 'operator+=()' into scope

#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

#include <boost/variant.hpp>
#include <boost/variant/get.hpp>
#include <CGAL/Point_3.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Spherical_kernel_3.h>
#include <CGAL/Exact_spherical_kernel_3.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>

using Kernel = CGAL::Cartesian<double>;
using Point = Kernel::Point_3;
using Point_2 = Kernel::Point_2;
using Plane_3 = Kernel::Plane_3;
using Vector_3 = Kernel::Vector_3;
using Vector_2 = Kernel::Vector_2;
using Line_2 = Kernel::Line_2;
using Line_3 = Kernel::Line_3;
using Triangle = Kernel::Triangle_3;
using Segment = Kernel::Segment_3;
using Ray = Kernel::Ray_3;
using Intersection_p2p = boost::optional<boost::variant<Line_3, Plane_3>>;

using Iterator = std::list<Triangle>::iterator;
using Primitive = CGAL::AABB_triangle_primitive<Kernel, Iterator>;
using AABB_triangle_traits = CGAL::AABB_traits<Kernel, Primitive>;
using Tree = CGAL::AABB_tree<AABB_triangle_traits>;

using SparseMatrix = Eigen::SparseMatrix<double>;


enum class GridCoordType {
  front_grid,
  front_cell,
  back_grid,
  back_cell,
  back_flat_x,
  back_flat_y
};
int gridCoord2Idx(int x, int y, GridCoordType type, int columns, int rows) {
  const int fullColumns = columns*2 - 1;
  const int backOffset = columns*rows + (columns-1)*(rows-1);
  switch (type) {
    case GridCoordType::front_grid:
      assert(x < columns && y < rows);
      return y*fullColumns + x;
    case GridCoordType::front_cell:
      assert(x < columns-1 && y < rows-1);
      return y*fullColumns + x + columns;
    case GridCoordType::back_grid:
      assert(x < columns && y < rows);
      return backOffset + y*fullColumns + x;
    case GridCoordType::back_cell:
      assert(x < columns-1 && y < rows-1);
      return backOffset + y*fullColumns + x + columns;
    case GridCoordType::back_flat_x:
      assert(x < columns && y < 2);
      if (y == 0) {
        return backOffset + x;
      } else {
        return backOffset + (columns+rows-2) + (columns-1)-x;
      }
    case GridCoordType::back_flat_y:
      assert(x < 2 && y < rows);
      if (x == 0) {
        return backOffset + ((columns-1 + (columns+rows-2) + (rows-1)-y) % (2 * (columns - 1) + 2 * (rows - 1)));
      } else {
        return backOffset + columns-1 + y;
      }
  }
  return 0;
}

static std::shared_ptr<AbstractNode> builtin_heightmap(const ModuleInstantiation *inst, Arguments arguments, const Children& children)
{
  if (!children.empty()) {
    LOG(message_group::Warning, inst->location(), arguments.documentRoot(),
        "module %1$s() does not support child modules", inst->name());
  }

  auto node = std::make_shared<HeightMapNode>(inst);

  Parameters parameters =
      Parameters::parse(std::move(arguments), inst->location(),
                        {"file", "size", "center", "convexity"}, {"invert", "doubleSided", "thickness", "pixelStep", "fadeTo", "fadeWidth", "shape", "cutout", "tiles"});

  std::string fileval = parameters["file"].isUndefined() ? "" : parameters["file"].toString();
  auto filename = lookup_file(fileval, inst->location().filePath().parent_path().string(), parameters.documentRoot());
  node->filename = filename;
  handle_dep(fs::path(filename).generic_string());

  if (parameters["size"].type() == Value::Type::VECTOR) {
    parameters["size"].getVec3(node->size.x(), node->size.y(), node->size.z());
  }

  if (parameters["center"].type() == Value::Type::BOOL) {
    node->center = parameters["center"].toBool();
  }

  if (parameters["convexity"].type() == Value::Type::NUMBER) {
    node->convexity = static_cast<int>(parameters["convexity"].toDouble());
  }

  if (parameters["doubleSided"].type() == Value::Type::BOOL) {
    node->doubleSided = parameters["doubleSided"].toBool();
  }

  if (parameters["thickness"].type() == Value::Type::NUMBER) {
    node->thickness = parameters["thickness"].toDouble();
  }

  if (parameters["pixelStep"].type() == Value::Type::NUMBER) {
    node->pixelStep = static_cast<int>(parameters["pixelStep"].toDouble());
    if (node->pixelStep < 1) {
      LOG(message_group::Warning, inst->location(), parameters.documentRoot(), "heightmap(..., pixelStep=%1$s) pixelStep parameter can not be less than 1, reset to 1.", parameters["pixelStep"].toEchoStringNoThrow());
      node->pixelStep = 1;
    }
  }

  if (parameters["fadeTo"].type() == Value::Type::VECTOR) {
    if (parameters["fadeTo"].toVector().size() != 4)
      LOG(message_group::Error, inst->location(), parameters.documentRoot(), "heightmap(..., fadeTo=%1$s) has to be either a vector with 4 elements or a number.", parameters["fadeTo"].toEchoStringNoThrow());
    auto& vec = parameters["fadeTo"].toVector();
    for(int i=0; i<4; ++i)
      node->fadeTo[i] = vec[i].toDouble();
  } else if (parameters["fadeTo"].type() == Value::Type::NUMBER) {
    double v = parameters["fadeTo"].toDouble();
    for(int i=0; i<4; ++i)
      node->fadeTo[i] = v;
  }

  if (parameters["fadeWidth"].type() == Value::Type::VECTOR) {
    if (parameters["fadeWidth"].toVector().size() != 4)
      LOG(message_group::Error, inst->location(), parameters.documentRoot(), "heightmap(..., fadeWidth=%1$s) has to be either a vector with 4 elements or a number.", parameters["fadeTo"].toEchoStringNoThrow());
    auto& vec = parameters["fadeWidth"].toVector();
    for(int i=0; i<4; ++i)
      node->fadeWidth[i] = vec[i].toDouble();
  } else if (parameters["fadeWidth"].type() == Value::Type::NUMBER) {
    double v = parameters["fadeWidth"].toDouble();
    for(int i=0; i<4; ++i)
      node->fadeWidth[i] = v;
  } else {
    for(int i=0; i<4; ++i)
      node->fadeWidth[i] = 0;
  }

  if (parameters["tiles"].type() == Value::Type::VECTOR) {
    parameters["tiles"].getVec2(node->tiles[0], node->tiles[1]);
  } else {
    node->tiles[0] = 1.0;
    node->tiles[1] = 1.0;
  }

  auto size = node->getDataSize(filename);
  auto width = size.first;
  auto height = size.second;

  const unsigned int columns = ceil(width / node->pixelStep) + 1;
  const unsigned int rows = ceil(height / node->pixelStep) + 1;

  const int m = (rows * columns) + (rows - 1) * (columns - 1);
  node->shapeGridPoints.resize(m);
  node->shapeGridNormals.resize(m);

  if (parameters["shape"].type() == Value::Type::FUNCTION) {
    const FunctionType& func = parameters["shape"].toFunction();

    std::stringstream exprStream;
    func.getExpr()->print(exprStream, "  ");
    node->shapeExpr = exprStream.str();

    if (func.getParameters()->size() != 2) {
      LOG(message_group::Error, inst->location(), parameters.documentRoot(), "heightmap(..., shape=%1$s) has to take exactly 2 arguments: u and v.", parameters["shape"].toEchoStringNoThrow());
    }

    ContextHandle<Context> c = Context::create<Context>(func.getContext());
    c->set_variable("u", 0.0);
    c->set_variable("v", 0.0);
    const Value ret = func.getExpr()->evaluate(*c);
    
    if (ret.type() != Value::Type::VECTOR || 
          (ret.toVector().size() != 3 && !(ret.toVector().size() == 2 && 
              ret.toVector()[0].type() == Value::Type::VECTOR && 
              ret.toVector()[0].toVector().size() == 3 && 
              ret.toVector()[1].type() == Value::Type::VECTOR && 
              ret.toVector()[1].toVector().size() == 3))) {
      LOG(message_group::Error, inst->location(), parameters.documentRoot(), "heightmap(..., shape=%1$s) has to return either a 3-vector (location) or 2 3-vectors (location+normal).", parameters["shape"].toEchoStringNoThrow());
    } else {

      if (ret.toVector().size() == 3) {
        for (unsigned int y=0; y<rows; ++y) {
          double v = y / static_cast<double>(rows-1);
          for (unsigned int x=0; x<columns; ++x) {
            double u = x / static_cast<double>(columns-1);

            c->set_variable("u", u);
            c->set_variable("v", v);

            Vector3d vec;
            func.getExpr()->evaluate(*c).getVec3(vec.x(), vec.y(), vec.z());
            node->shapeGridPoints[gridCoord2Idx(x, y, GridCoordType::front_grid, columns, rows)] = vec * node->size.x();
          }
        }
        for (unsigned int y=0; y<rows; ++y) {
          for (unsigned int x=0; x<columns; ++x) {
            const Vector3d& left  = node->shapeGridPoints[gridCoord2Idx(std::max(x, 1u)-1,        y, GridCoordType::front_grid, columns, rows)];
            const Vector3d& right = node->shapeGridPoints[gridCoord2Idx(std::min(x+1, columns-1), y, GridCoordType::front_grid, columns, rows)];
            const Vector3d xVec = (right - left).normalized();

            const Vector3d& top    = node->shapeGridPoints[gridCoord2Idx(x, std::max(y,   1u)-1,   GridCoordType::front_grid, columns, rows)];
            const Vector3d& bottom = node->shapeGridPoints[gridCoord2Idx(x, std::min(y+1, rows-1), GridCoordType::front_grid, columns, rows)];
            const Vector3d yVec = (bottom - top).normalized();
            node->shapeGridNormals[gridCoord2Idx(x, y, GridCoordType::front_grid, columns, rows)] = xVec.cross(yVec).normalized();
          }
        }

        for (unsigned int y=0; y<rows-1; ++y) {
          double v = (y+0.5) / static_cast<double>(rows-1);
          for (unsigned int x=0; x<columns-1; ++x) {
            double u = (x+0.5) / static_cast<double>(columns-1);

            c->set_variable("u", u);
            c->set_variable("v", v);

            Vector3d vec;
            func.getExpr()->evaluate(*c).getVec3(vec.x(), vec.y(), vec.z());
            node->shapeGridPoints[gridCoord2Idx(x, y, GridCoordType::front_cell, columns, rows)] = vec * node->size.x();
          }
        }
        for (unsigned int y=0; y<rows-1; ++y) {
          for (unsigned int x=0; x<columns-1; ++x) {
            const Vector3d& left  = node->shapeGridPoints[gridCoord2Idx(std::max(x, 1u)-1,        y, GridCoordType::front_cell, columns, rows)];
            const Vector3d& right = node->shapeGridPoints[gridCoord2Idx(std::min(x+1, columns-2), y, GridCoordType::front_cell, columns, rows)];
            const Vector3d xVec = (right - left).normalized();

            const Vector3d& top    = node->shapeGridPoints[gridCoord2Idx(x, std::max(y,   1u)-1,   GridCoordType::front_cell, columns, rows)];
            const Vector3d& bottom = node->shapeGridPoints[gridCoord2Idx(x, std::min(y+1, rows-2), GridCoordType::front_cell, columns, rows)];
            const Vector3d yVec = (bottom - top).normalized();
            node->shapeGridNormals[gridCoord2Idx(x, y, GridCoordType::front_cell, columns, rows)] = xVec.cross(yVec).normalized();
          }
        }
      } else {
        for (unsigned int y=0; y<rows; ++y) {
          double v = y / static_cast<double>(rows-1);
          for (unsigned int x=0; x<columns; ++x) {
            double u = x / static_cast<double>(columns-1);
            c->set_variable("u", u);
            c->set_variable("v", v);
            const auto val = func.getExpr()->evaluate(*c);
            const auto& vec = val.toVector();
            const auto& pos = vec[0].toVector();
            const auto& normal = vec[1].toVector();
            node->shapeGridPoints[gridCoord2Idx(x, y, GridCoordType::front_grid, columns, rows)] = Vector3d(pos[0].toDouble(), pos[1].toDouble(), pos[2].toDouble()) * node->size.x();
            node->shapeGridNormals[gridCoord2Idx(x, y, GridCoordType::front_grid, columns, rows)] = Vector3d(normal[0].toDouble(), normal[1].toDouble(), normal[2].toDouble()).normalized();
          }
        }
        for (unsigned int y=0; y<rows-1; ++y) {
          double v = (y+0.5) / static_cast<double>(rows-1);
          for (unsigned int x=0; x<columns-1; ++x) {
            double u = (x+0.5) / static_cast<double>(columns-1);
            c->set_variable("u", u);
            c->set_variable("v", v);
            const auto val = func.getExpr()->evaluate(*c);
            const auto& vec = val.toVector();
            const auto& pos = vec[0].toVector();
            const auto& normal = vec[1].toVector();
            node->shapeGridPoints[gridCoord2Idx(x, y, GridCoordType::front_cell, columns, rows)] = Vector3d(pos[0].toDouble(), pos[1].toDouble(), pos[2].toDouble()) * node->size.x();
            node->shapeGridNormals[gridCoord2Idx(x, y, GridCoordType::front_cell, columns, rows)] = Vector3d(normal[0].toDouble(), normal[1].toDouble(), normal[2].toDouble()).normalized();
          }
        }
      }
    }
  } else {
    for (unsigned int y=0; y<rows; ++y) {
      const double v = y / static_cast<double>(rows-1);
      for (unsigned int x=0; x<columns; ++x) {
        const double u = x / static_cast<double>(columns-1);
        node->shapeGridPoints[gridCoord2Idx(x, y, GridCoordType::front_grid, columns, rows)] = Vector3d(u, v, 0.5).cwiseProduct(node->size);
        node->shapeGridNormals[gridCoord2Idx(x, y, GridCoordType::front_grid, columns, rows)] = Vector3d(0, 0, 1);
      }
    }
    for (unsigned int y=0; y<rows-1; ++y) {
      const double v = (y+0.5) / static_cast<double>(rows-1);
      for (unsigned int x=0; x<columns-1; ++x) {
        const double u = (x+0.5) / static_cast<double>(columns-1);
        node->shapeGridPoints[gridCoord2Idx(x, y, GridCoordType::front_cell, columns, rows)] = Vector3d(u, v, 0.5).cwiseProduct(node->size);
        node->shapeGridNormals[gridCoord2Idx(x, y, GridCoordType::front_cell, columns, rows)] = Vector3d(0, 0, 1);
      }
    }
  }


  return node;
}

void HeightMapNode::convert_image(map_data_t& data, std::vector<uint8_t>& img, unsigned int width, unsigned int height) const
{
  data.width = width;
  data.height = height;
  data.resize( (size_t)width * height);
  for (unsigned int y = 0; y < height; ++y) {
    for (unsigned int x = 0; x < width; ++x) {
      const long idx = 4l * (y * width + x);
      const double pixel = 0.2126 * img[idx] + 0.7152 * img[idx + 1] + 0.0722 * img[idx + 2];
      const double z = (pixel / 255.0) - 0.5;
      data[ x + (width * (height - 1 - y)) ] = z;
    }
  }
}

bool HeightMapNode::is_png(std::vector<uint8_t>& png) const
{
  const size_t pngHeaderLength = 8;
  const uint8_t pngHeader[pngHeaderLength] = {0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a};
  return (png.size() >= pngHeaderLength &&
          std::memcmp(png.data(), pngHeader, pngHeaderLength) == 0);
}

map_data_t HeightMapNode::read_png_or_dat(std::string filename) const
{
  map_data_t data;
  std::vector<uint8_t> png;
  int ret_val = 0;
  try{
    ret_val = lodepng::load_file(png, filename);
  } catch (std::bad_alloc& ba) {
    LOG(message_group::Warning, "bad_alloc caught for '%1$s'.", ba.what());
    return data;
  }

  if (ret_val == 78) {
    LOG(message_group::Warning, "The file '%1$s' couldn't be opened.", filename);
    return data;
  }

  if (!is_png(png)) {
    png.clear();
    return read_dat(filename);
  }

  unsigned int width, height;
  std::vector<uint8_t> img;
  auto error = lodepng::decode(img, width, height, png);
  if (error) {
    LOG(message_group::Warning, "Can't read PNG image '%1$s'", filename);
    data.clear();
    return data;
  }

  convert_image(data, img, width, height);

  return data;
}

map_data_t HeightMapNode::read_dat(std::string filename) const
{
  map_data_t data;
  std::ifstream stream(filename.c_str());

  if (!stream.good()) {
    LOG(message_group::Warning, "Can't open DAT file '%1$s'.", filename);
    return data;
  }

  int lines = 0, columns = 0;

  using tokenizer = boost::tokenizer<boost::char_separator<char>>;
  boost::char_separator<char> sep(" \t");

  // We use an unordered map because the data file may not be rectangular,
  // and we may need to fill in some bits.
  using unordered_image_data_t = std::unordered_map<std::pair<int, int>, double, boost::hash<std::pair<int, int>>>;
  unordered_image_data_t unordered_data;

  while (!stream.eof()) {
    std::string line;
    while (!stream.eof() && (line.size() == 0 || line[0] == '#')) {
      std::getline(stream, line);
      boost::trim(line);
    }
    if (line.size() == 0 && stream.eof()) break;

    int col = 0;
    tokenizer tokens(line, sep);
    try {
      for (const auto& token : tokens) {
        auto v = boost::lexical_cast<double>(token);
        unordered_data[ std::make_pair(lines, col++) ] = v;
        if (col > columns) columns = col;
      }
    } catch (const boost::bad_lexical_cast& blc) {
      if (!stream.eof()) {
        LOG(message_group::Warning, "Illegal value in '%1$s': %2$s", filename, blc.what());
      }
      return data;
    }

    lines++;
  }

  data.width = columns;
  data.height = lines;

  // Now convert the unordered, possibly non-rectangular data into a well ordered vector
  // for faster access and reduced memory usage.
  data.resize( (size_t)lines * columns);
  for (int i = 0; i < lines; ++i)
    for (int j = 0; j < columns; ++j)
      data[ i * columns + j ] = unordered_data[std::make_pair(i, j)];

  return data;
}

std::pair<unsigned int, unsigned int> HeightMapNode::getDataSize(std::string filename) const
{
  auto data = read_png_or_dat(std::move(filename));
  return std::make_pair(data.width, data.height);
}

std::pair<Point, Point> getClosestPoints(const Line_3& l1, const Line_3& l2) {
  Plane_3 p1(l1.point(), CGAL::cross_product(CGAL::cross_product(l1.direction().to_vector(), l2.direction().to_vector()), l1.direction().to_vector()));
  Plane_3 p2(l2.point(), CGAL::cross_product(CGAL::cross_product(l1.direction().to_vector(), l2.direction().to_vector()), l2.direction().to_vector()));
  return std::make_pair(boost::get<Point>(*CGAL::intersection(p1, l2)), boost::get<Point>(*CGAL::intersection(p2, l1)));
}

boost::optional<Point> intersectWithTolerance(const Line_3& l1, const Line_3& l2, double tolerance = 0.001) {
  auto points = getClosestPoints(l1,l2);
  if(CGAL::squared_distance(points.first, points.second) < tolerance * tolerance)
    return CGAL::midpoint(points.first, points.second);
  return {};
}

std::vector<Point> cutPolygonAtPlane(std::vector<Point> polygon, const Plane_3& plane) {
  std::vector<bool> isInFront(polygon.size());
  std::transform(polygon.begin(), polygon.end(), isInFront.begin(), [&plane](const Point& p) { return plane.has_on_positive_side(p); });

  int numPointsBehind = std::count_if(isInFront.begin(), isInFront.end(), [&plane](bool b) {
    return !b;
  });

  if (numPointsBehind == 0)
    return polygon;

  if (numPointsBehind == polygon.size()) {
    return polygon;
  }

  for(int i=1; i<polygon.size(); ++i) {
    if (isInFront[i-1] != isInFront[i]) {
      auto intersection = CGAL::intersection(plane, Segment(polygon[i-1], polygon[i]));
      if (intersection && intersection->type() == typeid(Point)) {
        int indexToChange = isInFront[i-1] ? i : i-1;
        polygon[indexToChange] = boost::get<Point>(*intersection);
        isInFront[indexToChange] = true;
      }
    } else if(!isInFront[i]) {
      polygon.erase(polygon.begin()+i);
      --i;
    }
  }

  return polygon;
}

Point toPoint(const Vector3d& p) {
  return {p.x(), p.y(), p.z()};
}

Vector3d toVector3d(const Point& p) {
  return {p.x(), p.y(), p.z()};
}

std::vector<Vector3d> toPoly(const Triangle& t) {
  return {toVector3d(t.vertex(0)), toVector3d(t.vertex(1)), toVector3d(t.vertex(2))};
}

Vector3d getBottomSphereZLineIntersection(const Vector3d& center, double r, Vector3d lineOrigin) {
  const double distX = center.x() - lineOrigin.x();
  const double distY = center.y() - lineOrigin.y();
  const double squaredDistance = distX * distX + distY * distY;
  if (squaredDistance >= r * r)
    return lineOrigin;

  const double squaredHeight = (r*r) - squaredDistance;
  const double height = std::sqrt(squaredHeight);

  lineOrigin.z() = std::min(center.z() - height, lineOrigin.z());
  return lineOrigin;
}

Vector3d getLastSphereLineIntersection(const Vector3d& center, double r, Vector3d lineOrigin, const Vector3d& direction) {
  const Vector3d dist = center - lineOrigin;
  const Vector3d zComponent = direction.dot(dist) * direction;
  const double squaredDistance = (dist - zComponent).squaredNorm();
  if (squaredDistance >= r * r)
    return lineOrigin;

  const double squaredHeight = (r*r) - squaredDistance;
  const double height = std::sqrt(squaredHeight);

  Vector3d bottomIntersection = lineOrigin + zComponent + height * direction;
  if (direction.dot(bottomIntersection - lineOrigin) > 0)
    return bottomIntersection;
  else
    return lineOrigin;
}

Vector3d getTriangleZLineIntersection(const Triangle& t, Vector3d p) {
  const Ray r(toPoint(p), Vector_3(0,0,-1));
  auto intersect = CGAL::intersection(r, t);

  if (intersect && intersect->type() == typeid(Point)) {
    p = toVector3d(boost::get<Point>(*intersect));
  }
  return p;
}

Vector3d getTriangleLineIntersection(const Triangle& t, Vector3d p, const Vector3d& direction) {
  const Ray r(toPoint(p), Vector_3(direction.x(),direction.y(),direction.z()));
  auto intersect = CGAL::intersection(r, t);

  if (intersect && intersect->type() == typeid(Point)) {
    p = toVector3d(boost::get<Point>(*intersect));
  }
  return p;
}

double getSweepingSphereZLineIntersectionZAtStart(const Vector3d& start, const Vector3d& end, double thickness) {
  const double lenX = end.x() - start.x();
  const double lenY = end.y() - start.y();
  const double len = sqrt(lenX * lenX + lenY * lenY);
  const double up = end.z() - start.z();
  Vector_2 n(up,len);
  const double lenN = sqrt(up * up + len * len);
  n /= lenN;
  n *= thickness;
  const Line_2 line(Point_2(-n.x(),n.y()+start.z()), Vector_2(len, up));
  return line.y_at_x(0.0);
}

Vector3d getPillLineIntersectionAtPillStart(const Vector3d& start, const Vector3d& end, double r, const Vector3d& direction) {
  const Vector3d line = end - start;
  const double up = direction.dot(line);

  if (up <= 0.0)
    return start + r * direction;

  const Vector3d vertComponent = up * direction;
  const Vector3d horComponent = line - vertComponent;
  const double len = horComponent.norm();

  Vector_2 n(up,len);
  const double lenN = sqrt(up * up + len * len);
  n /= lenN;
  n *= r;
  const Line_2 line2(Point_2(-n.x(),n.y()), Vector_2(len, up));
  return start + direction * line2.y_at_x(0.0);
}

double getZAt(int row, int column, double xDataStep, double yDataStep, const map_data_t& data)
{
  return data[round(column * xDataStep) + round(row * yDataStep) * data.width];
}

// FIXME: Look for faster way to generate PolySet directly
std::unique_ptr<const Geometry> HeightMapNode::createGeometry() const
{
  if (filename.empty())
    return {};
  auto data = read_png_or_dat(filename);

  const int columns = ceil(data.width / pixelStep) + 1;
  const int rows = ceil(data.height / pixelStep) + 1;

  const double width = size.x();
  const double height = size.y();

  const double xStep = width / static_cast<double>(columns-1);
  const double yStep = height / static_cast<double>(rows-1);

  const double xDataStep = (data.width - 1.0) / (columns - 1.0);
  const double yDataStep = (data.height - 1.0) / (rows - 1.0);

  const double ox = center ? -width / 2.0 : 0;
  const double oy = center ? -height / 2.0 : 0;
  const Vector3d centerOffset(ox, oy, 0);

  const int m = (rows * columns) + (rows - 1) * (columns - 1);

  
  for (unsigned int y=0; y<std::ceil(fadeWidth[0]*data.height); ++y) {
    const double a = sqrt(y / std::ceil(fadeWidth[0]*data.height));
    for (unsigned int x=0; x<data.width; ++x) {
      data[x + y * data.width] = a * data[x + y * data.width] + (1-a) * fadeTo[0] * size.z();
    }
  }
  for (unsigned int y=0; y<std::ceil(fadeWidth[1]*data.height); ++y) {
    const double a = sqrt(y / std::ceil(fadeWidth[1]*data.height));
    for (unsigned int x=0; x<data.width; ++x) {
      data[x + (data.height - y - 1) * data.width] = a * data[x + (data.height - y - 1) * data.width] + (1-a) * fadeTo[1] * size.z();
    }
  }
  
  double relWidth = std::ceil(fadeWidth[2]*data.width);
  for (unsigned int x=0; x<relWidth; ++x) {
    const double a = sqrt(x / relWidth);
    for (unsigned int y=0; y<data.height; ++y) {
      data[x + y * data.width] = a * data[x + y * data.width] + (1-a) * fadeTo[2] * size.z();
    }
  }
  relWidth = std::ceil(fadeWidth[3]*data.width);
  for (unsigned int x=0; x<relWidth; ++x) {
    const double a = sqrt(x / relWidth);
    for (unsigned int y=0; y<data.height; ++y) {
      data[data.width - x - 1 + y * data.width] = a * data[data.width - x - 1 + y * data.width] + (1-a) * fadeTo[3] * size.z();
    }
  }

  // reserve the polygon vector size so we don't have to reallocate as often
  int numIndices = ((rows - 1) * (columns - 1) * 8); // heightmap (on top and bottom)
  numIndices += ((rows - 1) * 2 + (columns - 1) * 2); // sides

  int numVertices = 2 * ((rows * columns) + (rows - 1) * (columns - 1));
                      
  auto polyset = std::make_unique<PolySet>(3);
  polyset->setConvexity(convexity);
  polyset->vertices.resize(numVertices);
  polyset->indices.reserve(numIndices);

  for (unsigned int i = 0; i < rows; ++i) {
    for (unsigned int j = 0; j < columns; ++j) {
      const double z = getZAt(i, j, xDataStep, yDataStep, data);
      const int idx = gridCoord2Idx(j, i, GridCoordType::front_grid, columns, rows);
      const auto& pos = shapeGridPoints[idx] + centerOffset;
      const auto& normal = shapeGridNormals[idx];
      polyset->vertices[idx] = pos + normal * z * size.z();
    }
  }
  for (unsigned int i = 0; i < rows-1; ++i) {
    for (unsigned int j = 0; j < columns-1; ++j) {
      double z = getZAt(i, j, xDataStep, yDataStep, data);
      z += getZAt(i+0, j+1, xDataStep, yDataStep, data);
      z += getZAt(i+1, j+0, xDataStep, yDataStep, data);
      z += getZAt(i+1, j+1, xDataStep, yDataStep, data);
      z /= 4.0;
      const int idx = gridCoord2Idx(j, i, GridCoordType::front_cell, columns, rows);
      const auto& pos = shapeGridPoints[idx] + centerOffset;
      const auto& normal = shapeGridNormals[idx];
      polyset->vertices[idx] = (pos + normal * z * size.z());
    }
  }

  // the bulk of the heightmap
  for (unsigned int i = 0; i < rows-1; ++i) {
      const int top = i;
      const int bottom = i + 1;

      for (unsigned int j = 0; j < columns-1; ++j) {
        const int left = j;
        const int right = j + 1;

        const int topLeft     = gridCoord2Idx(left,  top,    GridCoordType::front_grid, columns, rows);
        const int topRight    = gridCoord2Idx(right, top,    GridCoordType::front_grid, columns, rows);
        const int bottomLeft  = gridCoord2Idx(left,  bottom, GridCoordType::front_grid, columns, rows);
        const int bottomRight = gridCoord2Idx(right, bottom, GridCoordType::front_grid, columns, rows);
        const int center      = gridCoord2Idx(left,  top,    GridCoordType::front_cell, columns, rows);

        polyset->indices.push_back(IndexedFace({topLeft, topRight, center}));
        polyset->indices.push_back(IndexedFace({bottomLeft, topLeft, center}));
        polyset->indices.push_back(IndexedFace({bottomRight, bottomLeft, center}));
        polyset->indices.push_back(IndexedFace({topRight, bottomRight, center}));
    }
  }

  if (doubleSided) {
    for (int y = 0; y < rows; ++y) {
      for (int x = 0; x < columns; ++x) {
        const int frontIdx = gridCoord2Idx(x, y, GridCoordType::front_grid, columns, rows);
        const auto& normal = shapeGridNormals[frontIdx];
        polyset->vertices[gridCoord2Idx(x, y, GridCoordType::back_grid, columns, rows)]  = polyset->vertices[frontIdx] - thickness * normal;
      }
    }
    for (int y = 0; y < rows-1; ++y) {
      for (int x = 0; x < columns-1; ++x) {
        const int frontIdx = gridCoord2Idx(x, y, GridCoordType::front_cell, columns, rows);
        const auto& normal = shapeGridNormals[frontIdx];
        polyset->vertices[gridCoord2Idx(x, y, GridCoordType::back_cell, columns, rows)]  = polyset->vertices[frontIdx] - thickness * normal;
      }
    }

    // std::vector<Eigen::Triplet<double>> tripletList;

    // for (int y1 = 0; y1 < rows; ++y1) {
    //   for (int x1 = 0; x1 < columns; ++x1) {
    //     const int idx1 = gridCoord2Idx(x1, y1, GridCoordType::front_grid, columns, rows);
    //     const Vector3d& currentPoint = polyset->vertices[idx1];

    //     for (int y2 = 0; y2 < rows; ++y2) {
    //       for (int x2 = 0; x2 < columns; ++x2) {
    //         const int idx2 = gridCoord2Idx(x2, y2, GridCoordType::front_grid, columns, rows);
    //         const double dist = (currentPoint - polyset->vertices[idx2+m]).norm();
    //         if (dist < thickness * 2)
    //           tripletList.emplace_back(idx1, idx2, dist);
    //       }
    //     }
        
    //     for (int y2 = 0; y2 < rows-1; ++y2) {
    //       for (int x2 = 0; x2 < columns-1; ++x2) {
    //         const int idx2 = gridCoord2Idx(x2, y2, GridCoordType::front_cell, columns, rows);
    //         const double dist = (currentPoint - polyset->vertices[idx2+m]).norm();
    //         if (dist < thickness * 2)
    //           tripletList.emplace_back(idx1, idx2, dist);
    //       }
    //     }
    //   }
    // }

    // for (int y1 = 0; y1 < rows-1; ++y1) {
    //   for (int x1 = 0; x1 < columns-1; ++x1) {
    //     const int idx1 = gridCoord2Idx(x1, y1, GridCoordType::front_cell, columns, rows);
    //     const Vector3d& currentPoint = polyset->vertices[idx1];

    //     for (int y2 = 0; y2 < rows; ++y2) {
    //       for (int x2 = 0; x2 < columns; ++x2) {
    //         const int idx2 = gridCoord2Idx(x2, y2, GridCoordType::front_grid, columns, rows);
    //         const double dist = (currentPoint - polyset->vertices[idx2+m]).norm();
    //         if (dist < thickness * 2)
    //           tripletList.emplace_back(idx1, idx2, dist);
    //       }
    //     }
        
    //     for (int y2 = 0; y2 < rows-1; ++y2) {
    //       for (int x2 = 0; x2 < columns-1; ++x2) {
    //         const int idx2 = gridCoord2Idx(x2, y2, GridCoordType::front_cell, columns, rows);
    //         const double dist = (currentPoint - polyset->vertices[idx2+m]).norm();
    //         if (dist < thickness * 2)
    //           tripletList.emplace_back(idx1, idx2, dist);
    //       }
    //     }
    //   }
    // }

    // create sparse matrix with the distance of the front vertices (rows) to the back vertices (columns)
    // SparseMatrix influenceMatrix(m, m);
    // influenceMatrix.setFromTriplets(tripletList.begin(), tripletList.end());

    const int rangeOfInfluenceX = std::ceil(thickness / xStep);
    const int rangeOfInfluenceY = std::ceil(thickness / yStep);

    std::list<Triangle> triangles;
    for(const auto& t : polyset->indices) {
      if (t.size() == 3) {
        const Point p1 = toPoint(polyset->vertices[t[0]]);
        const Point p2 = toPoint(polyset->vertices[t[1]]);
        const Point p3 = toPoint(polyset->vertices[t[2]]);
        triangles.emplace_back(p1, p2, p3);
      }
    }

    Tree tree(triangles.begin(), triangles.end());
    tree.accelerate_distance_queries();

    for (int k=0; k<m; ++k) {
      Vector3d& currentPoint = polyset->vertices[k+m];
      const Vector3d& normal = -shapeGridNormals[k];
      auto closest = tree.closest_point_and_primitive(toPoint(currentPoint));
      double sqrDist = (toVector3d(closest.first)-currentPoint).squaredNorm();
      while (sqrDist < thickness * thickness) {
        currentPoint = getLastSphereLineIntersection(toVector3d(closest.first), thickness*1.01, currentPoint, normal);
        for (int i=0; i<3; ++i) {
          currentPoint = getLastSphereLineIntersection(toVector3d(closest.second->vertex(i)), thickness*1.01, currentPoint, normal);
        }
        if (rangeOfInfluenceX < 5) {
          const Point& p1 = closest.second->vertex(0);
          const Point& p2 = closest.second->vertex(1);
          const Point& p3 = closest.second->vertex(2);
          const Vector_3 n = CGAL::unit_normal(p1, p2, p3);
          const auto o = thickness * 1.01 * n;
          currentPoint = getTriangleLineIntersection(Triangle(p1 + o, p2 + o, p3 + o), currentPoint, normal);
        }

        closest = tree.closest_point_and_primitive(toPoint(currentPoint), closest);
        sqrDist = (toVector3d(closest.first)-currentPoint).squaredNorm();
      }
    }

    // // per-point operations
    // if (rangeOfInfluenceX > 1) {
    //   for (int k=0; k<influenceMatrix.outerSize(); ++k) { // iterate over columns aka back vertices
    //     const auto& normal = -shapeGridNormals[k];
    //     Vector3d& currentPoint = polyset->vertices[k+m];

    //     // InnerIterator iterates over the rows in one column, i.e. the front vertices close enough to the kth back vertex
    //     for (SparseMatrix::InnerIterator it(influenceMatrix,k); it; ++it) { 
    //       const Vector3d& otherPoint = polyset->vertices[it.row()];
    //       currentPoint = getLastSphereLineIntersection(otherPoint, thickness, currentPoint, normal);
    //     }
    //   }
    // }

    // if (rangeOfInfluenceX < 5) {
    //   std::vector<Triangle> triangles;
    //   triangles.reserve((rows-1) * (columns-1) * 4);
    //   for(const auto& t : polyset->indices) {
    //     if (t.size() == 3) {
    //       const Point p1 = toPoint(polyset->vertices[t[0]]);
    //       const Point p2 = toPoint(polyset->vertices[t[1]]);
    //       const Point p3 = toPoint(polyset->vertices[t[2]]);
          // const Vector_3 n = CGAL::unit_normal(p1, p2, p3);
          // const auto o = thickness * n;
          // triangles.emplace_back(p1 + o, p2 + o, p3 + o);
    //     }
    //   }
      
    //   const int fullColumns = columns*2 - 1;
    //   for (int k=0; k<influenceMatrix.outerSize(); ++k) { // iterate over columns aka back vertices
    //     const auto& normal = -shapeGridNormals[k];
    //     Vector3d& currentPoint = polyset->vertices[k+m];

    //     // InnerIterator iterates over the rows in one column, i.e. the front vertices close enough to the kth back vertex
    //     for (SparseMatrix::InnerIterator it(influenceMatrix,k); it; ++it) {
    //       if (it.row() % fullColumns > columns) {
    //         const int cellX = it.row() % fullColumns - columns;
    //         const int cellY = it.row() / fullColumns;
    //         const int cellIdx = cellY*(columns-1) + cellX;
    //         for (int t = 0; t < 4; ++t) {
    //           const auto& tri = triangles[cellIdx * 4 + t];
    //           currentPoint = getTriangleLineIntersection(tri, currentPoint, normal);
    //         }
    //       }
    //     }
    //   }
    // }

    // if (rangeOfInfluenceX == 1) {
    //   const auto gridNeighbors = std::array<std::pair<int, int>,4>{std::make_pair(0,1), {1,0}, {0,-1}, {-1,0}};
    //   const auto cellNeighbors = std::array<std::pair<int, int>,4>{std::make_pair(0,0), {-1,0}, {-1,-1}, {0,-1}};
    //   for (int i = 0; i < rows; ++i) {
    //     for (int j = 0; j < columns; ++j) {
    //       const int frontIdx = gridCoord2Idx(j, i, GridCoordType::front_grid, columns, rows);
    //       const auto& normal = shapeGridNormals[frontIdx];
    //       const Vector3d& frontPoint = polyset->vertices[frontIdx];
    //       Vector3d& currentPoint = polyset->vertices[gridCoord2Idx(j, i, GridCoordType::back_grid, columns, rows)];

    //       for (const auto& coord: gridNeighbors) {
    //         const int x = j + coord.first;
    //         if (x<0 || x>=columns) continue;
    //         const int y = i + coord.second;
    //         if (y<0 || y>=rows) continue;

    //         const Vector3d& otherPoint = polyset->vertices[gridCoord2Idx(x, y, GridCoordType::front_grid, columns, rows)];
    //         const Vector3d intersection = getPillLineIntersectionAtPillStart(frontPoint, otherPoint, thickness, -normal);
    //         if (normal.dot(currentPoint - intersection) > 0)
    //           currentPoint = intersection;
    //       }
    //       for (const auto& coord: cellNeighbors) {
    //         const int x = j + coord.first;
    //         if (x<0 || x>=columns-1) continue;
    //         const int y = i + coord.second;
    //         if (y<0 || y>=rows-1) continue;

    //         const Vector3d& otherPoint = polyset->vertices[gridCoord2Idx(x, y, GridCoordType::front_cell, columns, rows)];
    //         const Vector3d intersection = getPillLineIntersectionAtPillStart(frontPoint, otherPoint, thickness, -normal);
    //         if (normal.dot(currentPoint - intersection) > 0)
    //           currentPoint = intersection;
    //       }
    //     }
    //   }
    //   const auto cellGridNeighbors = std::array<std::pair<int, int>,4>{std::make_pair(0,0), {1,0}, {1,1}, {0,1}};
    //   for (int i = 0; i < rows-1; ++i) {
    //     for (int j = 0; j < columns-1; ++j) {
    //       const int frontIdx = gridCoord2Idx(j, i, GridCoordType::front_cell, columns, rows);
    //       const auto& normal = shapeGridNormals[frontIdx];
    //       const Vector3d& frontPoint = polyset->vertices[frontIdx];
    //       Vector3d& currentPoint = polyset->vertices[gridCoord2Idx(j, i, GridCoordType::back_cell, columns, rows)];

    //       for (const auto& coord: cellGridNeighbors) {
    //         const int x = j + coord.first;
    //         const int y = i + coord.second;

    //         const Vector3d& otherPoint = polyset->vertices[gridCoord2Idx(x, y, GridCoordType::front_grid, columns, rows)];
    //         const Vector3d intersection = getPillLineIntersectionAtPillStart(frontPoint, otherPoint, thickness, -normal);
    //         if (normal.dot(currentPoint - intersection) > 0)
    //           currentPoint = intersection;
    //       }
    //     }
    //   }
    // }
  }
  else if (columns > 1 && rows > 1) {
    for (unsigned int i = 0; i < rows; ++i) {
      for (unsigned int j = 0; j < columns; ++j) {
        const int idx = gridCoord2Idx(j, i, GridCoordType::front_grid, columns, rows);
        const auto& pos = shapeGridPoints[idx] + centerOffset;
        const auto& normal = shapeGridNormals[idx];
        polyset->vertices[idx+m] = pos - normal * (size.z()/2 + thickness);
      }
    }
    for (unsigned int i = 0; i < rows-1; ++i) {
      for (unsigned int j = 0; j < columns-1; ++j) {
        const int idx = gridCoord2Idx(j, i, GridCoordType::front_cell, columns, rows);
        const auto& pos = shapeGridPoints[idx] + centerOffset;
        const auto& normal = shapeGridNormals[idx];
        polyset->vertices[idx+m] = pos - normal * (size.z()/2 + thickness);
      }
    }
  }

  for (unsigned int i = 1; i < rows; ++i) {
    const int top = (i - 1);
    const int bottom = i;

    for (unsigned int j = 1; j < columns; ++j) {
      const int left = (j - 1);
      const int right = j;

      const int topLeft     = gridCoord2Idx(left,  top,    GridCoordType::back_grid, columns, rows);
      const int topRight    = gridCoord2Idx(right, top,    GridCoordType::back_grid, columns, rows);
      const int bottomLeft  = gridCoord2Idx(left,  bottom, GridCoordType::back_grid, columns, rows);
      const int bottomRight = gridCoord2Idx(right, bottom, GridCoordType::back_grid, columns, rows);
      const int center      = gridCoord2Idx(left,  top,    GridCoordType::back_cell, columns, rows);

      polyset->indices.push_back(IndexedFace({topRight, topLeft, center}));
      polyset->indices.push_back(IndexedFace({topLeft, bottomLeft, center}));
      polyset->indices.push_back(IndexedFace({bottomLeft, bottomRight, center}));
      polyset->indices.push_back(IndexedFace({bottomRight, topRight, center}));
    }
  }

  const GridCoordType backType_x = GridCoordType::back_grid;
  const GridCoordType backType_y = GridCoordType::back_grid;
  const int back_maxX = columns-1;
  const int back_maxY = rows-1;

  // edges along Y
  for (int i = 0; i < rows-1; ++i) {
    polyset->indices.push_back(IndexedFace({
      gridCoord2Idx(0, i, GridCoordType::front_grid, columns, rows), 
      gridCoord2Idx(0, i+1, GridCoordType::front_grid, columns, rows), 
      gridCoord2Idx(0, i+1, backType_y, columns, rows), 
      gridCoord2Idx(0, i, backType_y, columns, rows)}));

    polyset->indices.push_back(IndexedFace({
      gridCoord2Idx(columns - 1, i+1, GridCoordType::front_grid, columns, rows), 
      gridCoord2Idx(columns - 1, i,   GridCoordType::front_grid, columns, rows), 
      gridCoord2Idx(back_maxX,   i, backType_y, columns, rows), 
      gridCoord2Idx(back_maxX,   i+1, backType_y, columns, rows)}));
  }

  // edges along X
  for (int i = 0; i < columns-1; ++i) {
    polyset->indices.push_back(IndexedFace({
      gridCoord2Idx(i+1, 0, GridCoordType::front_grid, columns, rows), 
      gridCoord2Idx(i,   0, GridCoordType::front_grid, columns, rows),
      gridCoord2Idx(i,   0, backType_x, columns, rows), 
      gridCoord2Idx(i+1, 0, backType_x, columns, rows)}));

    polyset->indices.push_back(IndexedFace({
      gridCoord2Idx(i,   rows - 1, GridCoordType::front_grid, columns, rows), 
      gridCoord2Idx(i+1, rows - 1, GridCoordType::front_grid, columns, rows), 
      gridCoord2Idx(i+1, back_maxY, backType_x, columns, rows), 
      gridCoord2Idx(i,   back_maxY, backType_x, columns, rows)}));
  }

  return polyset;
}

std::string HeightMapNode::toString() const
{
  std::ostringstream stream;
  fs::path path{static_cast<std::string>(this->filename)}; // gcc-4.6

  stream << this->name() << "(file = " << this->filename
         << ", size = " << (this->size)
         << ", center = " << (this->center ? "true" : "false")
         << ", doubleSided = " << (this->doubleSided ? "true" : "false")
         << ", thickness = " << this->thickness
         << ", pixelStep = " << this->pixelStep
         << ", fadeTo = [" << this->fadeTo[0] << ", " << this->fadeTo[1] << ", " << this->fadeTo[2] << ", " << this->fadeTo[3] << "]"
         << ", fadeWidth = [" << this->fadeWidth[0] << ", " << this->fadeWidth[1] << ", " << this->fadeWidth[2] << ", " << this->fadeWidth[3] << "]"
         << ", shape = " << shapeExpr
         << ", "
            "timestamp = "
         << (fs::exists(path) ? fs::last_write_time(path) : 0) << ")";

  return stream.str();
}

void register_builtin_heightmap()
{
  Builtins::init("heightmap", new BuiltinModule(builtin_heightmap),
                 {
                     "heightmap(string, size = [1,1,1], center = false, invert = false, convexity = number, doubleSided = false, thickness = 1, pixelStep = 1, fadeTo = 0, fadeWidth = 0, shape = 0)",
                 });
}
