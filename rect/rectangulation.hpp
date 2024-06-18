/*
 * Copyright (c) 2021 Arturo Merino and Torsten Muetze
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include "edge.hpp"
#include "rectangle.hpp"
#include "vertex.hpp"
#include "wall.hpp"
#include <array>
#include <vector>

enum class RectangulationType { generic, baligned, diagonal };
enum class RectangulationDirection { left, right, None };
enum class RectangulationPattern {
  wmill_clockwise,
  wmill_counterclockwise,
  brick_leftright,
  brick_rightleft,
  brick_topbottom,
  brick_bottomtop,
  H_vertical,
  H_horizontal
};

class Rectangulation {
public:
private:
  int n_;
  RectangulationType type_;
  std::vector<RectangulationPattern> patterns_;
  std::vector<RectangulationDirection> o_;
  std::vector<int> s_;
  void set_all_vertical();

public:
  std::vector<Vertex> vertices_;
  std::vector<Wall> walls_;
  std::vector<Edge> edges_;
  std::vector<Rectangle> rectangles_;
  Rectangulation(int, RectangulationType, std::vector<RectangulationPattern> &);
  void init(std::vector<Vertex> &, std::vector<Wall> &, std::vector<Edge> &,
            std::vector<Rectangle> &);
  void print_data();
  void print_coordinates();
  void print_coordinates_generic();
  void print_coordinates_diagonal();
  void DFS_BL(int, int &, std::vector<int> &, std::vector<int> &);
  void DFS_TR(int, int &, std::vector<int> &, std::vector<int> &);
  bool next();
  bool is_bottom_based(int);
  bool is_right_based(int);
  void remHead(int);
  void remTail(int);
  void insBefore(int, int, int);
  void insAfter(int, int, int);
  void Wjump_hor(int, RectangulationDirection, int);
  void Wjump_ver(int, RectangulationDirection, int);
  void Sjump(int, RectangulationDirection, int);
  void Tjump_hor(int, RectangulationDirection, int);
  void Tjump_ver(int, RectangulationDirection, int);
  void next_generic(int, RectangulationDirection);
  void next_diagonal(int, RectangulationDirection);
  void next_baligned(int, RectangulationDirection);
  void lock(int, EdgeDir);
  void unlock(int, RectangulationDirection);
  bool contains_pattern(int);
  bool contains_brick_leftright(int);
  bool contains_brick_rightleft(int);
  bool contains_brick_bottomtop(int);
  bool contains_brick_topbottom(int);
  bool contains_wmill_clockwise(int);
  bool contains_wmill_counterclockwise(int);
  bool contains_H_vertical(int);
  bool contains_H_horizontal(int);
};
