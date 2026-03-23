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

/**
 * @brief Type of rectangulation.
 *
 * - generic: Any rectangulation with n rectangles
 * - diagonal: Rectangulations where all vertices lie on the main diagonal
 * - baligned: Block-aligned rectangulations (subset of diagonal)
 */
enum class RectangulationType { generic, baligned, diagonal };

/**
 * @brief Direction for jump operations in rectangulation traversal.
 *
 * Used to specify the direction of S-jump, T-jump, and W-jump operations.
 */
enum class RectangulationDirection { left, right, None };

/**
 * @brief Pattern types that can be forbidden in rectangulations.
 *
 * - wmill_clockwise: Clockwise windmill pattern
 * - wmill_counterclockwise: Counter-clockwise windmill pattern
 * - brick_leftright: Left-right brick pattern
 * - brick_rightleft: Right-left brick pattern
 * - brick_topbottom: Top-bottom brick pattern
 * - brick_bottomtop: Bottom-top brick pattern
 * - H_vertical: Vertical H pattern
 * - H_horizontal: Horizontal H pattern
 */
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

/**
 * @brief Represents a rectangulation of n rectangles.
 *
 * Implements the memoryless Gray code algorithm for generating all
 * rectangulations of n rectangles, as described in the paper by
 * Merino and Muetze. The class maintains the data structures for
 * vertices, edges, walls, and rectangles that form the rectangulation.
 */
class Rectangulation {
public:
private:
  /** @brief Number of rectangles in the rectangulation */
  int n_;
  /** @brief Type of rectangulation (generic, diagonal, or baligned) */
  RectangulationType type_;
  /** @brief Forbidden patterns to exclude from generation */
  std::vector<RectangulationPattern> patterns_;
  /** @brief Orientation vector for each rectangle (left or right) */
  std::vector<RectangulationDirection> o_;
  /** @brief State vector for the Gray code algorithm */
  std::vector<int> s_;
  /** @brief Initialize the rectangulation with all vertical edges */
  void set_all_vertical();

 public:
  /** @brief Vector of all vertices in the rectangulation */
  std::vector<Vertex> vertices_;
  /** @brief Vector of all walls in the rectangulation */
  std::vector<Wall> walls_;
  /** @brief Vector of all edges in the rectangulation */
  std::vector<Edge> edges_;
  /** @brief Vector of all rectangles in the rectangulation */
  std::vector<Rectangle> rectangles_;

  /**
   * @brief Construct a rectangulation with n rectangles.
   * @param n Number of rectangles
   * @param type Type of rectangulation (generic, diagonal, or baligned)
   * @param patterns List of forbidden patterns
   */
  Rectangulation(int n, RectangulationType type, std::vector<RectangulationPattern> &patterns);

  /**
   * @brief Initialize rectangulation with pre-built data structures.
   * @param vertices Vector of vertices
   * @param walls Vector of walls
   * @param edges Vector of edges
   * @param rectangles Vector of rectangles
   */
  void init(std::vector<Vertex> &vertices, std::vector<Wall> &walls, std::vector<Edge> &edges,
            std::vector<Rectangle> &rectangles);

  /**
   * @brief Print internal data structures (edges, vertices, walls, rectangles).
   */
  void print_data();

  /**
   * @brief Print rectangle coordinates to stdout.
   *
   * Delegates to print_coordinates_generic or print_coordinates_diagonal
   * based on the rectangulation type.
   */
  void print_coordinates();

  /**
   * @brief Print coordinates using a greedy grid placement algorithm.
   *
   * For generic rectangulations, computes equispaced grid coordinates
   * by sweeping from west to east and south to north.
   */
  void print_coordinates_generic();

  /**
   * @brief Print coordinates using diagonal-based placement.
   *
   * For diagonal and baligned rectangulations, computes coordinates
   * by growing binary trees from diagonal corners.
   */
  void print_coordinates_diagonal();

  /**
   * @brief DFS from bottom-left corner for coordinate assignment.
   * @param vertex_id Starting vertex ID
   * @param val Current value counter
   * @param vertex_x_coord Vector to store x-coordinates
   * @param vertex_y_coord Vector to store y-coordinates
   */
  void DFS_BL(int vertex_id, int &val, std::vector<int> &vertex_x_coord,
              std::vector<int> &vertex_y_coord);

  /**
   * @brief DFS from top-right corner for coordinate assignment.
   * @param vertex_id Starting vertex ID
   * @param val Current value counter
   * @param vertex_x_coord Vector to store x-coordinates
   * @param vertex_y_coord Vector to store y-coordinates
   */
  void DFS_TR(int vertex_id, int &val, std::vector<int> &vertex_x_coord,
              std::vector<int> &vertex_y_coord);

  /**
   * @brief Generate the next rectangulation in Gray code order.
   * @return true if next rectangulation exists, false if at end
   */
  bool next();

  /**
   * @brief Check if rectangle j is bottom-based.
   * @param j Rectangle index
   * @return true if rectangle j is bottom-based
   */
  bool is_bottom_based(int j);

  /**
   * @brief Check if rectangle j is right-based.
   * @param j Rectangle index
   * @return true if rectangle j is right-based
   */
  bool is_right_based(int j);

  // ========== Edge manipulation functions ==========

  /**
   * @brief Remove edge from the head of its current position.
   * @param beta Edge ID to remove
   */
  void remHead(int beta);

  /**
   * @brief Remove edge from the tail of its current position.
   * @param beta Edge ID to remove
   */
  void remTail(int beta);

  /**
   * @brief Insert edge before a given edge in the list.
   * @param beta Edge ID to insert
   * @param a Vertex ID for the new edge's head
   * @param gamma Edge ID to insert before
   */
  void insBefore(int beta, int a, int gamma);

  /**
   * @brief Insert edge after a given edge in the list.
   * @param alpha Edge ID to insert after
   * @param a Vertex ID for the new edge's head
   * @param beta Edge ID to insert
   */
  void insAfter(int alpha, int a, int beta);

  // ========== Jump operations ==========

  /**
   * @brief Perform horizontal W-jump on rectangle j.
   * @param j Rectangle index
   * @param dir Direction (left or right)
   * @param alpha Reference edge ID
   */
  void Wjump_hor(int j, RectangulationDirection dir, int alpha);

  /**
   * @brief Perform vertical W-jump on rectangle j.
   * @param j Rectangle index
   * @param dir Direction (left or right)
   * @param alpha Reference edge ID
   */
  void Wjump_ver(int j, RectangulationDirection dir, int alpha);

  /**
   * @brief Perform S-jump on rectangle j.
   * @param j Rectangle index
   * @param dir Direction (left or right)
   * @param alpha Reference edge ID
   */
  void Sjump(int j, RectangulationDirection dir, int alpha);

  /**
   * @brief Perform horizontal T-jump on rectangle j.
   * @param j Rectangle index
   * @param dir Direction (left or right)
   * @param alpha Reference edge ID
   */
  void Tjump_hor(int j, RectangulationDirection dir, int alpha);

  /**
   * @brief Perform vertical T-jump on rectangle j.
   * @param j Rectangle index
   * @param dir Direction (left or right)
   * @param alpha Reference edge ID
   */
  void Tjump_ver(int j, RectangulationDirection dir, int alpha);

  // ========== Next rectangulation for different types ==========

  /**
   * @brief Generate next generic rectangulation.
   * @param j Rectangle index
   * @param dir Direction
   */
  void next_generic(int j, RectangulationDirection dir);

  /**
   * @brief Generate next diagonal rectangulation.
   * @param j Rectangle index
   * @param dir Direction
   */
  void next_diagonal(int j, RectangulationDirection dir);

  /**
   * @brief Generate next block-aligned rectangulation.
   * @param j Rectangle index
   * @param dir Direction
   */
  void next_baligned(int j, RectangulationDirection dir);

  // ========== Lock/unlock operations ==========

  /**
   * @brief Lock rectangle j in the given direction.
   * @param j Rectangle index
   * @param dir Edge direction to lock
   */
  void lock(int j, EdgeDir dir);

  /**
   * @brief Unlock rectangle j in the given direction.
   * @param j Rectangle index
   * @param dir Direction to unlock
   */
  void unlock(int j, RectangulationDirection dir);

  // ========== Pattern containment checks ==========

  /**
   * @brief Check if rectangle j contains any forbidden pattern.
   * @param j Rectangle index
   * @return true if rectangle j contains a forbidden pattern
   */
  bool contains_pattern(int j);

  /**
   * @brief Check if rectangle j contains a left-right brick pattern.
   * @param j Rectangle index
   * @return true if pattern found
   */
  bool contains_brick_leftright(int j);

  /**
   * @brief Check if rectangle j contains a right-left brick pattern.
   * @param j Rectangle index
   * @return true if pattern found
   */
  bool contains_brick_rightleft(int j);

  /**
   * @brief Check if rectangle j contains a bottom-top brick pattern.
   * @param j Rectangle index
   * @return true if pattern found
   */
  bool contains_brick_bottomtop(int j);

  /**
   * @brief Check if rectangle j contains a top-bottom brick pattern.
   * @param j Rectangle index
   * @return true if pattern found
   */
  bool contains_brick_topbottom(int j);

  /**
   * @brief Check if rectangle j contains a clockwise windmill pattern.
   * @param j Rectangle index
   * @return true if pattern found
   */
  bool contains_wmill_clockwise(int j);

  /**
   * @brief Check if rectangle j contains a counter-clockwise windmill pattern.
   * @param j Rectangle index
   * @return true if pattern found
   */
  bool contains_wmill_counterclockwise(int j);

  /**
   * @brief Check if rectangle j contains a vertical H pattern.
   * @param j Rectangle index
   * @return true if pattern found
   */
  bool contains_H_vertical(int j);

  /**
   * @brief Check if rectangle j contains a horizontal H pattern.
   * @param j Rectangle index
   * @return true if pattern found
   */
  bool contains_H_horizontal(int j);
};
