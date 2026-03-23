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

/**
 * @brief Type of a vertex in the rectangulation.
 *
 * The type indicates where the vertex is located:
 * - corner: Corner vertices of the rectangulation
 * - top: Vertex on the top boundary
 * - bottom: Vertex on the bottom boundary
 * - left: Vertex on the left boundary
 * - right: Vertex on the right boundary
 * - None: Invalid/uninitialized vertex
 */
enum class VertexType { corner, top, bottom, left, right, None };

/**
 * @brief Represents a vertex in a rectangulation.
 *
 * Each vertex is defined by its incident edges (north, east, south, west)
 * and its type indicating its position in the rectangulation structure.
 */
class Vertex {
  public:
    /** @brief Edge ID to the north */
    int north_;
    /** @brief Edge ID to the east */
    int east_;
    /** @brief Edge ID to the south */
    int south_;
    /** @brief Edge ID to the west */
    int west_;
    /** @brief Type of this vertex (corner, top, bottom, left, right, or None) */
    VertexType type_;

    /** @brief Default constructor */
    Vertex();

    /**
     * @brief Initialize vertex with edge IDs.
     * @param north Edge ID to the north
     * @param east Edge ID to the east
     * @param south Edge ID to the south
     * @param west Edge ID to the west
     */
    void init(int north, int east, int south, int west);
};