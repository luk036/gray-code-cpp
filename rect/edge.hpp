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
 * @brief Direction of an edge in the rectangulation.
 *
 * Used to distinguish between:
 * - Hor: Horizontal edge
 * - Ver: Vertical edge
 * - None: Undefined/invalid edge
 */
enum class EdgeDir { Hor, Ver, None };

/**
 * @brief Represents an edge in a rectangulation.
 *
 * Each edge connects two vertices and maintains connectivity information
 * for traversing the rectangulation structure. Edges are stored in a
 * doubly-linked list structure for efficient traversal.
 */
class Edge {
  public:
    /** @brief Direction of the edge (horizontal, vertical, or none) */
    EdgeDir dir_;
    /** @brief Tail vertex ID of this edge */
    int tail_;
    /** @brief Head vertex ID of this edge */
    int head_;
    /** @brief Previous edge ID in the doubly-linked list */
    int prev_;
    /** @brief Next edge ID in the doubly-linked list */
    int next_;
    /** @brief Left rectangle/wall ID */
    int left_;
    /** @brief Right rectangle ID */
    int right_;
    /** @brief Wall ID this edge belongs to */
    int wall_;

    /** @brief Default constructor */
    Edge();

    /**
     * @brief Initialize edge with all properties.
     * @param dir Edge direction (Hor, Ver, or None)
     * @param tail Tail vertex ID
     * @param head Head vertex ID
     * @param prev Previous edge ID in the list
     * @param next Next edge ID in the list
     * @param left Left rectangle/wall ID
     * @param right Right rectangle ID
     * @param wall Wall ID
     */
    void init(EdgeDir dir, int tail, int head, int prev, int next, int left, int right, int wall);
};