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
 * @brief Represents a wall in a rectangulation.
 *
 * A wall is a connected sequence of edges forming a boundary
 * between rectangles or the outer boundary of the rectangulation.
 * Walls maintain the first and last vertices of the wall edge sequence.
 */
class Wall {
  public:
    /** @brief First vertex ID in the wall */
    int first_;
    /** @brief Last vertex ID in the wall */
    int last_;

    /** @brief Default constructor */
    Wall();

    /**
     * @brief Initialize wall with boundary vertices.
     * @param first First vertex ID in the wall
     * @param last Last vertex ID in the wall
     */
    void init(int first, int last);
};