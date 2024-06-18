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

#include "vertex.hpp"

Vertex::Vertex() { this->type_ = VertexType::None; }

void Vertex::init(int north, int east, int south, int west) {
    this->north_ = north;
    this->east_ = east;
    this->south_ = south;
    this->west_ = west;
    const int zeros = int(this->north_ == 0) + int(this->south_ == 0) +
                int(this->west_ == 0) + int(this->east_ == 0);
    if (zeros >= 3 or zeros <= 0)
        this->type_ = VertexType::None;
    else if (zeros == 2)
        this->type_ = VertexType::corner;
    else if (this->south_ == 0)
        this->type_ = VertexType::top;
    else if (this->north_ == 0)
        this->type_ = VertexType::bottom;
    else if (this->east_ == 0)
        this->type_ = VertexType::left;
    else if (this->west_ == 0)
        this->type_ = VertexType::right;
}