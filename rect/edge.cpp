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

#include "edge.hpp"

Edge::Edge() {}

void Edge::init(EdgeDir dir, int tail, int head, int prev, int next, int left,
                int right, int wall) {
    this->dir_ = dir;
    this->tail_ = tail;
    this->head_ = head;
    this->left_ = left;
    this->right_ = right;
    this->wall_ = wall;
    this->prev_ = prev;
    this->next_ = next;
}