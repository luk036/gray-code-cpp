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
 * @brief Represents a rectangle in a rectangulation.
 *
 * Each rectangle is defined by four vertices representing its corners:
 * - nwest: Northwest corner vertex ID
 * - neast: Northeast corner vertex ID
 * - swest: Southwest corner vertex ID
 * - seast: Southeast corner vertex ID
 */
class Rectangle {
  public:
    /** @brief Northwest corner vertex ID */
    int nwest_;
    /** @brief Northeast corner vertex ID */
    int neast_;
    /** @brief Southwest corner vertex ID */
    int swest_;
    /** @brief Southeast corner vertex ID */
    int seast_;

    /** @brief Default constructor */
    Rectangle();

    /**
     * @brief Initialize rectangle with corner vertex IDs.
     * @param neast Northeast corner vertex ID
     * @param seast Southeast corner vertex ID
     * @param swest Southwest corner vertex ID
     * @param nwest Northwest corner vertex ID
     */
    void init(int neast, int seast, int swest, int nwest);
};
