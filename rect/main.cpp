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

#include "rectangulation.hpp"
#include <cstring>
#include <getopt.h>
#include <iostream>

// display help
void help() {
    std::cout << "./rect [options]   generate various classes of "
                 "rectangulations as described in [Merino,Muetze]"
              << std::endl;
    std::cout << "-h                 display this help" << std::endl;
    std::cout << "-n{1,2,...}        number of rectangles" << std::endl;
    std::cout << "-t{1,2,3}          base type of rectangulations: *1=generic, "
                 "2=diagonal, 3=block-aligned"
              << std::endl;
    std::cout << "-p{1,2,..,8}       forbidden patterns: 1=cw windmill, 2=ccw "
                 "windmill, 3=left/right brick,"
              << std::endl;
    std::cout << "                     4=bottom/top brick, 5=right/left brick, "
                 "6=top/bottom brick, 7=vertical H,"
              << std::endl;
    std::cout << "                     8=horizontal H (see the paper for "
                 "definitions; -p3,...,-p8 unavailable for -t3)"
              << std::endl;
    std::cout
        << "-l{-1,0,1,2,...}   number of rectangulations to list; *-1 for all"
        << std::endl;
    std::cout << "-q                 quiet output" << std::endl;
    std::cout << "-c                 output number of rectangles" << std::endl;
    std::cout << "examples:  ./rect -n5 -c" << std::endl;
    std::cout << "           ./rect -n5 -t2 -c" << std::endl;
    std::cout << "           ./rect -n5 -p3456 -c" << std::endl;
    std::cout << "           ./rect -n10 -t3 -l30" << std::endl;
    std::cout << "           ./rect -n10 -q -c" << std::endl;
}

int main(int argc, char *argv[]) {
    int n;
    bool n_set = false; // flag whether option -n is present
    int t;
    long long steps = -1; // compute all rectangulations by default
    bool quiet = false;   // print output by default
    int c;
    bool output_counts = false; // omit counts by default
    bool patterns3to8 = false;
    const int quiet_dot = 10000000; // print one dot every 10^7 rectangulations in
                              // quiet output mode

    RectangulationType type =
        RectangulationType::generic; // compute generic rectangulations by
                                       // default
    std::vector<RectangulationPattern> patterns;
    while ((c = getopt(argc, argv, "hn:t:p:l:qc")) != -1) {
        switch (c) {
        case 'h':
            help();
            return 0;
        case 'n':
            n = atoi(optarg);
            if (n < 1) {
                std::cerr
                    << "option -n must be followed by an integer from {1,2,...}"
                    << std::endl;
                return 1;
            }
            n_set = true;
            break;
        case 't':
            t = atoi(optarg);
            if ((t < 1) || (t > 3)) {
                std::cerr
                    << "option -t must be followed by an integer from {1,2,3}"
                    << std::endl;
                return 1;
            }
            switch (t) {
            case 1:
                type = RectangulationType::generic;
                break;
            case 2:
                type = RectangulationType::diagonal;
                break;
            case 3:
                type = RectangulationType::baligned;
                break;
            }
            break;
        case 'p':
            for (int i = 0; i < strlen(optarg); i++) {
                const int p = optarg[i] - 48;
                if ((p < 1) || (p > 8)) {
                    std::cerr << "option -p must be followed by an integer "
                                 "from {1,2,...,8}"
                              << std::endl;
                    return 1;
                }
                switch (p) {
                case 1:
                    patterns.push_back(
                        RectangulationPattern::wmill_clockwise);
                    break;
                case 2:
                    patterns.push_back(
                        RectangulationPattern::wmill_counterclockwise);
                    break;
                case 3:
                    patterns.push_back(
                        RectangulationPattern::brick_leftright);
                    patterns3to8 = true;
                    break;
                case 4:
                    patterns.push_back(
                        RectangulationPattern::brick_bottomtop);
                    patterns3to8 = true;
                    break;
                case 5:
                    patterns.push_back(
                        RectangulationPattern::brick_rightleft);
                    patterns3to8 = true;
                    break;
                case 6:
                    patterns.push_back(
                        RectangulationPattern::brick_topbottom);
                    patterns3to8 = true;
                    break;
                case 7:
                    patterns.push_back(RectangulationPattern::H_vertical);
                    patterns3to8 = true;
                    break;
                case 8:
                    patterns.push_back(RectangulationPattern::H_horizontal);
                    patterns3to8 = true;
                    break;
                }
            }
            break;
        case 'l':
            steps = atoi(optarg);
            if (steps < -1) {
                std::cerr << "option -l must be followed by an integer from "
                             "{-1,0,1,2,...}"
                          << std::endl;
                return 1;
            }
            break;
        case 'q':
            quiet = true;
            break;
        case 'c':
            output_counts = true;
            break;
        }
    }
    if (!n_set) {
        std::cerr << "option -n is mandatory" << std::endl;
        help();
        return 1;
    }
    if ((type == RectangulationType::baligned) && patterns3to8) {
        std::cerr << "patterns -p3 to -p8 unavailable for -t3" << std::endl;
        return 1;
    }

    int num_rectangulations = 0;
    Rectangulation rect = Rectangulation(n, type, patterns);

    if (steps == 0) {
        std::cout << "output limit reached" << std::endl;
        return 0;
    }

    bool next;
    do {
        num_rectangulations++;
        if (!quiet) {
            rect.print_coordinates();
        } else if (num_rectangulations % quiet_dot == 0) {
            std::cout << "." << std::flush;
        }

        next = rect.next();
        if (next && (steps >= 0) && (num_rectangulations >= steps)) {
            std::cout << "output limit reached" << std::endl;
            break;
        }
    } while (next);
    if (output_counts) {
        if (quiet && num_rectangulations >= quiet_dot)
            std::cout << std::endl;
        std::cout << "number of rectangulations: " << num_rectangulations
                  << std::endl;
    }

    return 0;
}