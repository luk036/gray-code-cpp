# AGENTS.md - Gray Code C++ Project

Coding guidelines for AI agents in this repository.

## Project Overview

Two independent modules implementing Gray code algorithms:
- **middle**: Middle-level Gray code (Muetze, Nummenpalo research)
- **rect**: Rectangulation generation (Merino, Muetze research)

## Build Commands

### Middle Module
```bash
cd middle
make              # Compile
./middle -h      # Show help
./middle -n2     # Calculate n=2
make clean
```

### Rect Module
```bash
cd rect
make
./rect -h
./rect -n5 -c    # Count solutions
make clean
```

### Run Single Test
```bash
cd middle && make && ./middle -n2 -l10
cd rect && make && ./rect -n3 -c
```

### Build Options
- `-O3`: Optimization (both)
- `-static`: Static linking
- middle: `-std=c++0x`, rect: `-std=c++11`

---

## Code Style

### Indentation
- 4 spaces (NO tabs)
- Brace on same line as declaration

```cpp
class HamCycle {
  public:
    explicit HamCycle(const Vertex &x);
  private:
    Vertex x_;
};

while (true) {
    if (cond) { /* code */ }
}
```

### Headers
- middle: `#ifndef NAME_HPP / #define ... / #endif`
- rect: `#pragma once`
- Include order: local headers first, then stdlib

### Documentation
Doxygen-style for public APIs:
```cpp
/**
 * @brief Construct a Hamilton cycle.
 * @param x Starting vertex
 */
HamCycle::HamCycle(const Vertex &x) : x_(x) { }
```

---

## Naming Conventions

| Type | Style | Example |
|------|-------|---------|
| Classes | PascalCase | `HamCycle`, `Vertex` |
| Members | trailing_underscore | `x_`, `bits_` |
| Functions | camelCase | `is_first_vertex()` |
| Enums | lowercase_underscores | `RectangulationType::generic` |
| Locals | short | `i`, `j`, `k`, `a`, `b` |
| Constants | UPPER_CASE | `NDEBUG` |

---

## Error Handling

- Use `assert()` for preconditions/invariants
- No try-catch (use assertions instead)
- Performance code uses `#ifndef NVISIT`:
```cpp
#ifndef NVISIT
    visit_f_(this->y_.get_bits(), i);
#endif
```

---

## Performance

- Member initializer lists
- Mark functions/params `const` when appropriate
- Default: `-O3` optimization, `-static` linking

---

## Files

### Required Copyright Header (every source file)
```cpp
/*
 * Copyright (c) YYYY Author Name
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 */
```

### License
- middle: GPL v2
- rect: GPL v3

### Adding New Files
1. Create `.cpp` in `middle/` or `rect/`
2. Add copyright header
3. Add to Makefile SOURCES variable

---

## Patterns

### Reference Parameters
```cpp
void init(std::vector<Vertex> &vertices, std::vector<Wall> &walls,
          std::vector<Edge> &edges);
```

### Switch Statement
```cpp
switch (this->type_) {
case (RectangulationType::generic): print_generic(); break;
case (RectangulationType::diagonal): print_diagonal(); break;
}
```

---

## Not Configured

- No clang-format
- No clang-tidy
- No test framework (each module's main.cpp is its own test driver)
- No CI/CD except GitHub Pages