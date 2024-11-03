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
};

void Rectangulation::set_all_vertical() {
  // initialize the 3n+1 edges, 2n+2 vertices, n rectangles and n+3 walls
  std::vector<Edge> edges(3 * this->n_ + 2, Edge());
  std::vector<Vertex> vertices(2 * this->n_ + 3, Vertex());
  std::vector<Rectangle> rectangles(this->n_ + 1, Rectangle());
  std::vector<Wall> walls(this->n_ + 4, Wall());
  // set edge, vertex, rectangle and wall "0"
  edges[0].init(EdgeDir::None, 0, 0, 0, 0, 0, 0, 0);
  vertices[0].init(0, 0, 0, 0);
  rectangles[0].init(0, 0, 0, 0);
  walls[0].init(0, 0);
  // set edge properties on the top side of the rectangulation
  for (int i = 1; i <= this->n_ - 1; ++i)
    edges[i].init(EdgeDir::Hor, i, i + 1, i - 1, i + 1, 0, i, 1);
  edges[n_].init(EdgeDir::Hor, this->n_, this->n_ + 1, this->n_ - 1, 0, 0,
                 this->n_, 1);
  // set vertex properties on the top side of the rectangulation
  vertices[1].init(0, 1, 2 * this->n_ + 1, 0);
  vertices[this->n_ + 1].init(0, 0, 3 * this->n_ + 1, this->n_);
  for (int i = 2; i <= this->n_; ++i)
    vertices[i].init(0, i, 2 * this->n_ + i, i - 1);
  // set edge properties on the bottom side of the rectangulation
  for (int i = 2; i <= this->n_ - 1; ++i)
    edges[this->n_ + i].init(EdgeDir::Hor, this->n_ + i + 1, this->n_ + i + 2,
                             this->n_ + i - 1, this->n_ + i + 1, i, 0, 2);
  edges[this->n_ + 1].init(EdgeDir::Hor, this->n_ + 2, this->n_ + 3, 0,
                           this->n_ + 2, 1, 0, 2);
  edges[2 * this->n_].init(EdgeDir::Hor, 2 * this->n_ + 1, 2 * this->n_ + 2,
                           2 * this->n_ - 1, 0, this->n_, 0, 2);
  // set vertex properties on the bottom side of the rectangulation
  vertices[this->n_ + 2].init(2 * this->n_ + 1, this->n_ + 1, 0, 0);
  vertices[2 * this->n_ + 2].init(3 * this->n_ + 1, 0, 0, 2 * this->n_);
  for (int i = 2; i <= this->n_; ++i)
    vertices[this->n_ + i + 1].init(2 * this->n_ + i, this->n_ + i, 0,
                                    this->n_ + i - 1);
  // set edge properties of the vertical edges
  edges[2 * this->n_ + 1].init(EdgeDir::Ver, this->n_ + 2, 1, 0, 0, 0, 1, 3);
  edges[3 * this->n_ + 1].init(EdgeDir::Ver, 2 * this->n_ + 2, this->n_ + 1, 0,
                               0, this->n_, 0, this->n_ + 3);
  for (int i = 2; i <= this->n_; ++i)
    edges[2 * this->n_ + i].init(EdgeDir::Ver, this->n_ + i + 1, i, 0, 0, i - 1,
                                 i, i + 2);
  // set rectangle properties
  for (int i = 1; i <= this->n_; ++i)
    rectangles[i].init(i + 1, this->n_ + i + 2, this->n_ + i + 1, i);
  // set wall parameters
  walls[1].init(1, this->n_ + 1);
  walls[2].init(this->n_ + 2, 2 * this->n_ + 2);
  for (int i = 1; i <= this->n_ + 1; ++i)
    walls[i + 2].init(this->n_ + i + 1, i);
  this->init(vertices, walls, edges, rectangles);
}

void Rectangulation::print_data() {
  // output data structures as in the corresponding figure in the paper
  std::cout << "edges:" << std::endl;
  int i = 0;
  for (auto e : this->edges_) {
    std::cout << "\t" << i << ". ";
    switch (e.dir_) {
    case (EdgeDir::Hor):
      std::cout << "Hor ";
      break;
    case (EdgeDir::Ver):
      std::cout << "Ver ";
      break;
    case (EdgeDir::None):
      std::cout << "None ";
      break;
    }
    std::cout << e.tail_ << " " << e.head_ << " " << e.prev_ << " " << e.next_
              << " " << e.left_ << " " << e.right_ << " " << e.wall_
              << std::endl;
    i++;
  }
  i = 0;
  std::cout << "vertices:" << std::endl;
  for (auto v : this->vertices_) {
    std::cout << "\t" << i << ". " << v.north_ << " " << v.east_ << " "
              << v.south_ << " " << v.west_ << " ";
    switch (v.type_) {
    case (VertexType::None):
      std::cout << "None" << std::endl;
      break;
    case (VertexType::corner):
      std::cout << "corner" << std::endl;
      break;
    case (VertexType::bottom):
      std::cout << "bottom" << std::endl;
      break;
    case (VertexType::top):
      std::cout << "top" << std::endl;
      break;
    case (VertexType::left):
      std::cout << "left" << std::endl;
      break;
    case (VertexType::right):
      std::cout << "right" << std::endl;
      break;
    }
    i++;
  }
  i = 0;
  std::cout << "walls:" << std::endl;
  for (auto w : this->walls_) {
    std::cout << "\t" << i << ". " << w.first_ << " " << w.last_ << std::endl;
    i++;
  }
  i = 0;
  std::cout << "rectangles:" << std::endl;
  for (auto r : this->rectangles_) {
    std::cout << "\t" << i << ". " << r.neast_ << " " << r.seast_ << " "
              << r.swest_ << " " << r.nwest_ << std::endl;
    i++;
  }
}

void Rectangulation::print_coordinates() {
  switch (this->type_) {
  case (RectangulationType::generic):
    print_coordinates_generic();
    break;
  case (RectangulationType::diagonal):
    print_coordinates_diagonal();
    break;
  case (RectangulationType::baligned):
    print_coordinates_diagonal();
    break; // block-aligned rectangulations are a subset of diagonal
           // rectangulations
  }
}

void Rectangulation::print_coordinates_generic() {
  // Run the greedy algorithm for finding an equispaced grid to place the
  // vertices of the rectangles. For finding the x-coordinates of the grid, we
  // do a sweep from west to east.
  std::vector<int> active_vertices;
  std::vector<int> vertex_x_coord(2 * this->n_ + 3, -1);
  // start with every vertex which lies on the western side as an active
  // vertex
  int a = 0;
  for (auto v : this->vertices_) {
    int side_edge_id;
    if (v.type_ == VertexType::right)
      side_edge_id = v.north_;
    else if (v.type_ == VertexType::corner)
      side_edge_id = std::max(v.north_, v.south_);
    else {
      a++;
      continue;
    }
    if (this->edges_[side_edge_id].left_ == 0) {
      active_vertices.push_back(a);
    }
    a++;
  }
  // propagate active vertices
  int x_value = 0;
  while (active_vertices.size() != 0) {
    // set x-coordinate for active vertices
    for (auto idx : active_vertices)
      vertex_x_coord[idx] = x_value;
    x_value++;
    std::vector<int> new_active_vertices;
    // propagate east
    for (auto idx : active_vertices) {
      if (this->vertices_[idx].east_ != 0) {
        const int alpha = this->vertices_[idx].east_;
        new_active_vertices.push_back(this->edges_[alpha].head_);
      }
    }
    // propagate north and south
    std::vector<int> new_active_vertices_copy(new_active_vertices);
    for (auto idx : new_active_vertices_copy) {
      Vertex propagate_from = this->vertices_[idx];
      // propagate north
      while (propagate_from.north_ != 0) {
        Edge e = this->edges_[propagate_from.north_];
        propagate_from = this->vertices_[e.head_];
        new_active_vertices.push_back(e.head_);
      }
      propagate_from = this->vertices_[idx];
      // propagate south
      while (propagate_from.south_ != 0) {
        Edge e = this->edges_[propagate_from.south_];
        propagate_from = this->vertices_[e.tail_];
        new_active_vertices.push_back(e.tail_);
      }
    }
    active_vertices = new_active_vertices;
  }
  // For finding the y-coordinates of the grid, we do the same from south to
  // north.
  active_vertices.clear();
  std::vector<int> vertex_y_coord(2 * this->n_ + 3, -1);
  // start with every vertex which lies on the southern side as an active
  // vertex
  a = 0;
  for (auto v : this->vertices_) {
    int side_edge_id;
    if (v.type_ == VertexType::top)
      side_edge_id = v.east_;
    else if (v.type_ == VertexType::corner)
      side_edge_id = std::max(v.east_, v.west_);
    else {
      a++;
      continue;
    }
    if (this->edges_[side_edge_id].right_ == 0) {
      active_vertices.push_back(a);
    }
    a++;
  }
  // propagate active vertices
  int y_value = 0;
  while (active_vertices.size() != 0) {
    // set y-coordinate for active vertices
    for (auto idx : active_vertices)
      vertex_y_coord[idx] = y_value;
    y_value++;
    std::vector<int> new_active_vertices;
    // propagate north
    for (auto idx : active_vertices) {
      if (this->vertices_[idx].north_ != 0) {
        const int alpha = this->vertices_[idx].north_;
        new_active_vertices.push_back(this->edges_[alpha].head_);
      }
    }
    // propagate west and east
    std::vector<int> new_active_vertices_copy(new_active_vertices);
    for (auto idx : new_active_vertices_copy) {
      Vertex propagate_from = this->vertices_[idx];
      // propagate east
      while (propagate_from.east_ != 0) {
        Edge e = this->edges_[propagate_from.east_];
        propagate_from = this->vertices_[e.head_];
        new_active_vertices.push_back(e.head_);
      }
      propagate_from = this->vertices_[idx];
      // propagate west
      while (propagate_from.west_ != 0) {
        Edge e = this->edges_[propagate_from.west_];
        propagate_from = this->vertices_[e.tail_];
        new_active_vertices.push_back(e.tail_);
      }
    }
    active_vertices = new_active_vertices;
  }
  bool is_first = true;
  bool is_second = true;
  for (auto r : this->rectangles_) {
    if (is_first) {
      is_first = false;
      continue;
    }
    if (!is_second)
      std::cout << " | "; // separator between rectangles
    else
      is_second = false;
    std::cout << vertex_x_coord[r.swest_] << " " << vertex_y_coord[r.swest_]
              << " "; // bottom-left corner
    std::cout << vertex_x_coord[r.neast_] << " "
              << vertex_y_coord[r.neast_]; // top-right corner
  }
  std::cout << std::endl;
}

void Rectangulation::DFS_BL(int vertex_id, int &val,
                            std::vector<int> &vertex_x_coord,
                            std::vector<int> &vertex_y_coord) {
  // We grow the binary tree rooted in the bottom-left corner by a DFS-like
  // procedure. We then check the order in which the diagonal is intersected
  // by the leaves of this tree, and then propagate the value of the vertex
  // coordinates upwards.
  Vertex v = this->vertices_[vertex_id];
  const int top_vertex = this->edges_[v.north_].head_;
  const int right_vertex = this->edges_[v.east_].head_;
  // check if the DFS hits the diagonal while moving up, in that case fix the
  // x-coordinate and increment the value by 1
  if (this->vertices_[top_vertex].type_ == VertexType::corner ||
      this->vertices_[top_vertex].type_ == VertexType::bottom ||
      this->vertices_[top_vertex].type_ == VertexType::left) {
    vertex_x_coord[vertex_id] = val;
    val++;
  }
  // otherwise recurse on the next top vertex
  else {
    DFS_BL(top_vertex, val, vertex_x_coord, vertex_y_coord);
    vertex_x_coord[vertex_id] = vertex_x_coord[top_vertex];
  }
  // check if the DFS hits the diagonal while moving right, in that case fix
  // the y-coordinate and increment the value by 1
  if (this->vertices_[right_vertex].type_ == VertexType::corner ||
      this->vertices_[right_vertex].type_ == VertexType::bottom ||
      this->vertices_[right_vertex].type_ == VertexType::left) {
    vertex_y_coord[vertex_id] = this->n_ - val;
    val++;
  }
  // otherwise recurse on the next right vertex
  else {
    DFS_BL(right_vertex, val, vertex_x_coord, vertex_y_coord);
    vertex_y_coord[vertex_id] = vertex_y_coord[right_vertex];
  }
}
void Rectangulation::DFS_TR(int vertex_id, int &val,
                            std::vector<int> &vertex_x_coord,
                            std::vector<int> &vertex_y_coord) {
  // This function is a symmetric version of DFS_BL, but starting in the
  // top-right corner.
  Vertex v = this->vertices_[vertex_id];
  const int left_vertex = this->edges_[v.west_].tail_;
  const int bottom_vertex = this->edges_[v.south_].tail_;
  if (this->vertices_[bottom_vertex].type_ == VertexType::corner ||
      this->vertices_[bottom_vertex].type_ == VertexType::right ||
      this->vertices_[bottom_vertex].type_ == VertexType::top) {
    vertex_x_coord[vertex_id] = this->n_ - val;
    val++;
  } else {
    DFS_TR(bottom_vertex, val, vertex_x_coord, vertex_y_coord);
    vertex_x_coord[vertex_id] = vertex_x_coord[bottom_vertex];
  }
  if (this->vertices_[left_vertex].type_ == VertexType::corner ||
      this->vertices_[left_vertex].type_ == VertexType::right ||
      this->vertices_[left_vertex].type_ == VertexType::top) {
    vertex_y_coord[vertex_id] = val;
    val++;
  } else {
    DFS_TR(left_vertex, val, vertex_x_coord, vertex_y_coord);
    vertex_y_coord[vertex_id] = vertex_y_coord[left_vertex];
  }
}

void Rectangulation::print_coordinates_diagonal() {
  std::vector<int> vertex_x_coord(2 * this->n_ + 3, -1);
  std::vector<int> vertex_y_coord(2 * this->n_ + 3, -1);
  // find bottom-left corner and top-right corner ids
  int BL = -1;
  int TR = -1;
  for (int i = 1; i < 2 * this->n_ + 3; ++i) {
    Vertex v = this->vertices_[i];
    if (v.north_ == 0 && v.east_ == 0 && v.type_ == VertexType::corner)
      TR = i;
    else if (v.south_ == 0 && v.west_ == 0 && v.type_ == VertexType::corner)
      BL = i;
  }
  assert(BL != -1 && TR != -1);
  int val = 0;
  DFS_BL(BL, val, vertex_x_coord, vertex_y_coord);
  val = 0;
  DFS_TR(TR, val, vertex_x_coord, vertex_y_coord);
  bool is_first = true;
  bool is_second = true;
  for (auto r : this->rectangles_) {
    if (is_first) {
      is_first = false;
      continue; // skip non-used 0th rectangle
    }
    if (!is_second)
      std::cout << " | "; // don't print separator before first rectangle
    else
      is_second = false;
    std::cout << vertex_x_coord[r.swest_] << " " << vertex_y_coord[r.swest_]
              << " ";
    std::cout << vertex_x_coord[r.neast_] << " " << vertex_y_coord[r.neast_];
  }
  std::cout << std::endl;
}

bool Rectangulation::next() {
  // Run one iteration of the memoryless algorithm, return true if the next
  // rectangulation is not the identity, false otherwise. M3 - Select
  // rectangle
  const int j = this->s_[n_];
  if ((j == 1) ||
      ((this->n_ == 2) && (this->type_ == RectangulationType::baligned)))
    return false;
  // M4 - Jump rectangle
  switch (this->type_) {
  case (RectangulationType::generic):
    next_generic(j, this->o_[j]);
    while (contains_pattern(j))
      next_generic(j, this->o_[j]);
    break;
  case (RectangulationType::diagonal):
    next_diagonal(j, this->o_[j]);
    while (contains_pattern(j))
      next_diagonal(j, this->o_[j]);
    break;
  case (RectangulationType::baligned):
    next_baligned(j, this->o_[j]);
    while (contains_pattern(j))
      next_baligned(j, this->o_[j]);
    break;
  }
  // M5 - Update o and s
  this->s_[n_] = n_;
  if ((this->type_ == RectangulationType::baligned) &&
      (this->o_[j - 1] == RectangulationDirection::left) &&
      is_bottom_based(j - 1)) {
    this->o_[j - 1] = RectangulationDirection::right;
    this->s_[j - 1] = this->s_[j - 2];
    this->s_[j - 2] = j - 2;
  }
  if (this->o_[j] == RectangulationDirection::left && is_bottom_based(j)) {
    this->o_[j] = RectangulationDirection::right;
    this->s_[j] = this->s_[j - 1];
    this->s_[j - 1] = j - 1;
  }
  if ((this->type_ == RectangulationType::baligned) &&
      (this->o_[j - 1] == RectangulationDirection::right) &&
      is_right_based(j - 1)) {
    this->o_[j - 1] = RectangulationDirection::left;
    this->s_[j - 1] = this->s_[j - 2];
    this->s_[j - 2] = j - 2;
  }
  if ((this->o_[j] == RectangulationDirection::right) && is_right_based(j)) {
    this->o_[j] = RectangulationDirection::left;
    this->s_[j] = this->s_[j - 1];
    this->s_[j - 1] = j - 1;
  }
  return true;
}

bool Rectangulation::is_bottom_based(int j) {
  const int a = this->rectangles_[j].nwest_;
  const int alpha = this->vertices_[a].south_;
  const int b = this->rectangles_[j].swest_;
  if (this->edges_[alpha].left_ == 0)
    return true;
  if (this->type_ == RectangulationType::baligned &&
      a == this->rectangles_[j - 1].neast_ &&
      b == this->rectangles_[j - 1].seast_) {
    const int c = this->rectangles_[j - 1].nwest_;
    const int gamma = this->vertices_[c].south_;
    if (this->edges_[gamma].left_ == 0)
      return true;
  }
  return false;
}

bool Rectangulation::is_right_based(int j) {
  const int a = this->rectangles_[j].nwest_;
  const int alpha = this->vertices_[a].east_;
  const int b = this->rectangles_[j].neast_;
  if (this->edges_[alpha].left_ == 0)
    return true;
  if (this->type_ == RectangulationType::baligned &&
      a == this->rectangles_[j - 1].swest_ &&
      b == this->rectangles_[j - 1].seast_) {
    const int c = this->rectangles_[j].nwest_;
    const int gamma = this->vertices_[c].east_;
    if (this->edges_[gamma].left_ == 0)
      return true;
  }
  return false;
}

void Rectangulation::remHead(int beta) {
  // Step 1 - Prepare
  const int alpha = this->edges_[beta].prev_;
  const int gamma = this->edges_[beta].next_;
  const int a = this->edges_[beta].tail_;
  // Step 2 - Update edges/vertices
  if (alpha != 0)
    this->edges_[alpha].next_ = gamma;
  if (gamma != 0)
    this->edges_[gamma].prev_ = alpha;
  this->edges_[gamma].tail_ = a;
  if (this->edges_[beta].dir_ == EdgeDir::Hor)
    this->vertices_[a].east_ = gamma;
  else
    this->vertices_[a].north_ = gamma;
  // Step 3 - Update wall
  const int x = this->edges_[beta].wall_;
  if (this->edges_[beta].head_ == this->walls_[x].last_)
    this->walls_[x].last_ = a;
}

void Rectangulation::remTail(int beta) {
  // Step 1 - Prepare
  const int alpha = this->edges_[beta].prev_;
  const int gamma = this->edges_[beta].next_;
  const int a = this->edges_[beta].head_;
  // Step 2 - Update edges/vertices
  if (alpha != 0)
    this->edges_[alpha].next_ = gamma;
  if (gamma != 0)
    this->edges_[gamma].prev_ = alpha;
  this->edges_[alpha].head_ = a;
  if (this->edges_[beta].dir_ == EdgeDir::Hor)
    this->vertices_[a].west_ = alpha;
  else
    this->vertices_[a].south_ = alpha;
  // Step 3 - Update wall
  const int x = this->edges_[beta].wall_;
  if (this->edges_[beta].tail_ == this->walls_[x].first_)
    this->walls_[x].first_ = a;
}

void Rectangulation::insBefore(int beta, int a, int gamma) {
  // Step 1 - Prepare
  const int alpha = this->edges_[gamma].prev_;
  const int b = this->edges_[gamma].tail_;
  // Step 2 - Update edges/vertices
  this->edges_[beta].tail_ = b;
  this->edges_[beta].head_ = a;
  this->edges_[beta].prev_ = alpha;
  this->edges_[beta].next_ = gamma;
  this->edges_[gamma].tail_ = a;
  this->edges_[gamma].prev_ = beta;
  if (alpha != 0)
    this->edges_[alpha].next_ = beta;
  if (this->edges_[gamma].dir_ == EdgeDir::Hor) {
    this->edges_[beta].dir_ = EdgeDir::Hor;
    this->vertices_[a].west_ = beta;
    this->vertices_[a].east_ = gamma;
    this->vertices_[b].east_ = beta;
  } else {
    assert(this->edges_[gamma].dir_ == EdgeDir::Ver);
    this->edges_[beta].dir_ = EdgeDir::Ver;
    this->vertices_[a].south_ = beta;
    this->vertices_[a].north_ = gamma;
    this->vertices_[b].north_ = beta;
  }
  // Step 3 - Update wall
  this->edges_[beta].wall_ = this->edges_[gamma].wall_;
}

void Rectangulation::insAfter(int alpha, int a, int beta) {
  // Step 1 - Prepare
  const int gamma = this->edges_[alpha].next_;
  const int b = this->edges_[alpha].head_;
  // Step 2 - Updates edges/vertices
  this->edges_[beta].tail_ = a;
  this->edges_[beta].head_ = b;
  this->edges_[beta].prev_ = alpha;
  this->edges_[beta].next_ = gamma;
  this->edges_[alpha].head_ = a;
  this->edges_[alpha].next_ = beta;
  if (gamma != 0)
    this->edges_[gamma].prev_ = beta;
  if (this->edges_[alpha].dir_ == EdgeDir::Hor) {
    this->edges_[beta].dir_ = EdgeDir::Hor;
    this->vertices_[a].west_ = alpha;
    this->vertices_[a].east_ = beta;
    this->vertices_[b].west_ = beta;
  } else {
    assert(this->edges_[alpha].dir_ == EdgeDir::Ver);
    this->edges_[beta].dir_ = EdgeDir::Ver;
    this->vertices_[a].south_ = alpha;
    this->vertices_[a].north_ = beta;
    this->vertices_[b].south_ = beta;
  }
  // Step 3 - Update wall
  this->edges_[beta].wall_ = this->edges_[alpha].wall_;
}

void Rectangulation::Wjump_hor(int j, RectangulationDirection dir, int alpha) {
  if (dir == RectangulationDirection::left) {
    // Step 1 - Prepare
    const int a = this->rectangles_[j].nwest_;
    const int beta = this->vertices_[a].west_;
    assert(alpha == this->edges_[beta].prev_);
    const int k = this->edges_[alpha].left_;
    // Step 2 - Flip and update rectangles
    remHead(beta);
    insAfter(alpha, a, beta);
    this->edges_[beta].left_ = k;
    this->edges_[beta].right_ = j;
  } else {
    assert(dir == RectangulationDirection::right);
    // Step 1 - Prepare
    const int a = this->rectangles_[j].nwest_;
    const int beta = this->vertices_[a].east_;
    assert(alpha == this->edges_[beta].next_);
    const int k = this->edges_[alpha].left_;
    const int alpha_prime = this->edges_[beta].prev_;
    const int l = this->edges_[alpha_prime].right_;
    // Step 2 - Flip and update rectangles
    remTail(beta);
    insBefore(beta, a, alpha);
    this->edges_[beta].left_ = k;
    this->edges_[beta].right_ = l;
  }
}

void Rectangulation::Wjump_ver(int j, RectangulationDirection dir, int alpha) {
  if (dir == RectangulationDirection::right) {
    // Step 1 - Prepare
    const int a = this->rectangles_[j].nwest_;
    const int beta = this->vertices_[a].north_;
    assert(alpha == this->edges_[beta].next_);
    const int k = this->edges_[alpha].left_;
    // Step 2 - Flip and update rectangles
    remTail(beta);
    insBefore(beta, a, alpha);
    this->edges_[beta].left_ = k;
    this->edges_[beta].right_ = j;
  } else {
    assert(dir == RectangulationDirection::left);
    // Step 1 - Prepare
    const int a = this->rectangles_[j].nwest_;
    const int beta = this->vertices_[a].south_;
    assert(alpha == this->edges_[beta].prev_);
    const int k = this->edges_[alpha].left_;
    const int alpha_prime = this->edges_[beta].next_;
    const int l = this->edges_[alpha_prime].right_;
    // Step 2 - Flip and update rectangles
    remHead(beta);
    insAfter(alpha, a, beta);
    this->edges_[beta].left_ = k;
    this->edges_[beta].right_ = l;
  }
}

void Rectangulation::Sjump(int j, RectangulationDirection d, int alpha) {
  if (d == RectangulationDirection::left) {
    // Step 1 - Prepare
    const int a = this->rectangles_[j].nwest_;
    const int b = this->rectangles_[j].swest_;
    const int c = this->rectangles_[j].neast_;
    const int alpha_prime = this->vertices_[a].west_;
    const int beta = this->vertices_[a].east_;
    const int beta_prime = this->vertices_[b].west_;
    const int gamma = this->vertices_[c].south_;
    const int delta = this->vertices_[a].south_;
    const int c_prime = this->edges_[beta_prime].tail_;
    const int k = this->edges_[alpha].left_;
    const int l = this->edges_[gamma].right_;
    const int x = this->edges_[delta].wall_;
    // Step 2 - Flip
    remTail(beta);
    remHead(beta_prime);
    insBefore(beta, a, alpha);
    insAfter(gamma, b, beta_prime);
    this->edges_[delta].dir_ = EdgeDir::Hor;
    this->edges_[delta].tail_ = a;
    this->edges_[delta].head_ = b;
    this->vertices_[a].east_ = delta;
    this->vertices_[a].west_ = 0;
    this->vertices_[a].type_ = VertexType::right;
    this->vertices_[b].east_ = 0;
    this->vertices_[b].west_ = delta;
    this->vertices_[b].type_ = VertexType::left;
    this->walls_[x].first_ = a;
    this->walls_[x].last_ = b;
    // Step 3 - Update rectangles
    this->rectangles_[j].neast_ = b;
    this->rectangles_[j].swest_ = c_prime;
    this->rectangles_[j - 1].neast_ = c;
    this->rectangles_[j - 1].swest_ = a;
    int nu = this->vertices_[c].west_;
    while (nu != alpha_prime) {
      this->edges_[nu].right_ = j - 1;
      nu = this->edges_[nu].prev_;
    }
    nu = this->vertices_[c_prime].north_;
    while (nu != alpha) {
      this->edges_[nu].right_ = j;
      nu = this->edges_[nu].next_;
    }
    this->edges_[beta].left_ = k;
    this->edges_[beta].right_ = j;
    this->edges_[beta_prime].left_ = j - 1;
    this->edges_[beta_prime].right_ = l;
  } else {
    assert(d == RectangulationDirection::right);
    // this code is obtained by mirroring the case RectangulationDirection::left
    // along the main diagonal Step 1 - Prepare
    const int a = this->rectangles_[j].nwest_;
    const int b = this->rectangles_[j].neast_;
    const int c = this->rectangles_[j].swest_;
    const int alpha_prime = this->vertices_[a].north_;
    const int beta = this->vertices_[a].south_;
    const int beta_prime = this->vertices_[b].north_;
    const int gamma = this->vertices_[c].east_;
    const int delta = this->vertices_[a].east_;
    const int c_prime = this->edges_[beta_prime].head_;
    const int k = this->edges_[alpha].left_;
    const int l = this->edges_[gamma].right_;
    const int x = this->edges_[delta].wall_;
    // Step 2 - Flip
    remHead(beta);
    remTail(beta_prime);
    insAfter(alpha, a, beta);
    insBefore(beta_prime, b, gamma);
    this->edges_[delta].dir_ = EdgeDir::Ver;
    this->edges_[delta].head_ = a;
    this->edges_[delta].tail_ = b;
    this->vertices_[a].south_ = delta;
    this->vertices_[a].north_ = 0;
    this->vertices_[a].type_ = VertexType::bottom;
    this->vertices_[b].south_ = 0;
    this->vertices_[b].north_ = delta;
    this->vertices_[b].type_ = VertexType::top;
    this->walls_[x].last_ = a;
    this->walls_[x].first_ = b;
    // Step 3 - Update rectangles
    this->rectangles_[j].swest_ = b;
    this->rectangles_[j].neast_ = c_prime;
    this->rectangles_[j - 1].swest_ = c;
    this->rectangles_[j - 1].neast_ = a;
    int nu = this->vertices_[c].north_;
    while (nu != alpha_prime) {
      this->edges_[nu].right_ = j - 1;
      nu = this->edges_[nu].next_;
    }
    nu = this->vertices_[c_prime].west_;
    while (nu != alpha) {
      this->edges_[nu].right_ = j;
      nu = this->edges_[nu].prev_;
    }
    this->edges_[beta].left_ = k;
    this->edges_[beta].right_ = j;
    this->edges_[beta_prime].left_ = j - 1;
    this->edges_[beta_prime].right_ = l;
  }
}

void Rectangulation::Tjump_hor(int j, RectangulationDirection dir, int alpha) {
  if (dir == RectangulationDirection::left) {
    // Step 1 - Prepare
    const int a = this->rectangles_[j].nwest_;
    const int b = this->edges_[alpha].head_;
    const int c = this->rectangles_[j].neast_;
    const int alpha_prime = this->vertices_[a].west_;
    const int beta = this->vertices_[a].east_;
    const int beta_prime = this->vertices_[a].south_;
    const int gamma = this->vertices_[c].south_;
    const int gamma_prime = this->vertices_[b].south_;
    const int k = this->edges_[beta_prime].left_;
    const int l = this->edges_[gamma].right_;
    const int m = this->edges_[alpha].right_;
    const int x = this->edges_[alpha].wall_;
    const int y = this->edges_[gamma_prime].wall_;
    // Step 2 - Flip
    remTail(beta);
    remTail(beta_prime);
    insAfter(alpha, a, beta);
    insAfter(gamma, b, beta_prime);
    this->edges_[beta].head_ = b;
    this->edges_[gamma_prime].head_ = a;
    this->vertices_[a].south_ = gamma_prime;
    this->vertices_[b].west_ = beta;
    this->walls_[x].last_ = b;
    this->walls_[y].last_ = a;
    // Step 3 - Update rectangles
    this->rectangles_[j].neast_ = b;
    this->rectangles_[k].neast_ = c;
    this->rectangles_[m].neast_ = a;
    int nu = this->vertices_[c].west_;
    while (nu != alpha_prime) {
      this->edges_[nu].right_ = k;
      nu = this->edges_[nu].prev_;
    }
    this->edges_[beta].left_ = k;
    this->edges_[beta_prime].right_ = l;
  } else {
    assert(dir == RectangulationDirection::right);
    // Step 1 - Prepare
    const int a = this->rectangles_[j].nwest_;
    const int b = this->rectangles_[j].neast_;
    const int alpha_prime = this->vertices_[a].west_;
    const int gamma_prime = this->vertices_[a].south_;
    const int beta = this->vertices_[a].east_;
    const int beta_prime = this->vertices_[b].north_;
    const int c = this->edges_[beta_prime].head_;
    const int k = this->edges_[beta].left_;
    const int l = this->edges_[alpha].left_;
    const int m = this->edges_[alpha_prime].right_;
    const int x = this->edges_[alpha_prime].wall_;
    const int y = this->edges_[gamma_prime].wall_;
    // Step 2 - Flip
    remTail(beta);
    remTail(beta_prime);
    insAfter(alpha, a, beta);
    insAfter(gamma_prime, b, beta_prime);
    this->edges_[alpha_prime].head_ = b;
    this->edges_[beta_prime].head_ = a;
    this->vertices_[a].south_ = beta_prime;
    this->vertices_[b].west_ = alpha_prime;
    this->walls_[x].last_ = b;
    this->walls_[y].last_ = a;
    // Step 3 - Update rectangles
    this->rectangles_[j].neast_ = c;
    this->rectangles_[k].neast_ = a;
    this->rectangles_[m].neast_ = b;
    int nu = this->vertices_[c].west_;
    while (nu != alpha) {
      this->edges_[nu].right_ = j;
      nu = this->edges_[nu].prev_;
    }
    this->edges_[beta].left_ = l;
    this->edges_[beta_prime].right_ = j;
  }
}

void Rectangulation::Tjump_ver(int j, RectangulationDirection dir, int alpha) {
  // this code is obtained by mirroring Tjump_hor along the main diagonal
  if (dir == RectangulationDirection::right) {
    // Step 1 - Prepare
    const int a = this->rectangles_[j].nwest_;
    const int b = this->edges_[alpha].tail_;
    const int c = this->rectangles_[j].swest_;
    const int alpha_prime = this->vertices_[a].north_;
    const int beta = this->vertices_[a].south_;
    const int beta_prime = this->vertices_[a].east_;
    const int gamma = this->vertices_[c].east_;
    const int gamma_prime = this->vertices_[b].east_;
    const int k = this->edges_[beta_prime].left_;
    const int l = this->edges_[gamma].right_;
    const int m = this->edges_[alpha].right_;
    const int x = this->edges_[alpha].wall_;
    const int y = this->edges_[gamma_prime].wall_;
    // Step 2 - Flip
    remHead(beta);
    remHead(beta_prime);
    insBefore(beta, a, alpha);
    insBefore(beta_prime, b, gamma);
    this->edges_[beta].tail_ = b;
    this->edges_[gamma_prime].tail_ = a;
    this->vertices_[a].east_ = gamma_prime;
    this->vertices_[a].type_ = VertexType::right;
    this->vertices_[b].north_ = beta;
    this->vertices_[b].type_ = VertexType::top;
    this->walls_[x].first_ = b;
    this->walls_[y].first_ = a;
    // Step 3 - Update rectangles
    this->rectangles_[j].swest_ = b;
    this->rectangles_[k].swest_ = c;
    this->rectangles_[m].swest_ = a;
    int nu = this->vertices_[c].north_;
    while (nu != alpha_prime) {
      this->edges_[nu].right_ = k;
      nu = this->edges_[nu].next_;
    }
    this->edges_[beta].left_ = k;
    this->edges_[beta_prime].right_ = l;
  } else {
    assert(dir == RectangulationDirection::left);
    // Step 1 - Prepare
    const int a = this->rectangles_[j].nwest_;
    const int b = this->rectangles_[j].swest_;
    const int alpha_prime = this->vertices_[a].north_;
    const int gamma_prime = this->vertices_[a].east_;
    const int beta = this->vertices_[a].south_;
    const int beta_prime = this->vertices_[b].west_;
    const int c = this->edges_[beta_prime].tail_;
    const int k = this->edges_[beta].left_;
    const int l = this->edges_[alpha].left_;
    const int m = this->edges_[alpha_prime].right_;
    const int x = this->edges_[alpha_prime].wall_;
    const int y = this->edges_[gamma_prime].wall_;
    // Step 2 - Flip
    remHead(beta);
    remHead(beta_prime);
    insBefore(beta, a, alpha);
    insBefore(beta_prime, b, gamma_prime);
    this->edges_[alpha_prime].tail_ = b;
    this->edges_[beta_prime].tail_ = a;
    this->vertices_[a].east_ = beta_prime;
    this->vertices_[a].type_ = VertexType::right;
    this->vertices_[b].north_ = alpha_prime;
    this->vertices_[b].type_ = VertexType::top;
    this->walls_[x].first_ = b;
    this->walls_[y].first_ = a;
    // Step 3 - Update rectangles
    this->rectangles_[j].swest_ = c;
    this->rectangles_[k].swest_ = a;
    this->rectangles_[m].swest_ = b;
    int nu = this->vertices_[c].north_;
    while (nu != alpha) {
      this->edges_[nu].right_ = j;
      nu = this->edges_[nu].next_;
    }
    this->edges_[beta].left_ = l;
    this->edges_[beta_prime].right_ = j;
  }
}

void Rectangulation::next_generic(int j, RectangulationDirection dir) {
  // N1 - Prepare
  const int a = this->rectangles_[j].nwest_;
  if (dir == RectangulationDirection::left &&
      this->vertices_[a].type_ == VertexType::bottom) {
    // Preparing horizontal left jump
    const int alpha = this->vertices_[a].west_;
    const int beta = this->vertices_[a].south_;
    const int b = this->edges_[beta].tail_;
    const int c = this->edges_[alpha].tail_;
    // N2 - Execute horizontal left jump
    if (this->vertices_[c].type_ == VertexType::top) {
      const int gamma = this->vertices_[c].west_;
      Wjump_hor(j, RectangulationDirection::left, gamma);
    } else if (this->vertices_[b].type_ == VertexType::left) {
      const int gamma = this->vertices_[b].west_;
      Tjump_hor(j, RectangulationDirection::left, gamma);
    } else {
      assert(this->vertices_[b].type_ == VertexType::top);
      const int gamma = this->vertices_[c].south_;
      Sjump(j, RectangulationDirection::left, gamma);
    }
  } else if (dir == RectangulationDirection::right &&
             this->vertices_[a].type_ == VertexType::bottom) {
    // Prepare horizontal right jump
    const int alpha = this->vertices_[a].east_;
    const int b = this->edges_[alpha].head_;
    // N3 - Execute horizontal right jump
    if (this->vertices_[b].type_ == VertexType::top) {
      const int gamma = this->vertices_[b].east_;
      Wjump_hor(j, RectangulationDirection::right, gamma);
    } else {
      assert(this->vertices_[b].type_ == VertexType::left);
      const int k = this->edges_[alpha].left_;
      const int c = this->rectangles_[k].nwest_;
      const int gamma = this->vertices_[c].east_;
      Tjump_hor(j, RectangulationDirection::right, gamma);
    }
  } else if (dir == RectangulationDirection::right &&
             this->vertices_[a].type_ == VertexType::right) {
    // Prepare vertical right jump
    const int alpha = this->vertices_[a].north_;
    const int beta = this->vertices_[a].east_;
    const int b = this->edges_[beta].head_;
    const int c = this->edges_[alpha].head_;
    // N4 - Execute vertical right jump
    if (this->vertices_[c].type_ == VertexType::left) {
      const int gamma = this->vertices_[c].north_;
      Wjump_ver(j, RectangulationDirection::right, gamma);
    } else if (this->vertices_[b].type_ == VertexType::top) {
      const int gamma = this->vertices_[b].north_;
      Tjump_ver(j, RectangulationDirection::right, gamma);
    } else {
      assert(this->vertices_[b].type_ == VertexType::left);
      const int gamma = this->vertices_[c].east_;
      Sjump(j, RectangulationDirection::right, gamma);
    }
  } else {
    assert(dir == RectangulationDirection::left &&
           this->vertices_[a].type_ == VertexType::right);
    // Prepare vertical left jump
    const int alpha = this->vertices_[a].south_;
    const int b = this->edges_[alpha].tail_;
    // N5 - Execute vertical left jump
    if (this->vertices_[b].type_ == VertexType::left) {
      const int gamma = this->vertices_[b].south_;
      Wjump_ver(j, RectangulationDirection::left, gamma);
    } else {
      assert(this->vertices_[b].type_ == VertexType::top);
      const int k = this->edges_[alpha].left_;
      const int c = this->rectangles_[k].nwest_;
      const int gamma = this->vertices_[c].south_;
      Tjump_ver(j, RectangulationDirection::left, gamma);
    }
  }
}

void Rectangulation::next_diagonal(int j, RectangulationDirection dir) {
  // N1 - Prepare
  const int a = this->rectangles_[j].nwest_;
  if (dir == RectangulationDirection::left &&
      this->vertices_[a].type_ == VertexType::bottom) {
    // Prepare horizontal left jump
    const int alpha = this->vertices_[a].south_;
    const int b = this->edges_[alpha].tail_;
    // N2 - Horizontal left jump
    if (this->vertices_[b].type_ == VertexType::left) {
      const int gamma = this->vertices_[b].west_;
      Tjump_hor(j, RectangulationDirection::left, gamma);
    } else {
      assert(this->vertices_[b].type_ == VertexType::top);
      const int c = this->rectangles_[j - 1].swest_;
      const int gamma = this->vertices_[c].north_;
      Sjump(j, RectangulationDirection::left, gamma);
    }
  } else if (dir == RectangulationDirection::right &&
             this->vertices_[a].type_ == VertexType::bottom) {
    // Prepare horizontal right jump
    const int alpha = this->vertices_[a].east_;
    // N3 - Horizontal right jump
    const int k = this->edges_[alpha].left_;
    const int b = this->rectangles_[k].neast_;
    const int gamma = this->vertices_[b].west_;
    Tjump_hor(j, RectangulationDirection::right, gamma);
  } else if (dir == RectangulationDirection::right &&
             this->vertices_[a].type_ == VertexType::right) {
    // Prepare vertical right jump
    const int alpha = this->vertices_[a].east_;
    const int b = this->edges_[alpha].head_;
    // N4 - Vertical right jump
    if (this->vertices_[b].type_ == VertexType::top) {
      const int gamma = this->vertices_[b].north_;
      Tjump_ver(j, RectangulationDirection::right, gamma);
    } else {
      assert(this->vertices_[b].type_ == VertexType::left);
      const int c = this->rectangles_[j - 1].neast_;
      const int gamma = this->vertices_[c].west_;
      Sjump(j, RectangulationDirection::right, gamma);
    }
  } else {
    assert(dir == RectangulationDirection::left &&
           this->vertices_[a].type_ == VertexType::right);
    // Prepare vertical left jump
    const int alpha = this->vertices_[a].south_;
    // N5 - Vertical left jump
    const int k = this->edges_[alpha].left_;
    const int b = this->rectangles_[k].swest_;
    const int gamma = this->vertices_[b].north_;
    Tjump_ver(j, RectangulationDirection::left, gamma);
  }
}

void Rectangulation::next_baligned(int j, RectangulationDirection dir) {
  // N1 - Prepare
  int a = this->rectangles_[j].nwest_;
  unlock(j, dir);
  if (dir == RectangulationDirection::left &&
      this->vertices_[a].type_ == VertexType::bottom) {
    int alpha = this->vertices_[a].south_;
    int b = this->edges_[alpha].tail_;
    if (this->vertices_[b].type_ == VertexType::left) {
      // N2 - Horizontal left jump (T/TS)
      int gamma = this->vertices_[b].west_;
      Tjump_hor(j, RectangulationDirection::left, gamma);
      a = this->rectangles_[j].nwest_;
      alpha = this->vertices_[a].south_;
      b = this->edges_[alpha].tail_;
      const int c = this->rectangles_[j - 1].swest_;
      gamma = this->vertices_[c].north_;
      const int c_prime = this->rectangles_[j].seast_;
      if (this->vertices_[b].type_ == VertexType::top &&
          (this->vertices_[c_prime].type_ == VertexType::left ||
           (j == this->n_ && this->edges_[gamma].left_ == 0)))
        Sjump(j, RectangulationDirection::left, gamma);
      lock(j, EdgeDir::Hor);
    } else {
      assert(this->vertices_[b].type_ == VertexType::top);
      // N3 - Horizontal left jump (ST/D)
      int c = this->rectangles_[j - 1].swest_;
      int gamma = this->vertices_[c].north_;
      Sjump(j, RectangulationDirection::left, gamma);
      gamma = this->vertices_[c].north_;
      const int k = this->edges_[gamma].left_;
      int c_prime = this->rectangles_[k].swest_;
      int gamma_prime = this->vertices_[c_prime].north_;
      Tjump_ver(j, RectangulationDirection::left, gamma_prime);
      c = this->rectangles_[j - 1].swest_;
      gamma = this->vertices_[c].north_;
      const int a = this->edges_[gamma].head_;
      if (this->vertices_[a].type_ == VertexType::bottom) {
        assert(k == j - 2);
        c_prime = this->rectangles_[j - 2].swest_;
        gamma_prime = this->vertices_[c_prime].north_;
        Sjump(j - 1, RectangulationDirection::left, gamma_prime);
      }
      lock(j - 1, EdgeDir::Hor);
    }
  } else if (dir == RectangulationDirection::right &&
             this->vertices_[a].type_ == VertexType::bottom) {
    int alpha = this->vertices_[a].east_;
    // N4 - Horizontal right jump (T/TS)
    const int k = this->edges_[alpha].left_;
    int b = this->rectangles_[k].neast_;
    int gamma = this->vertices_[b].west_;
    Tjump_hor(j, RectangulationDirection::right, gamma);
    a = this->rectangles_[j].nwest_;
    alpha = this->vertices_[a].south_;
    b = this->edges_[alpha].tail_;
    const int beta = this->vertices_[b].south_;
    gamma = this->vertices_[b].west_;
    const int c = this->edges_[beta].tail_;
    const int c_prime = this->edges_[gamma].tail_;
    if (this->vertices_[c].type_ == VertexType::top &&
        this->vertices_[c_prime].type_ == VertexType::right) {
      assert(k == j - 2);
      const int gamma_prime = this->vertices_[a].west_;
      Sjump(j - 1, RectangulationDirection::right, gamma_prime);
    }
    lock(j, EdgeDir::Ver);
  } else if (dir == RectangulationDirection::right &&
             this->vertices_[a].type_ == VertexType::right) {
    int alpha = this->vertices_[a].east_;
    int b = this->edges_[alpha].head_;
    if (this->vertices_[b].type_ == VertexType::top) {
      // N5 - Vertical right jump (T/TS)
      int gamma = this->vertices_[b].north_;
      Tjump_ver(j, RectangulationDirection::right, gamma);
      a = this->rectangles_[j].nwest_;
      alpha = this->vertices_[a].east_;
      b = this->edges_[alpha].head_;
      const int c = this->rectangles_[j - 1].neast_;
      gamma = this->vertices_[c].west_;
      const int c_prime = this->rectangles_[j].seast_;
      const int e = this->rectangles_[j - 1].nwest_;
      if (this->vertices_[b].type_ == VertexType::left &&
          (this->vertices_[c_prime].type_ == VertexType::top ||
           (j == this->n_ and
            !(this->vertices_[e].type_ == VertexType::right and
              this->edges_[gamma].tail_ == e))))
        Sjump(j, RectangulationDirection::right, gamma);
      lock(j, EdgeDir::Ver);
    } else {
      assert(this->vertices_[b].type_ == VertexType::left);
      // N6 - Vertical right jump (ST/D)
      int c = this->rectangles_[j - 1].neast_;
      int gamma = this->vertices_[c].west_;
      Sjump(j, RectangulationDirection::right, gamma);
      gamma = this->vertices_[c].west_;
      const int k = this->edges_[gamma].left_;
      int c_prime = this->rectangles_[k].neast_;
      int gamma_prime = this->vertices_[c_prime].west_;
      Tjump_hor(j, RectangulationDirection::right, gamma_prime);
      c = this->rectangles_[j - 1].neast_;
      gamma = this->vertices_[c].west_;
      const int a = this->edges_[gamma].tail_;
      if (this->vertices_[a].type_ == VertexType::right) {
        assert(k == j - 2);
        c_prime = this->rectangles_[j - 2].neast_;
        gamma_prime = this->vertices_[c_prime].west_;
        Sjump(j - 1, RectangulationDirection::right, gamma_prime);
      }
      lock(j - 1, EdgeDir::Ver);
    }
  } else {
    assert(dir == RectangulationDirection::left &&
           this->vertices_[a].type_ == VertexType::right);
    int alpha = this->vertices_[a].south_;
    // N7 - Vertical left jump (T/TS)
    const int k = this->edges_[alpha].left_;
    int b = this->rectangles_[k].swest_;
    int gamma = this->vertices_[b].north_;
    Tjump_ver(j, RectangulationDirection::left, gamma);
    a = this->rectangles_[j].nwest_;
    alpha = this->vertices_[a].east_;
    b = this->edges_[alpha].head_;
    const int beta = this->vertices_[b].east_;
    gamma = this->vertices_[b].north_;
    const int c = this->edges_[beta].head_;
    const int c_prime = this->edges_[gamma].head_;
    if (this->vertices_[c].type_ == VertexType::left &&
        this->vertices_[c_prime].type_ == VertexType::bottom) {
      assert(k == j - 2);
      const int gamma_prime = this->vertices_[a].north_;
      Sjump(j - 1, RectangulationDirection::left, gamma_prime);
    }
    lock(j, EdgeDir::Hor);
  }
}

void Rectangulation::lock(int j, EdgeDir dir) {
  if (dir == EdgeDir::Hor) {
    // L1 - Prepare
    const int a = this->rectangles_[j].neast_;
    const int b = this->rectangles_[j].swest_;
    const int c = this->rectangles_[j].seast_;
    const int alpha = this->vertices_[a].west_;
    const int beta = this->vertices_[b].east_;
    if ((this->vertices_[b].type_ != VertexType::right) ||
        (this->vertices_[c].type_ != VertexType::left) ||
        (this->edges_[beta].head_ != c))
      return;
    // L2 - Lock if necessary
    const int d = this->rectangles_[j + 1].seast_;
    if (this->vertices_[d].type_ == VertexType::top)
      Sjump(j + 1, RectangulationDirection::right, alpha);
  } else {
    assert(dir == EdgeDir::Ver);
    // this code is obtained from the previous case by mirroring along the
    // main diagonal L1 - Prepare
    const int a = this->rectangles_[j].swest_;
    const int b = this->rectangles_[j].neast_;
    const int c = this->rectangles_[j].seast_;
    const int alpha = this->vertices_[a].north_;
    const int beta = this->vertices_[b].south_;
    if ((this->vertices_[b].type_ != VertexType::bottom) ||
        (this->vertices_[c].type_ != VertexType::top) ||
        (this->edges_[beta].tail_ != c))
      return;
    // L2 - Lock if necessary
    const int d = this->rectangles_[j + 1].seast_;
    if (this->vertices_[d].type_ == VertexType::left)
      Sjump(j + 1, RectangulationDirection::left, alpha);
  }
}

void Rectangulation::unlock(int j, RectangulationDirection dir) {
  if (dir == RectangulationDirection::right) {
    // U1 - Prepare
    const int a = this->rectangles_[j].neast_;
    const int b = this->rectangles_[j].seast_;
    const int c = this->rectangles_[j].swest_;
    const int gamma = this->vertices_[c].north_;
    // U2 - Unlock if necessary
    if (this->vertices_[a].type_ == VertexType::bottom &&
        this->vertices_[b].type_ == VertexType::top)
      Sjump(j + 1, RectangulationDirection::left, gamma);
  } else {
    // this code is obtained from the previous case by mirroring along the
    // main diagonal U1 - Prepare
    const int a = this->rectangles_[j].swest_;
    const int b = this->rectangles_[j].seast_;
    const int c = this->rectangles_[j].neast_;
    const int gamma = this->vertices_[c].west_;
    // U2 - Unlock if necessary
    if (this->vertices_[a].type_ == VertexType::right &&
        this->vertices_[b].type_ == VertexType::left)
      Sjump(j + 1, RectangulationDirection::right, gamma);
  }
}

bool Rectangulation::contains_pattern(int j) {
  // check containment of each possible pattern
  for (auto p : this->patterns_) {
    switch (p) {
    case RectangulationPattern::brick_leftright:
      if (contains_brick_leftright(j))
        return true;
      break;
    case RectangulationPattern::brick_rightleft:
      if (contains_brick_rightleft(j))
        return true;
      break;
    case RectangulationPattern::brick_bottomtop:
      if (contains_brick_bottomtop(j))
        return true;
      break;
    case RectangulationPattern::brick_topbottom:
      if (contains_brick_topbottom(j))
        return true;
      break;
    case RectangulationPattern::wmill_clockwise:
      if (contains_wmill_clockwise(j))
        return true;
      break;
    case RectangulationPattern::wmill_counterclockwise:
      if (contains_wmill_counterclockwise(j))
        return true;
      break;
    case RectangulationPattern::H_vertical:
      if (contains_H_vertical(j))
        return true;
      break;
    case RectangulationPattern::H_horizontal:
      if (contains_H_horizontal(j))
        return true;
      break;
    }
  }
  return false;
}

bool Rectangulation::contains_brick_leftright(int j) {
  // C1 - Prepare
  const int a = this->rectangles_[j].nwest_;
  if (this->vertices_[a].type_ == VertexType::bottom)
    return false;
  assert(this->vertices_[a].type_ == VertexType::right);
  // C2 - Check
  const int alpha = this->vertices_[a].south_;
  const int b = this->edges_[alpha].tail_;
  if (this->vertices_[b].type_ == VertexType::left)
    return true;
  else
    return false;
}

bool Rectangulation::contains_brick_rightleft(int j) {
  // C1 - Prepare
  const int a = this->rectangles_[j].nwest_;
  if (this->vertices_[a].type_ == VertexType::bottom)
    return false;
  assert(this->vertices_[a].type_ == VertexType::right);
  // C2 - Check
  const int alpha = this->vertices_[a].north_;
  const int b = this->edges_[alpha].head_;
  if (this->vertices_[b].type_ == VertexType::left)
    return true;
  else
    return false;
}

bool Rectangulation::contains_brick_bottomtop(int j) {
  // this code is obtained by mirroring contains_brick_leftright along the
  // main diagonal C1 - Prepare
  const int a = this->rectangles_[j].nwest_;
  if (this->vertices_[a].type_ == VertexType::right)
    return false;
  assert(this->vertices_[a].type_ == VertexType::bottom);
  // C2 - Check
  const int alpha = this->vertices_[a].east_;
  const int b = this->edges_[alpha].head_;
  if (this->vertices_[b].type_ == VertexType::top)
    return true;
  else
    return false;
}

bool Rectangulation::contains_brick_topbottom(int j) {
  // this code is obtained by mirroring contains_brick_rightleft along the
  // main diagonal C1 - Prepare
  const int a = this->rectangles_[j].nwest_;
  if (this->vertices_[a].type_ == VertexType::right)
    return false;
  assert(this->vertices_[a].type_ == VertexType::bottom);
  // C2 - Check
  const int alpha = this->vertices_[a].west_;
  const int b = this->edges_[alpha].tail_;
  if (this->vertices_[b].type_ == VertexType::top)
    return true;
  else
    return false;
}

bool Rectangulation::contains_wmill_clockwise(int j) {
  // C1 - Prepare
  const int a = this->rectangles_[j].nwest_;
  if (this->vertices_[a].type_ == VertexType::bottom)
    return false;
  assert(this->vertices_[a].type_ == VertexType::right);
  // C2 - Check
  const int alpha = this->vertices_[a].north_;
  const int x = this->edges_[alpha].wall_;
  const int b = this->walls_[x].last_;
  const int beta = this->vertices_[b].east_;
  const int y = this->edges_[beta].wall_;
  const int c = this->walls_[y].last_;
  const int gamma = this->vertices_[c].south_;
  const int z = this->edges_[gamma].wall_;
  const int d = this->walls_[z].first_;
  const int delta = this->vertices_[d].west_;
  if (this->edges_[delta].right_ == j)
    return true;
  else
    return false;
}

bool Rectangulation::contains_wmill_counterclockwise(int j) {
  // this code is obtained by mirroring contains_wmill_clockwise along the
  // main diagonal C1 - Prepare
  const int a = this->rectangles_[j].nwest_;
  if (this->vertices_[a].type_ == VertexType::right)
    return false;
  assert(this->vertices_[a].type_ == VertexType::bottom);
  // C2 - Check
  const int alpha = this->vertices_[a].west_;
  const int x = this->edges_[alpha].wall_;
  const int b = this->walls_[x].first_;
  const int beta = this->vertices_[b].south_;
  const int y = this->edges_[beta].wall_;
  const int c = this->walls_[y].first_;
  const int gamma = this->vertices_[c].east_;
  const int z = this->edges_[gamma].wall_;
  const int d = this->walls_[z].last_;
  const int delta = this->vertices_[d].north_;
  if (this->edges_[delta].right_ == j)
    return true;
  else
    return false;
}

bool Rectangulation::contains_H_vertical(int j) {
  // C1 - Prepare
  const int a = this->rectangles_[j].nwest_;
  if (this->vertices_[a].type_ == VertexType::bottom)
    return false;
  assert(this->vertices_[a].type_ == VertexType::right);
  int b = this->rectangles_[j].swest_;
  // C2 - Go up
  while ((this->vertices_[b].type_ != VertexType::bottom) &&
         (this->vertices_[b].type_ != VertexType::corner)) {
    // C3 - Go left
    int c = b;
    while ((this->vertices_[c].type_ != VertexType::right) &&
           (this->vertices_[c].type_ != VertexType::corner)) {
      if (this->vertices_[c].type_ == VertexType::top) {
        // C4 - Go up
        int d = c;
        while ((d != b) && (this->vertices_[d].type_ != VertexType::bottom) &&
               (this->vertices_[d].type_ != VertexType::corner)) {
          if (this->vertices_[d].type_ == VertexType::left)
            return true;
          const int delta = this->vertices_[d].north_;
          d = this->edges_[delta].head_;
        }
      }
      const int gamma = this->vertices_[c].west_;
      c = this->edges_[gamma].tail_;
    }
    const int beta = this->vertices_[b].north_;
    b = this->edges_[beta].head_;
  }
  return false;
}

bool Rectangulation::contains_H_horizontal(int j) {
  // this code is obtained by mirroring contains_H_vertical along the main
  // diagonal C1 - Prepare
  const int a = this->rectangles_[j].nwest_;
  if (this->vertices_[a].type_ == VertexType::right)
    return false;
  assert(this->vertices_[a].type_ == VertexType::bottom);
  int b = this->rectangles_[j].neast_;
  // C2 - Go left
  while ((this->vertices_[b].type_ != VertexType::right) &&
         (this->vertices_[b].type_ != VertexType::corner)) {
    // C3 - Go up
    int c = b;
    while ((this->vertices_[c].type_ != VertexType::bottom) &&
           (this->vertices_[c].type_ != VertexType::corner)) {
      if (this->vertices_[c].type_ == VertexType::left) {
        // C4 - Go left
        int d = c;
        while ((d != b) && (this->vertices_[d].type_ != VertexType::right) &&
               (this->vertices_[d].type_ != VertexType::corner)) {
          if (this->vertices_[d].type_ == VertexType::top)
            return true;
          const int delta = this->vertices_[d].west_;
          d = this->edges_[delta].tail_;
        }
      }
      const int gamma = this->vertices_[c].north_;
      c = this->edges_[gamma].head_;
    }
    const int beta = this->vertices_[b].west_;
    b = this->edges_[beta].tail_;
  }
  return false;
}
