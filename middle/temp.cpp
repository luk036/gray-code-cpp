#include <iostream>
#include <vector>

// A class to represent and manipulate a vertex of the cube of odd dimension
// 2n+1 for use within the Hamilton cycle algorithm.

class Vertex {
public:
  // constructor
  explicit Vertex(const std::vector<int> &x);

  const std::vector<int> &get_bits() const { return bits_; }
  // allow bit vector access directly via the vertex object
  int &operator[](int i) { return this->bits_[i]; }
  const int &operator[](int i) const { return this->bits_[i]; }
  // the size of the bit vector representing the vertex (=2n+1)
  int size() { return this->bits_.size(); }
  const int size() const { return this->bits_.size(); }

  // reverse and invert bitstring
  void rev_inv();

  // compute whether vertex is first or last vertex on a path
  bool is_first_vertex() const;
  bool is_last_vertex() const;

  // For the current vertex, compute the corresponding first/last vertex
  // on the unflipped path. Return the distance between the two vertices
  // along the path.
  int to_first_vertex();
  int to_last_vertex();

  // Given a first vertex of a path, compute the sequence of bit positions
  // to be flipped to follow the corresponding path all the way to the
  // last vertex. Flipping the bit positions one after the other, we
  // follow a path from P_n\circ 0 from its first vertex to its last
  // vertex in the graph H_n\circ 0. If flip == false, then follow the
  // unflipped path, otherwise the flipped path. This is exactly the
  // recursion rule sigma() defined in the paper. The current vertex is
  // modified temporarily, which is the only reason that this function is
  // not qualified as const.
  void compute_flip_seq_0(std::vector<int> &seq, bool flip);
  // Given a last vertex of a path, compute the sequence of bit positions
  // to be flipped to follow the corresponding path all the way to the
  // first vertex. The resulting positions are already transformed
  // according to the isomorphism \overline{\rev}
  // (complementation+reversal), so flipping the bit positions one after
  // the other, we follow a path from \overline{\rev}(P_n)\circ 1 from its
  // last vertex to its first vertex in the graph
  // \overline{\rev}(H_n)\circ 1. This recursion rule is derived from
  // sigma() defined in the paper. In the algorithm HamCycle() described
  // in the paper, we instead apply sigma() and transform the resulting
  // flip sequence afterwards. However, computing the recursion directly
  // via this modified recursion is faster in practice.
  void compute_flip_seq_1(std::vector<int> &seq) const;

private:
  // the bitstring representation of the vertex
  std::vector<int> bits_;

  // #### auxiliary functions ####

  // reverse and invert bitstring in the range [left,right]
  void rev_inv(int left, int right);
  // auxiliary function to compute the smallest index b
  // for which the Dyck path represented by this->bits_[a,...,b]
  // returns to the abscissa
  int first_touchdown(int a) const;
  // auxiliary function to compute the smallest index b
  // for which the Dyck path represented by this->bits_[0,...,b]
  // moves below the abscissa
  int first_dive() const;
  // auxiliary function to compute the positions of all downsteps and
  // upsteps below the abscissa and above the abscissa sorted by
  // increasing depth/height of those steps
  void steps_height(std::vector<std::vector<int>> &usteps_neg,
                    std::vector<std::vector<int>> &usteps_pos,
                    std::vector<std::vector<int>> &dsteps_neg,
                    std::vector<std::vector<int>> &dsteps_pos) const;

  // interpreting the bitstring representation as a Dyck path
  // (0=downstep, 1=upstep), count the number of downsteps below the
  // abscissa
  int count_flaws() const;

  // count the number of 1-bits in this->bits_ (when ignoring the last
  // bit)
  int count_ones() const;

  // Auxiliary recursive version of the function compute_flip_seq_0().
  // Compute the flip sequence for the given first vertex, interpreted
  // as a Dyck path, in the range [left,...,right]. The value idx is the
  // index into the array seq where to write the next bit positions.
  // The auxiliary array next_step is an array of bidirectional pointers
  // below the hills of the Dyck path between pairs of an upstep
  // and downstep that 'see each other' (this allows the canonical
  // decomposition to be performed in constant time).
  // This is the recursion rule sigma() defined in the paper.
  // The array next_step is allocated on the stack for speed reasons.
  void compute_flip_seq_0_rec(std::vector<int> &seq, int &idx, int left,
                              int right, int *next_step) const;
  // Auxiliary recursive version of the function compute_flip_seq_1().
  void compute_flip_seq_1_rec(std::vector<int> &seq, int &idx, int left,
                              int right, int *next_step) const;
  // Compute auxiliary bidirectional pointers for the Dyck path
  // in the range [a,b] below the hills of the Dyck path between
  // pairs of an upstep and downstep that 'see each other'.
  // We allocate this array on the stack for speed reasons.
  void aux_pointers(int a, int b, int *next_step) const;
};

// allow comparison between vertices (by comparing the corresponding bit
// vectors)
inline bool operator==(const Vertex &lhs, const Vertex &rhs) {
  return lhs.get_bits() == rhs.get_bits();
}
inline bool operator!=(const Vertex &lhs, const Vertex &rhs) {
  return !operator==(lhs, rhs);
}
std::ostream &operator<<(std::ostream &os, const Vertex &v);

// compare two bitstrings of the same length (arrays of 0s and 1s)
// lexicographically
bool bitstrings_less_than(int *x, int *y, int length);
bool bitstrings_equal(int *x, int *y, int length);

#include <algorithm>
#include <cassert>
#include <vector>

Vertex::Vertex(const std::vector<int> &x) : bits_(x) {
  assert(x.size() % 2 == 1);
  assert(x.size() >= 3);
}

void Vertex::rev_inv() {
  rev_inv(0, this->bits_.size() - 2); // ignore the last bit
}

void Vertex::rev_inv(int left, int right) {
  for (int i = left; i <= right; ++i) {
    this->bits_[i] = 1 - this->bits_[i];
  }
  std::reverse(this->bits_.begin() + left, this->bits_.begin() + right + 1);
}

int Vertex::first_touchdown(int a) const {
  int height = 0;
  for (int i = a; i < this->bits_.size() - 1; ++i) { // ignore the last bit
    height += (2 * this->bits_[i] - 1);
    if (height == 0) {
      return i;
    }
  }
  return -1; // error
}

int Vertex::first_dive() const {
  int height = 0;
  for (int i = 0; i < this->bits_.size() - 1; ++i) { // ignore the last bit
    height += (2 * this->bits_[i] - 1);
    if (height == -1) {
      return i;
    }
  }
  return -1; // error
}

void Vertex::steps_height(std::vector<std::vector<int>> &usteps_neg,
                          std::vector<std::vector<int>> &usteps_pos,
                          std::vector<std::vector<int>> &dsteps_neg,
                          std::vector<std::vector<int>> &dsteps_pos) const {
  usteps_neg.resize(0);
  usteps_pos.resize(0);
  dsteps_neg.resize(0);
  dsteps_pos.resize(0);
  int height = 0;
  int min_height = 0;
  int max_height = 0;
  for (int i = 0; i < this->bits_.size() - 1; ++i) { // ignore the last bit
    if ((this->bits_[i] == 0) && (height <= 0)) {
      if (height == min_height) {
        usteps_neg.push_back(std::vector<int>());
        dsteps_neg.push_back(std::vector<int>());
      }
      dsteps_neg[-height].push_back(i);
    }
    if ((this->bits_[i] == 1) && (height >= 0)) {
      if (height == max_height) {
        usteps_pos.push_back(std::vector<int>());
        dsteps_pos.push_back(std::vector<int>());
      }
      usteps_pos[height].push_back(i);
    }
    height += (2 * this->bits_[i] - 1); // update height
    min_height = std::min(height, min_height);
    max_height = std::max(height, max_height);
    if ((this->bits_[i] == 0) && (height >= 0)) {
      dsteps_pos[height].push_back(i);
      assert(dsteps_pos[height].size() == usteps_pos[height].size());
    }
    if ((this->bits_[i] == 1) && (height <= 0)) {
      usteps_neg[-height].push_back(i);
      assert(usteps_neg[-height].size() == dsteps_neg[-height].size());
    }
  }
  assert(usteps_neg.size() == dsteps_neg.size());
}

int Vertex::count_flaws() const {
  int c = 0;
  int height = 0;
  for (int i = 0; i < this->bits_.size() - 1; ++i) { // ignore the last bit
    if ((height <= 0) && (this->bits_[i] == 0)) {
      ++c;
    }
    height += (2 * this->bits_[i] - 1);
  }
  return c;
}

int Vertex::count_ones() const {
  int c = 0;
  for (int i = 0; i < this->bits_.size() - 1; ++i) { // ignore last bit
    if (this->bits_[i] == 1) {
      ++c;
    }
  }
  return c;
}

bool Vertex::is_first_vertex() const {
  return ((count_flaws() == 0) && (count_ones() == this->bits_.size() / 2));
}

bool Vertex::is_last_vertex() const {
  return ((count_flaws() == 1) && (count_ones() == this->bits_.size() / 2));
}

int Vertex::to_first_vertex() {
  if (is_last_vertex()) {
    // This case is encountered during the Hamilton cycle
    // computation.
    const int b = first_dive();
    // shift this->bits_[0,...,b-1] to the right
    std::copy(this->bits_.begin(), this->bits_.begin() + b,
              this->bits_.begin() + 1);
    // reset boundary steps around shifted subpath
    this->bits_[0] = 1;
    this->bits_[b + 1] = 0;
    return (2 * b + 2);
  } else {
    // This case is encountered only during initialization.
    // The rule how to move from an intermediate path vertex back to
    // the first vertex is  described in detail in the paper.
    std::vector<std::vector<int>> usteps_neg, dsteps_neg, usteps_pos,
        dsteps_pos;
    steps_height(usteps_neg, usteps_pos, dsteps_neg, dsteps_pos);
    bool min_zero = (usteps_neg.size() ==
                     0); // minimum of lattice path lies on the abscissa y=0
    bool unique_min;     // minimum of lattice path is unique
    unique_min = (min_zero ? (usteps_pos.front().size() == 1)
                           : (usteps_neg.back().size() == 1));
    bool middle_level = (2 * count_ones() + 1 ==
                         this->bits_.size()); // vertex has same number of 0s
                                              // and 1s (ignoring the last bit)
    int to;
    if ((!unique_min && middle_level) || (unique_min && !middle_level)) {
      // take position of first upstep starting at minimum
      // height
      to = (min_zero ? usteps_pos.front().front() : usteps_neg.back().front()) -
           1;
    } else {
      // take position of last upstep starting at minimum
      // height
      to =
          (min_zero ? usteps_pos.front().back() : usteps_neg.back().back()) - 1;
    }
    // shift this->bits_[0,to] to the right by one
    std::copy(this->bits_.begin(), this->bits_.begin() + to + 1,
              this->bits_.begin() + 1);
    this->bits_[0] = 1;
    // flip downsteps in the left half to upsteps
    for (int d = 0;
         d < ((int)dsteps_neg.size()) - ((unique_min && middle_level) ? 1 : 0);
         ++d) {
      this->bits_[dsteps_neg[d].front() + 1] = 1;
    }
    // flip upsteps in the right half to downsteps
    for (int d = 0;
         d < ((int)usteps_neg.size()) - ((unique_min && !middle_level) ? 1 : 0);
         ++d) {
      this->bits_[usteps_neg[d].back()] = 0;
    }
    if (!middle_level) {
      // flip further upsteps above the line y=0 to downsteps
      for (int d = ((min_zero && unique_min) ? 1 : 0); d <= 1; ++d) {
        this->bits_[usteps_pos[d].back()] = 0;
      }
    }
    return 2 * (to + 1) + (middle_level ? 0 : 1);
  }
}

int Vertex::to_last_vertex() {
  int d = 0;
  if (!is_first_vertex()) {
    d = -to_first_vertex();
  }
  assert(is_first_vertex());

  const int b = first_touchdown(0);
  // shift this->bits_[1,...,b-1] to the left
  std::copy(this->bits_.begin() + 1, this->bits_.begin() + b,
            this->bits_.begin());
  // reset two steps that form a valley to the right of the shifted
  // subpath
  this->bits_[b - 1] = 0;
  this->bits_[b] = 1;
  d += 2 * (b - 1) + 2;
  return d;
}

void Vertex::compute_flip_seq_0(std::vector<int> &seq, bool flip) {
  assert(is_first_vertex());

  if (!flip) {
    // The length of the flip sequence is 2*(b-1)+2
    // where b is the index of the first time the Dyck
    // path corresponding to the vertex returns to the
    // abscissa.
    const int b = first_touchdown(0);
    const int length = 2 * (b - 1) + 2;
    seq.resize(length, 0);

    int next_step[b + 1];
    aux_pointers(0, b, next_step);

    int idx = 0;
    seq[idx++] = b;
    seq[idx++] = 0;
    compute_flip_seq_0_rec(seq, idx, 1, b - 1, next_step);
    return;
  }

  assert(flip);
  // From the length of the bitstring representation we can recover in
  // which of the two cases tau_image/tau_preimage we are.
  assert(this->bits_[0] == 1);
  if (this->bits_[1] == 1) {
    // create path of length 2
    assert(this->bits_[2] == 0);
    seq.resize(2);
    seq[0] = 2;
    seq[1] = 0;
  } else { // tau_image
    // We artificially change the bitstring into the other starting
    // vertex of the flippable pair, as this is the path we will be
    // modifying.
    this->bits_[1] = 1;
    this->bits_[2] = 0;

    const int b = first_touchdown(0);
    const int length = 2 * (b - 1) + 2;
    seq.resize(length, 0);

    int next_step[b + 1];
    aux_pointers(0, b, next_step);

    int idx = 0;
    seq[idx++] = b;
    seq[idx++] = 0;
    compute_flip_seq_0_rec(seq, idx, 1, b - 1, next_step);

    // revert artificial bitstring modification
    this->bits_[1] = 0;
    this->bits_[2] = 1;

    // modify flip sequence
    assert((seq[0] == b) && (seq[1] == 0) && (seq[2] == 2) && (seq[3] == 1) &&
           (seq[4] == 0) && (seq[5] == 2));
    seq[0] = b;
    seq[1] = 0;
    seq[2] = 1;
    seq[3] = 2;
    seq[4] = 0;
    seq[5] = 1;
  }
}

void Vertex::compute_flip_seq_0_rec(std::vector<int> &seq, int &idx, int left,
                                    int right, int *next_step) const {
  const int length =
      right - left + 1; // total length of the subpath under consideration
  if (length <= 0) {    // base case of the recursion
    return;
  }
  assert((this->bits_[left] == 1) && (this->bits_[right] == 0) &&
         (length % 2 == 0));

  const int m = next_step[left];
  assert((m <= right) && (this->bits_[m] == 0));
  seq[idx++] = m;
  seq[idx++] = left;
  // descend recursively into the left subpath
  compute_flip_seq_0_rec(seq, idx, left + 1, m - 1, next_step);
  seq[idx++] = left - 1;
  seq[idx++] = m;
  // descend recursively into the right subpath
  compute_flip_seq_0_rec(seq, idx, m + 1, right, next_step);
}

void Vertex::compute_flip_seq_1(std::vector<int> &seq) const {
  assert(is_last_vertex());

  const int b = first_dive();
  const int length = 2 * ((this->bits_.size() - 2) - (b + 2) + 1) + 2;
  seq.resize(length, 0);

  int next_step[this->bits_.size() - 1];
  aux_pointers(b + 2, this->bits_.size() - 2, next_step);

  int idx = 0;
  seq[idx++] = b + 1;
  compute_flip_seq_1_rec(seq, idx, b + 2, this->bits_.size() - 2, next_step);
  seq[idx++] = b;
}

void Vertex::compute_flip_seq_1_rec(std::vector<int> &seq, int &idx, int left,
                                    int right, int *next_step) const {
  const int length =
      right - left + 1; // total length of the subpath under consideration
  if (length <= 0) {    // base case of the recursion
    return;
  }
  assert((this->bits_[left] == 1) && (this->bits_[right] == 0) &&
         (length % 2 == 0));

  // index of the first return to the same height
  const int m = next_step[left];
  // length of the two subpaths u and v in the canonical decomposition

  seq[idx++] = m;
  seq[idx++] = left;
  // descend recursively into the left subpath
  compute_flip_seq_1_rec(seq, idx, left + 1, m - 1, next_step);
  seq[idx++] = left - 1;
  seq[idx++] = m;
  // descend recursively into the right subpath
  compute_flip_seq_1_rec(seq, idx, m + 1, right, next_step);
}

void Vertex::aux_pointers(int a, int b, int *next_step) const {
  assert((a == b + 1) || ((this->bits_[a] == 1) && (this->bits_[b] == 0)));
  // The array left_ustep_height[h] contains the index of the last upstep
  // starting at height h that has been encountered when moving
  // along the Dyck path from left to right.
  // We allocate this array on the stack for speed reasons.
  int left_ustep_height[b - a + 1]; // only the first (b-a)/2+1 many
                                    // entries should be needed
  int height = 0;
  for (int i = a; i <= b; ++i) {
    if (this->bits_[i] == 0) { // downstep (0-bit)
      assert(height >= 1);
      const int left = left_ustep_height[height - 1];
      assert((left >= 0) && (left < i));
      next_step[left] = i;
      next_step[i] = left;
    } else { // upstep (1-bit)
      assert(height >= 0);
      left_ustep_height[height] = i;
    }
    height += (2 * this->bits_[i] - 1);
  }
  assert(height == 0);
}

std::ostream &operator<<(std::ostream &os, const Vertex &v) {
  // bitstrings can be printed unambiguously without separation characters
  for (int i = 0; i < v.size(); ++i) {
    os << v[i];
  }
  return os;
}

bool bitstrings_less_than(int *x, int *y, int length) {
  for (int i = 0; i < length; ++i) {
    if (x[i] < y[i]) {
      return true;
    } else if (x[i] > y[i]) {
      return false;
    }
  }
  // bitstrings are the same
  return false;
}

bool bitstrings_equal(int *x, int *y, int length) {
  for (int i = 0; i < length; ++i) {
    if (x[i] != y[i]) {
      return false;
    }
  }
  // bitstrings are the same
  return true;
}

#include <list>
#include <vector>

// A class to represent and manipulate an ordered rooted tree in doubly linked
// adjacency list representation.

class Tree {
public:
  // Upon construction, we build the tree from the bitstring
  // representation of the given Vertex object. The first even number of
  // bits of the bitstring are interpreted as a lattice path that always
  // stays above the y-axis (=a Dyck path). Every 1-bit is interpreted as
  // attaching an edge/child to the current vertex of the tree, every
  // 0-bit is interpreted as going back to the parent vertex (without
  // adding any edges to the tree). For example, the bitstring
  // x=110111000010 corresponds to the following tree:
  //   /\
  //  /\
  //    \
  //     \
  //
  explicit Tree(const Vertex &x);

  // Decide if the tree is an end node of a directed edge in the
  // canonic spanning tree in the auxiliary graph defined in the paper
  // that defines a Hamilton cycle in the middle levels graph.
  // This function simultaneously applies tau() or tau_inverse()
  // to the tree.
  bool flip_tree();

  // rotate the tree by making the leftmost child of the root the new root
  void rotate();

  // Computes the bitstring representation of the tree.
  // The output array x should be allocated on the stack for speed
  // reasons.
  void to_bitstring(int x[]) const;

private:
  int num_vertices_;
  int root_;
  // Doubly linked adjacency lists where this->children_[i] is a list of
  // the children of vertex i starting with the leftmost child and ending
  // with the rightmost child. The class assumes that the vertices of the
  // tree are numbers from 0 up to num_vertices_-1.
  std::vector<std::list<int>> children_;
  // The entry parent[i] is the parent of vertex i. The value of the
  // parent of this->root_ is undefined.
  std::vector<int> parent_;

  // #### auxiliary functions ####

  // the degree of a vertex
  int deg(int u) const;
  // the number of children of vertex u
  int num_children(int u) const;

  // Return the i-th child of vertex u when going through the children
  // from left to right. The value i is an integer in the interval [0,
  // num_children(u)[.
  int ith_child(int u, int i) const;

  // decide if the tree is a preimage of the mapping tau defined in the
  // paper
  bool is_tau_preimage() const;
  // decide if the tree is an image of the mapping tau defined in the
  // paper
  bool is_tau_image() const;

  // Apply the mapping tau defined in the paper.
  // The function assumes that this->is_tau_preimage == true.
  void tau();
  // Compute the inverse of the mapping tau defined in the paper.
  // The function assumes that this->is_tau_image == true.
  void tau_inverse();

  // Remove the given leaf of the tree and attach it to the vertex
  // new_parent at position pos in new_parent's list of children. The
  // position must be in the interval [0,...,num_children[new_parent]].
  void move_leaf(int leaf, int new_parent, int pos);

  // rotate until the given vertex u is the root
  void rotate_to_vertex(int u);
  // make the leftmost subtree of the root the rightmost subtree
  void rotate_children();
  // apply the previous rotation k times
  void rotate_children(int k);

  // Compute a canonical way to root the tree, i.e., the resulting
  // rooted tree will be the same regardless from which rooted
  // version the function is called (apart from vertex names, so
  // the resulting tree will have the same bitstring representation).
  void root_canonically();

  // Compute the center vertices of the tree.
  // If there is only one center, then it is returned in c1 and c2 is set
  // to -1. If there are two centers, then they are returned in c1 and c2.
  void compute_center(int &c1, int &c2) const;

  // Decide if the the edge (tree,tau(tree)) in the auxiliary graph
  // defined in the paper belongs to the canonic spanning tree that
  // defines a Hamilton cycle in the middle levels graph. See the paper
  // for a detailed definition of this function. The function assumes that
  // is_tau_preimage(tree) == true.
  bool is_flip_tree_tau(); // We promise that we do not modify the tree in
                           // this function. For speed reasons however we
                           // modify it temporarily and restore the state
                           // at the end of the function rather than
                           // making a copy of ourselves.

  // checks if the tree is a star, rooted at the center of a leaf
  bool is_star() const;
  // checks if tree has the form 1(10)^k0(10)^l with 1<=l<k
  bool is_light_dumbbell() const;
  // check if the vertex u is a thin leaf, i.e., whether it has degree 1
  // and its neighbor degree 2
  bool is_thin_leaf(int u) const;
  // checks if the tree has a thin leaf somewhere
  bool has_thin_leaf() const;
  // For the subtree rooted at the given vertex u, count the number
  // of pending edges starting at u from left to right until the first
  // non-pending edge is encountered.
  int count_pending_edges(int u) const;

  // Auxiliary recursive function to compute the bitstring representation.
  // The output array x should be allocated on the stack for speed
  // reasons.
  void to_bitstring_rec(int x[], int u, int &pos) const;

  // Booth's algorithm to compute lexicographically smallest
  // rotation of the given string/array. Return value is the
  // index where the lexicographically smallest string (viewed
  // cyclically) starts. The running time is linear in the length
  // of the input array.
  // The input array x should be allocated on the stack for speed reasons.
  // The reference to the original paper is:
  // [K. Booth, Lexicographically least circular substrings,
  //  Inf. Proc. Letters, 10 (4-5): 240–242]
  int min_string_rotation(int x[], int length);
};

#include <algorithm>
#include <cassert>
#include <cstring>
#include <list>
#include <vector>

Tree::Tree(const Vertex &x) {
  std::vector<int> xv = x.get_bits(); // extract bitstring representation
  assert(xv.size() % 2 == 1);

  this->root_ = 0;
  this->num_vertices_ = (xv.size() - 1) / 2 + 1;
  this->children_.resize(this->num_vertices_);
  this->parent_.resize(num_vertices_, 0);

  int u = this->root_; // the current vertex
  int n = 1;           // the number of vertices created so far
  int height = 0;
  for (int i = 0; i < xv.size() - 1; ++i) {
    if (x[i] == 1) {
      // create a new child vertex and add an edge from the
      // current vertex leading to the new child
      this->children_[u].push_back(n);
      this->parent_[n] = u;
      u = n;
      ++n;
    } else {
      // move back to the parent of the current vertex
      u = this->parent_[u];
    }
    height += (2 * x[i] - 1);
    assert(height >= 0);
  }
  assert(n == this->num_vertices_);
}

int Tree::deg(int u) const {
  assert((0 <= u) && (u < this->num_vertices_));
  if (u == this->root_) {
    return this->children_[u].size();
  } else {
    return this->children_[u].size() + 1;
  }
}

int Tree::num_children(int u) const {
  assert((0 <= u) && (u < this->num_vertices_));
  return this->children_[u].size();
}

int Tree::ith_child(int u, int i) const {
  assert((0 <= u) && (u < this->num_vertices_));
  assert((0 <= i) && (i < num_children(u)));
  // move to i-th entry from the beginning
  std::list<int>::const_iterator it = this->children_[u].begin();
  std::advance(it, i);
  return *it;
}

// The tree has to look like this (x and y are arbitrary subtrees):
//   /x
//  /y
bool Tree::is_tau_preimage() const {
  if (this->num_vertices_ < 3) {
    return false;
  }
  // u is the leftmost child of the root
  const int u = ith_child(root_, 0);
  if (num_children(u) == 0) {
    return false;
  }
  // v is the leftmost child of u
  const int v = ith_child(u, 0);
  if (num_children(v) != 0) {
    return false;
  }
  return true;
}

// The tree has to look like this (x and y are arbitrary subtrees):
//   /|x
//    y
bool Tree::is_tau_image() const {
  if ((this->num_vertices_ < 3) || (num_children(this->root_) < 2) ||
      (num_children(ith_child(root_, 0)) > 0)) {
    return false;
  }
  return true;
}

void Tree::tau() {
  assert(is_tau_preimage());
  const int u = ith_child(root_, 0);
  const int v = ith_child(u, 0);
  move_leaf(v, root_, 0);
}

void Tree::tau_inverse() {
  assert(is_tau_image());
  const int v = ith_child(root_, 0);
  const int u = ith_child(root_, 1);
  move_leaf(v, u, 0);
}

void Tree::move_leaf(int leaf, int new_parent, int pos) {
  assert((0 <= leaf) && (leaf < this->num_vertices_));
  assert((0 <= new_parent) && (new_parent < this->num_vertices_));
  assert((0 <= pos) && (pos <= this->children_[new_parent].size()));
  assert(num_children(leaf) == 0);
  const int old_parent = this->parent_[leaf];
  // search through the children of the current parent
  for (std::list<int>::iterator it = this->children_[old_parent].begin();
       it != this->children_[old_parent].end(); ++it) {
    if (*it == leaf) { // remove this child
      this->children_[old_parent].erase(it);
      break;
    }
  }
  // add the leaf below the new parent vertex
  std::list<int>::iterator it = this->children_[new_parent].begin();
  std::advance(it, pos);
  this->children_[new_parent].insert(it, leaf);
  this->parent_[leaf] = new_parent;
}

void Tree::rotate() {
  assert(this->num_vertices_ >= 2);
  const int u = ith_child(root_, 0);
  this->parent_[root_] = u;
  // move first entry from root_'s list of children to end of u's list of
  // children (this is much faster than pop_front() and push_back())
  this->children_[u].splice(this->children_[u].end(), this->children_[root_],
                            this->children_[root_].begin());
  this->children_[u].back() = root_;
  this->root_ = u;
}

void Tree::rotate_to_vertex(int u) {
  while (this->root_ != u) {
    rotate();
  }
}

void Tree::rotate_children() { rotate_children(1); }

void Tree::rotate_children(int k) {
  std::list<int>::iterator it = this->children_[root_].begin();
  std::advance(it, k);
  std::rotate(this->children_[root_].begin(), it, this->children_[root_].end());
}

bool Tree::flip_tree() {
  if (is_tau_preimage() && is_flip_tree_tau()) {
    tau();
    return true;
  } else if (is_tau_image()) {
    tau_inverse();
    if (is_flip_tree_tau()) {
      return true;
    }
    tau(); // undo tau^{-1}
  }
  return false;
}

void Tree::root_canonically() {
  int c1, c2; // center vertices
  compute_center(c1, c2);
  if (c2 != -1) { // centers are different
    // compute bitstring representation x1 when rooting
    // the tree at c1 such that c2 is the leftmost child
    const int num_bits = 2 * (this->num_vertices_ - 1);
    int x1[num_bits];
    int x2[num_bits];
    rotate_to_vertex(c1);
    while (ith_child(root_, 0) != c2) {
      rotate_children();
    }
    to_bitstring(x1);

    // compute bitstring representation x2 when rooting
    // the tree at c2 such that c1 is the leftmost child
    rotate();
    rotate_children(num_children(this->root_) - 1);
    assert((this->root_ == c2) && (ith_child(root_, 0) == c1));
    to_bitstring(x2);

    // The canonical rooting of the tree is the one corresponding
    // to the lexicographically smallest bitstring representation
    if (bitstrings_less_than(x1, x2, num_bits)) {
      rotate();
      rotate_children(num_children(this->root_) - 1);
      assert((this->root_ == c1) && (ith_child(root_, 0) == c2));
    }
  } else { // centers are the same
    // root at the center and compute bitstring representation
    rotate_to_vertex(c1);
    const int num_bits = 2 * (this->num_vertices_ - 1);
    int x[num_bits];
    to_bitstring(x);

    // compute segments of the bitstring representation
    // belonging to the different subtrees
    int subtree_count[num_bits];
    int c = 0;
    int depth = 0;
    for (int i = 0; i < num_bits; ++i) {
      if (x[i] == 1) {
        ++depth;
      } else { // x[i] == 0
        --depth;
      }
      subtree_count[i] = c;
      if (depth == 0) {
        ++c;
      }
    }
    assert(c == num_children(this->root_));

    // Compute the lexicographically smallest rotation of the given
    // string. Note that this function returns slightly different
    // results than the function used in the paper which adds
    // additional -1s to between subtrees. We do not add these -1s,
    // but still obtain a canonically rooted tree.
    const int k = min_string_rotation(x, num_bits);

    // rotate children accordingly
    rotate_children(subtree_count[k]);
  }
}

void Tree::compute_center(int &c1, int &c2) const {
  // set vertex degrees and store leaves
  std::vector<int> degs(num_vertices_, 0);
  std::vector<int> leaves(num_vertices_,
                          0); // for sure this many entries will be enough
  int num_leaves = 0;
  for (int i = 0; i < num_vertices_; ++i) {
    degs[i] = deg(i);
    if (degs[i] == 1) {
      leaves[num_leaves++] = i;
    }
  }

  int num_vertices_remaining = num_vertices_;
  int num_new_leaves = 0;
  // cut away leaves in rounds until only 1 or 2 vertices are left,
  // they form the center vertices
  while (num_vertices_remaining > 2) {
    // remove leaves
    for (int i = 0; i < num_leaves; ++i) {
      const int u = leaves[i];
      for (std::list<int>::const_iterator it = this->children_[u].begin();
           it != this->children_[u].end(); ++it) {
        --degs[*it];
        if (degs[*it] == 1) {             // remember leaves for the next round
          leaves[num_new_leaves++] = *it; // we can fill the leaves for
                                          // the next round into the same
                                          // vector from the beginning,
                                          // as the number of leaves
                                          // decreases in every round
        }
      }
      if (u != this->root_) {
        --degs[this->parent_[u]];
        if (degs[this->parent_[u]] == 1) { // remember leaves for the next round
          leaves[num_new_leaves++] = this->parent_[u];
        }
      }
    }
    num_vertices_remaining -= num_leaves;
    num_leaves = num_new_leaves;
    num_new_leaves = 0;
  }
  assert((num_leaves >= 1) && (num_leaves <= 2));

  if (num_leaves == 1) {
    c1 = leaves[0];
    c2 = -1;
  } else {
    c1 = leaves[0];
    c2 = leaves[1];
  }
}

bool Tree::is_flip_tree_tau() {
  if (is_star()) {
    return false;
  }

  // remember root and its leftmost child
  const int r = root_;
  const int u = ith_child(root_, 0);

  const int num_bits = 2 * (this->num_vertices_ - 1);
  int this_bitstring[num_bits];
  int canon_bitstring[num_bits];

  int v = ith_child(root_, 0);
  if ((num_children(v) == 1) && (num_children(ith_child(v, 0)) == 0)) {
    // tree has the form 1100...
    // compute my bitstring representation
    to_bitstring(this_bitstring);
    // compute canonically rooted version of myself
    root_canonically();
    // rotate until tree has the form 1100... again
    v = ith_child(root_, 0);
    while ((num_children(v) != 1) || (num_children(ith_child(v, 0)) != 0)) {
      rotate();
      v = ith_child(root_, 0);
    };
  } else {
    if (has_thin_leaf()) { // tree should not have thin leaves
      return false;
    }
    v = ith_child(root_, 0);
    int c = count_pending_edges(v);
    if ((c < num_children(v)) || (c < 2) || (is_light_dumbbell())) {
      return false;
    }
    // tree has the form 1(10)^k0... with k>=2
    // compute my bitstring representation
    to_bitstring(this_bitstring);
    // compute canonically rooted version of myself
    root_canonically();
    // rotate until tree has the form 1(10)^k0... with k>=2 again
    v = ith_child(root_, 0);
    c = count_pending_edges(v);
    while ((c < num_children(v)) || (c < 2)) {
      rotate();
      rotate_children(c);
      v = ith_child(root_, 0);
      c = count_pending_edges(v);
    }
  }

  // compute bitstring representation of canonically rooted version
  to_bitstring(canon_bitstring);

  // restore tree to original state
  rotate_to_vertex(r);
  while (ith_child(root_, 0) != u) {
    rotate_children();
  }

  // compare bitstrings
  if (bitstrings_equal(this_bitstring, canon_bitstring, num_bits)) {
    return true;
  } else {
    return false;
  }
}

bool Tree::is_star() const {
  if ((this->num_vertices_ <= 3) ||
      (deg(this->root_) == this->num_vertices_ - 1) ||
      (deg(ith_child(root_, 0)) == this->num_vertices_ - 1)) {
    return true;
  }
  return false;
}

bool Tree::is_light_dumbbell() const {
  if (this->num_vertices_ < 5) {
    return false;
  }
  const int u = ith_child(root_, 0);
  const int k = num_children(u);
  const int l = num_children(this->root_) - 1;
  if ((k + l + 1 < this->num_vertices_ - 1) || (k <= l)) {
    return false;
  } else {
    return true;
  }
}

bool Tree::is_thin_leaf(int u) const {
  if (deg(u) > 1)
    return false;
  if (((u == this->root_) && (deg(ith_child(u, 0)) == 2)) ||
      ((u != this->root_) && (deg(this->parent_[u]) == 2))) {
    return true;
  } else {
    return false;
  }
}

bool Tree::has_thin_leaf() const {
  for (int i = 0; i < num_vertices_; ++i) {
    if (is_thin_leaf(i)) {
      return true;
    }
  }
  return false;
}

int Tree::count_pending_edges(int u) const {
  int c = 0;
  for (int i = 0; i < num_children(u); ++i) {
    const int v = ith_child(u, i);
    if (num_children(v) == 0) {
      ++c;
    } else {
      return c;
    }
  }
  return c;
}

void Tree::to_bitstring(int x[]) const {
  int pos = 0;
  to_bitstring_rec(x, root_, pos);
}

void Tree::to_bitstring_rec(int x[], int u, int &pos) const {
  if (num_children(u) == 0) {
    return;
  } else {
    for (std::list<int>::const_iterator it = this->children_[u].begin();
         it != this->children_[u].end(); ++it) {
      x[pos++] = 1;
      to_bitstring_rec(x, *it, pos);
      x[pos++] = 0;
    }
  }
}

// The source code for this implementation of Booth's
// algorithm was copied verbatim from the following Wikipedia
// site: "Lexicographically minimal string rotation"
int Tree::min_string_rotation(int x[], int length) {
  // concatenate array with itself to avoid modular arithmetic
  int xx[2 * length];
  std::memcpy(xx, x, sizeof(int) * length);
  std::memcpy(xx + length, x, sizeof(int) * length);
  // failure function
  std::vector<int> fail(2 * length, -1);
  int k = 0; // lexicographically smallest starting position found so far
  for (int j = 1; j < 2 * length; ++j) {
    const int xj = xx[j];
    int i = fail[j - k - 1];
    while ((i != -1) && (xj != xx[k + i + 1])) {
      if (xj < xx[k + i + 1]) {
        k = j - i - 1;
      }
      i = fail[i];
    }
    if (xj != xx[k + i + 1]) {
      if (xj < xx[k]) {
        k = j;
      }
      fail[j - k] = -1;
    } else {
      fail[j - k] = i + 1;
    }
  }
  return k;
}

#include <vector>

// type of user-defined visit function
typedef void (*visit_f_t)(const std::vector<int> &y, int i);

class HamCycle {
public:
  // Starting from a vertex x of the middle levels graph G_n (the
  // parameter n is implicit in the number 2n+1 of bits of x), compute the
  // next limit vertices on a Hamilton cycle (if the vertex x is reached
  // again, the computation is terminated prematurely). A value limit < 0
  // means computing until we are back at the starting vertex. If
  // store_vertices == true, then the vector of flipped bit positions is
  // stored along the way (otherwise not).
  explicit HamCycle(const Vertex &x, long long limit, visit_f_t visit_f);
  long long get_length() const { return length_; }

private:
  Vertex x_;          // starting vertex
  Vertex y_;          // current vertex
  long long limit_;   // the number of vertices to be visited
  visit_f_t visit_f_; // user-defined visit function
  long long length_;  // number of vertices visited so far

  // #### auxiliary functions ####

  // Execute the given flip sequence seq on the current vertex x, and do
  // the necessary termination checks along the way. Return value is true
  // if termination criterion is satisfied (break surrounding loop).
  bool flip_seq(const std::vector<int> &seq, int &dist_to_start,
                bool final_path);
};

#include <cassert>
#include <vector>

HamCycle::HamCycle(const Vertex &x, long long limit, visit_f_t visit_f)
    : x_(x), y_(x), limit_(limit), visit_f_(visit_f) {
  assert(this->x_.size() % 2 == 1);
  const int n = this->x_.size() / 2;

  // ############################
  // ### Initialization phase ###
  // ############################

  // The vertex xs is the first vertex of a path in H_n\circ 0 for which
  // either this path or the subsequent one in \overline{\rev}(H_n)\circ 1
  // contains the vertex this->x_.

  Vertex xs(this->x_);
  int skip = 0;         // number of vertices which are not visited (skipped)
                        // before reaching the first vertex x_
  if (xs[2 * n] == 1) { // last bit == 1
    // move backwards along the cycle to the first vertex of the
    // path in the graph \overline{\rev}(H_n)\circ 1 (this is the
    // last vertex of this oriented path, as the path is traversed
    // opposite to its orientation)
    xs.rev_inv();
    skip += xs.to_last_vertex();
    xs.rev_inv();
    xs[2 * n] = 0; // jump backwards to the graph H_n\circ 0 by
                   // flipping last bit
    skip++;
  }
  // move backwards along the cycle to the first vertex
  // of the path in the graph H_n\circ 0
  skip += xs.to_first_vertex();
  assert(xs.is_first_vertex());

  // maintain adjacency list representation of the tree
  // that corresponds to the current vertex this->y_ throughout the
  // algorithm
  this->y_ = xs;
  Tree y_tree(this->y_);

  // if initial cycle segment contains flipped paths,
  // we may need to start from the other path
  if ((skip > 0) && (y_tree.flip_tree())) {
    if ((xs[1] == 1) && (skip <= 5)) {
      skip = 6 - skip; // need to correct for the reverse traversal
                       // of the initial part of flipped path
    }
    int y_string[2 * n];
    y_tree.to_bitstring(y_string);
    std::vector<int> y_vec(y_string, y_string + 2 * n);
    y_vec.push_back(0); // add 0-bit at the end
    xs = Vertex(y_vec);
    this->y_ = xs;
  }

  this->length_ = 0;

  // ##################################
  // ### Hamilton cycle computation ###
  // ##################################
  std::vector<int> seq;
  std::vector<int> seq01;
  seq01.push_back(2 * n); // flip sequence that flips only the last bit
  int dist_to_start = skip;
  bool final_path = false; // back in the path that contains the starting vertex
  while (true) {
    // #################################################
    // follow the path in the graph H_n\circ 0
    // #################################################

    bool flip = y_tree.flip_tree(); // tau() or tau_inverse() is applied
                                    // inside the function call
    y_tree.rotate();

    // compute flip sequence
    this->y_.compute_flip_seq_0(seq, flip);

    // apply flip sequence
    assert(this->y_.is_first_vertex());
    if (flip_seq(seq, dist_to_start, final_path)) {
      break;
    }
    assert(this->y_.is_last_vertex());

    // flip last bit to jump to the graph \overline{\rev}(H_n)\circ
    // 1
    if (flip_seq(seq01, dist_to_start, final_path)) {
      break;
    }
    assert(this->y_[2 * n] == 1);

    // #################################################
    // follow the path in the graph \overline{\rev}(H_n)\circ 1
    // #################################################

    // ############# ORIGINAL SEQUENCE OF CODE AS DESCRIBED IN THE
    // PAPER
    // ############# compute flip sequence
    /*
    Vertex yp = y_;
    yp.rev_inv();
    yp.to_first_vertex();
    yp.compute_flip_seq_0(seq, false);
    // compute transformed flip sequence under the operations
    // of reversal and complementation (rev_inv)
    std::reverse(seq.begin(), seq.end());
    for (int j = 0; j < seq.size(); ++j) {
      seq[j] = 2*n-1 - seq[j];
    }
    */
    // ############# SPEED-OPTIMIZED VARIANT #############
    // compute transformed flip sequence directly
    this->y_.compute_flip_seq_1(seq);

    // apply flip sequence
    assert(this->y_.is_last_vertex());
    if (flip_seq(seq, dist_to_start, final_path)) {
      break;
    }
    assert(this->y_.is_first_vertex());

    // flip last bit to jump to the graph H_n\circ 0
    if (flip_seq(seq01, dist_to_start, final_path)) {
      break;
    }
    assert(this->y_[2 * n] == 0);

    // exit loop prematurely if starting vertex has been reached
    // again
    if (this->y_ == xs) {
      final_path = true;
      dist_to_start = skip;
    }
  }
}

bool HamCycle::flip_seq(const std::vector<int> &seq, int &dist_to_start,
                        bool final_path) {
  if ((dist_to_start > 0) || final_path ||
      ((this->limit_ >= 0) && (this->length_ + seq.size() >= this->limit_))) {
    // apply only part of the flip sequence
    for (int j = 0; j < seq.size(); ++j) {
      if ((final_path && (dist_to_start == 0)) ||
          ((this->limit_ >= 0) && (this->length_ == this->limit_))) {
        return true; // terminate Hamilton cycle
                     // computation prematurely
      }
      const int i = seq[j];
      if ((dist_to_start == 0) || final_path) {
        this->y_[i] = 1 - this->y_[i];
#ifndef NVISIT
        // the visit_f_() function is useful only when
        // nonempty
        visit_f_(this->y_.get_bits(), i);
#endif
        ++length_;
      } else {
        this->y_[i] = 1 - this->y_[i];
      }
      if (dist_to_start > 0) {
        dist_to_start--;
      }
    }
  } else {
    // highspeed loop without case distinctions
    // apply the entire flip sequence
    for (int j = 0; j < seq.size(); ++j) {
      const int i = seq[j];
      this->y_[i] = 1 - this->y_[i];
#ifndef NVISIT
      // the visit_f_() function is useful only when nonempty
      visit_f_(this->y_.get_bits(), i);
#endif
    }
    this->length_ += seq.size();
  }
  return false; // continue Hamilton cycle computation
}

#include <cstdlib>
#include <getopt.h>
#include <iostream>
#include <vector>

// display help
void help() {
  std::cout << "./middle [options]  compute middle levels Gray code from "
               "[Muetze,Nummenpalo]"
            << std::endl;
  std::cout << "-h                  display this help" << std::endl;
  std::cout << "-n{1,2,...}         list bitstrings of length 2n+1 with weight "
               "n or n+1"
            << std::endl;
  std::cout << "-l{-1,0,1,2,...}    number of bitstrings to list; -1 for "
               "full cycle"
            << std::endl;
  std::cout << "-v{0,1}^{2n+1}      initial bitstring (length 2n+1, "
               "weight n or n+1)"
            << std::endl;
  std::cout << "-s{0,1}             store and print all visited bitstrings "
               "(no=0, yes=1)"
            << std::endl;
  std::cout << "-p{0,1}             print the flip positions instead of "
               "bitstrings (no=0, yes=1)"
            << std::endl;
  std::cout << "examples:  ./middle -n2" << std::endl;
  std::cout << "           ./middle -n2 -v01010" << std::endl;
  std::cout << "           ./middle -n2 -p1" << std::endl;
  std::cout << "           ./middle -n10 -l50" << std::endl;
  std::cout << "           ./middle -n12 -s0" << std::endl;
}

void opt_n_missing() {
  std::cerr << "option -n is mandatory and must come before -v" << std::endl;
}

void opt_v_error() {
  std::cerr << "option -v must be followed by a bitstring of length 2n+1 with "
               "weight n or n+1"
            << std::endl;
}

// user-defined visit function
// y is the current vertex and i is the position of the last flip
void visit_f_empty(const std::vector<int> &y, int i) {
  // visit vertex y
}

std::vector<int> flip_seq_; // flip sequence

void visit_f_log(const std::vector<int> &y, int i) { flip_seq_.push_back(i); }

int main(int argc, char *argv[]) {
  int n;
  bool n_set = false;          // flag whether option -n is present
  long long limit = -1;        // compute all vertices by default
  std::vector<int> v;          // starting vertex
  bool v_set = false;          // flag whether option -v is present
  bool store_vertices = true;  // store vertices by default
  bool print_flip_pos = false; // print bitstrings by default

  // process command line options
  int c;
  while ((c = getopt(argc, argv, ":hn:l:v:s:p:")) != -1) {
    switch (c) {
    case 'h':
      help();
      return 0;
    case 'n':
      n = atoi(optarg);
      if (n < 1) {
        std::cerr << "option -n must be followed by an "
                     "integer from {1,2,...}"
                  << std::endl;
        return 1;
      }
      v.resize(2 * n + 1, 0);
      n_set = true;
      break;
    case 'l':
      limit = atoi(optarg);
      if (limit < -1) {
        std::cerr << "option -l must be followed by an "
                     "integer from {-1,0,1,2,...}"
                  << std::endl;
        return 1;
      }
      break;
    case 'v': {
      if (!n_set) {
        opt_n_missing();
        help();
        return 1;
      }
      char *p = optarg;
      int length = 0;
      int num_ones = 0;
      // parse bitstring
      while (*p != 0) {
        if ((*p == '0') || (*p == '1')) {
          v[length] = ((int)(*p)) - 48; // convert ASCII
                                        // character to bit
          length++;
          num_ones += ((*p) - 48);
          if (length > 2 * n + 1) {
            opt_v_error();
            return 1;
          }
        } else {
          opt_v_error();
          return 1;
        }
        p++;
      }
      // check length and weight of bitstring
      if ((length != 2 * n + 1) || (num_ones < n) || (num_ones > n + 1)) {
        opt_v_error();
        return 1;
      }
      v_set = true;
      break;
    }
    case 's': {
      const int arg = atoi(optarg);
      if ((arg < 0) || (arg > 1)) {
        std::cerr << "option -s must be followed by 0 or 1" << std::endl;
        return 1;
      }
      store_vertices = (bool)arg;
      break;
    }
    case 'p': {
      const int arg = atoi(optarg);
      if ((arg < 0) || (arg > 1)) {
        std::cerr << "option -p must be followed by 0 or 1" << std::endl;
        return 1;
      }
      print_flip_pos = (bool)arg;
      break;
    }
    case ':':
      std::cerr << "option -" << (char)optopt << " requires an operand"
                << std::endl;
      return 1;
    case '?':
      std::cerr << "unrecognized option -" << (char)optopt << std::endl;
      return 1;
    }
  }
  if (!n_set) {
    opt_n_missing();
    help();
    return 1;
  }

  // define a default starting vertex
  if (!v_set) {
    for (int i = 0; i < n; ++i) {
      v[i] = 1;
    }
  }
  Vertex x(v);

  // Starting from x, compute the next limit vertices on a Hamilton cycle.
  // A value limit < 0 means computing until we are back at the starting
  // vertex.
  visit_f_t visit_f = store_vertices ? visit_f_log : visit_f_empty;
  HamCycle hc(x, limit, visit_f);

  // print vertices encountered along the cycle
  if (store_vertices) {
    if (limit != 0) {
      std::cout << x << std::endl;
    }
    for (int j = 0; j < (int)flip_seq_.size() - 1; ++j) {
      const int i = flip_seq_[j];
      x[i] = 1 - x[i];      // flip bit
      if (print_flip_pos) { // print only flip positions
        std::cout << i << std::endl;
      } else { // print actual bitstring
        std::cout << x << std::endl;
      }
    }
    if (limit == flip_seq_.size()) {
      std::cout << "output limit reached" << std::endl;
    }
  }

  return 0;
}
