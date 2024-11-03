
#include "vertex.hpp"
#include <algorithm>
#include <cassert>
#include <cstring>
#include <list>
#include <vector>

// A class to represent and manipulate an ordered rooted tree in doubly linked
// adjacency list representation.

class Tree {
public:
  explicit Tree(const Vertex &x);
  bool flip_tree();
  void rotate();
  void to_bitstring(int x[]) const;

private:
  int num_vertices_;
  int root_;
  std::vector<std::list<int>> children_;
  std::vector<int> parent_;

  // #### auxiliary functions ####

  // the degree of a vertex
  int deg(int u) const;
  // the number of children of vertex u
  int num_children(int u) const;

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
  void root_canonically();
  void compute_center(int &c1, int &c2) const;
  bool is_flip_tree_tau(); // We promise that we do not modify the tree in

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

  int min_string_rotation(int x[], int length);
};

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
