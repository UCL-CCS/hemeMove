// -*- mode: c++; -*-
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#ifndef HLBGMYTOOL_OCT_OCT_HPP
#define HLBGMYTOOL_OCT_OCT_HPP

#include <ostream>
#include <stdexcept>
#include <string>
#include "Oct.h"

namespace hemelb::gmytool::oct {

template <class T>
Octree<T>::Node::Node(Int i, Int j, Int k, Int l)
    : x(i), y(j), z(k), level(l), value() {}
template <class T>
T& Octree<T>::Node::Data() {
  return value;
}
template <class T>
const T& Octree<T>::Node::Data() const {
  return value;
}

template <class T>
auto Octree<T>::Node::X() const -> Int {
  return x;
}
template <class T>
auto Octree<T>::Node::Y() const -> Int {
  return y;
}
template <class T>
auto Octree<T>::Node::Z() const -> Int {
  return z;
}
template <class T>
auto Octree<T>::Node::Level() const -> Int {
  return level;
}

template <class T>
bool Octree<T>::Node::IsNodeInRange(Int i, Int j, Int k, Int l) const {
  // the least significant this->level bits represent the relative index
  // the other bits must match my bits
  auto tl = this->level;
  return ((i >> tl) == (this->x >> tl)) && ((j >> tl) == (this->y >> tl)) &&
         ((k >> tl) == (this->z >> tl));
}

// Get a node without creating - returns null pointer if doensn't exist
template <class T>
auto Octree<T>::Branch::Get(Int i, Int j, Int k, Int l) -> NodePtr {
  if (l >= this->level) {
    // a parent - error
    throw std::out_of_range("trying to get a parent node");
  }

  if (!this->IsNodeInRange(i, j, k, l))
    throw std::out_of_range("requested node not in my range");

  // OK
  return get_nocreate_internal(i, j, k, l);
}

// Get a node without creating - returns null pointer if doensn't exist
template <class T>
auto Octree<T>::Branch::Get(Int i, Int j, Int k, Int l) const -> ConstNodePtr {
  if (l >= this->level) {
    // a parent - error
    throw std::out_of_range("trying to get a parent node");
  }

  if (!this->IsNodeInRange(i, j, k))
    throw std::out_of_range("requested node not in my range");

  // OK
  return get_nocreate_internal(i, j, k, l);
}

// Get a node without creating - returns null pointer if doensn't exist
template <class T>
auto Octree<T>::Branch::GetCreate(Int i, Int j, Int k, Int l) -> NodePtr {
  if (l > this->level) {
    // a parent - error
    throw std::out_of_range("trying to get a parent node");
  }

  if (!this->IsNodeInRange(i, j, k))
    throw std::out_of_range("requested node not in my range");

  // OK
  if (l == this->level) {
    return this->shared_from_this();
  }

  return get_create_internal(i, j, k, l);
}

template <class T>
auto Octree<T>::Branch::GetCreatePath(Int i, Int j, Int k, Int l) -> NodeList {
  // if (l >= this->level) {
  //   // a parent - error
  //   throw std::out_of_range("trying to get a parent node");
  // }

  // if (!this->IsNodeInRange(i, j, k))
  //   throw std::out_of_range("requested node not in my range");

  if (l == this->level) {
    return NodeList(1, this->shared_from_this());
  }

  auto child = GetCreate(i, j, k, this->level - 1);

  auto ans = child->GetCreatePath(i, j, k, l);
  ans.push_front(this->shared_from_this());
  return ans;
}

template <class T>
void Octree<T>::Branch::Set(Int i, Int j, Int k, Int l, NodePtr n) {
  Int pl = l + 1;
  Int mask = ~(1 << l);
  Int pi = i & mask;
  Int pj = j & mask;
  Int pk = k & mask;

  auto parent = std::dynamic_pointer_cast<Branch>(GetCreate(pi, pj, pk, pl));
  if (!parent)
    throw std::out_of_range("parent of request index not a branch node");

  parent->children[Octree<T>::Node::LocalIndex(i, j, k, l)] = n;
}

template <class T>
void Octree<T>::Branch::Accept(Visitor& v) {
  v.Arrive(this->shared_from_this());
  if (v.ShouldDescend(this->shared_from_this()))
    for (auto i = 0; i < 8; ++i) {
      auto child = children[i];
      if (child)
        child->Accept(v);
    }

  v.Depart(this->shared_from_this());
}
template <class T>
void Octree<T>::Branch::Accept(ConstVisitor& v) const {
  v.Arrive(this->shared_from_this());
  if (v.ShouldDescend(this->shared_from_this()))
    for (auto i = 0; i < 8; ++i) {
      auto child = children[i];
      if (child)
        child->Accept(v);
    }

  v.Depart(this->shared_from_this());
}

template <class T>
auto Octree<T>::Branch::get_create_internal(Int i, Int j, Int k, Int l)
    -> NodePtr {
  Int child_level = this->level - 1;
  // Get the local index
  auto i_child = Octree<T>::Node::LocalIndex(i, j, k, child_level);

  Int li = (i >> child_level) & 1;
  Int lj = (j >> child_level) & 1;
  Int lk = (k >> child_level) & 1;
  auto child = children[i_child];

  if (!child) {
    // Create the child node
    Int ci = this->x + (li << child_level);
    Int cj = this->y + (lj << child_level);
    Int ck = this->z + (lk << child_level);

    if (child_level) {
      // child is branch
      child = children[i_child] = NodePtr(new Branch(ci, cj, ck, child_level));
    } else {
      // child is leaf
      child = children[i_child] = NodePtr(new Leaf(ci, cj, ck, child_level));
    }
  }

  if (child_level == l) {
    return child;
  } else {
    auto branch = dynamic_cast<Branch*>(child.get());
    return branch->get_create_internal(i, j, k, l);
  }
}

template <class T>
auto Octree<T>::Branch::get_nocreate_internal(Int i, Int j, Int k, Int l)
    -> NodePtr {
  // Get the local index
  Int child_level = this->level - 1;
  auto child = children[Octree<T>::Node::LocalIndex(i, j, k, child_level)];

  if (!child) {
    return child;
  }

  if (child_level == l) {
    return child;
  } else {
    auto branch = dynamic_cast<Branch*>(child.get());
    return branch->get_nocreate_internal(i, j, k, l);
  }
}
template <class T>
auto Octree<T>::Branch::get_nocreate_internal(Int i, Int j, Int k, Int l) const
    -> ConstNodePtr {
  // Get the local index
  Int child_level = this->level - 1;
  auto child = children[Octree<T>::Node::LocalIndex(i, j, k, child_level)];

  if (!child) {
    return child;
  }

  if (child_level == l) {
    return child;
  } else {
    auto branch = dynamic_cast<Branch*>(child.get());
    return branch->get_nocreate_internal(i, j, k, l);
  }
}

// Get a node without creating
template <class T>
auto Octree<T>::Leaf::Get(Int i, Int j, Int k, Int l) -> NodePtr {
  throw std::out_of_range("trying to get a child of a leaf node");
}
template <class T>
auto Octree<T>::Leaf::Get(Int i, Int j, Int k, Int l) const -> ConstNodePtr {
  throw std::out_of_range("trying to get a child of a leaf node");
}
template <class T>
auto Octree<T>::Leaf::GetCreate(Int i, Int j, Int k, Int l) -> NodePtr {
  throw std::out_of_range("trying to get a child of a leaf node");
}
template <class T>
auto Octree<T>::Leaf::GetCreatePath(Int i, Int j, Int k, Int l) -> NodeList {
  return NodeList(1, this->shared_from_this());
}

template <class T>
void Octree<T>::Leaf::Set(Int i, Int j, Int k, Int l, NodePtr n) {
  throw std::out_of_range("trying to set a child of a leaf node");
}

template <class T>
void Octree<T>::Leaf::Accept(Visitor& v) {
  v.Arrive(this->shared_from_this());
  v.Depart(this->shared_from_this());
}

template <class T>
void Octree<T>::Leaf::Accept(ConstVisitor& v) const {
  v.Arrive(this->shared_from_this());
  v.Depart(this->shared_from_this());
}

template <class T>
Octree<T>::Octree(Int level_) : level(level_) {
  if (level) {
    root = NodePtr(new Branch(0, 0, 0, level));
  } else {
    root = NodePtr(new Leaf(0, 0, 0, 0));
  }
}
template <class T>
auto Octree<T>::Root() -> NodePtr {
  return root;
}
template <class T>
auto Octree<T>::Root() const -> ConstNodePtr {
  return root;
}

template <class T>
auto Octree<T>::Get(Int i, Int j, Int k, Int l) -> NodePtr {
  return root->Get(i, j, k, l);
}
template <class T>
auto Octree<T>::Get(Int i, Int j, Int k, Int l) const -> ConstNodePtr {
  return root->Get(i, j, k, l);
}

template <class T>
auto Octree<T>::GetCreate(Int i, Int j, Int k, Int l) -> NodePtr {
  return root->GetCreate(i, j, k, l);
}
template <class T>
void Octree<T>::Set(Int i, Int j, Int k, Int l, NodePtr n) {
  root->Set(i, j, k, l, n);
}

template <class T>
void WriteNodeData(std::ostream& os, const T& data);

template <class T>
class Writer : public Octree<T>::ConstVisitor {
 public:
  typedef Octree<T> Tree;
  typedef typename Tree::ConstNodePtr ConstNodePtr;

  Writer(std::ostream& stream) : os(stream), padding(0) {}

  virtual void Arrive(ConstNodePtr node) {
    os << std::string(padding, ' ') << node->X() << ", " << node->Y() << ", "
       << node->Z() << " = ";
    WriteNodeData(os, node->Data());
    os << std::endl;
    padding++;
  }
  virtual void Depart(ConstNodePtr node) { padding--; }

 private:
  std::ostream& os;
  int padding;
};

template <class T>
std::ostream& operator<<(std::ostream& os, Octree<T>& tree) {
  os << "Octree levels = " << tree.Level() << std::endl;
  Writer<T> w(os);
  tree.Root()->Accept(w);
  return os;
}
// template<class T>
// std::ostream& operator<<(std::ostream& os, typename Octree<T>::Node& node) {
//   Writer<T> w(os);
//   node.Accept(w);
//   return os;
// }

template <class T>
void WriteNodeData(std::ostream& os, const T& data) {
  os << data;
}

}  // namespace hemelb::gmytool::oct
#endif
