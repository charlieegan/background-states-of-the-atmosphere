#ifndef DVECTOR_HPP
#define DVECTOR_HPP

#include <vector>
#include <iostream>
#include <sstream>
#include <pybind11/pybind11.h>

// dynamic vector template allowing efficient deletion of elements (size of the structure will be maximal intermediate size)
template <typename T>
class PYBIND11_EXPORT dvector
{
public:
  std::vector<T> data;
  std::vector<bool> avail;
  std::vector<int> avail_idx;

  dvector() : data(), avail(), avail_idx() {}
  dvector(const int &n) : data(n), avail(n, false), avail_idx() {}
  
  
  // add an element, return it's index
  int add(const T &t);
  // remove element at index
  void remove(const int &i);
  // check whether index i is contained
  bool contains_idx(const int &i) const;
  // return number of elements currently stored
  int size() const;
  // return number of elements that can be stored without resizing (which is an upper bound for the largest index)
  int capacity() const;
  /// access element at index
  T& operator[](const int &i);
  const T& operator[](const int &i) const;
  // reserve at least n elements for more efficient incremental adding
  void reserve(const int &n);
  
  class iterator
  {
    friend dvector;
  protected:
    dvector *v;
    int i;

    iterator(dvector *v, const int &i) : v(v), i(i) {}
    
  public:
    iterator() : v(NULL), i(0) {}
    iterator(iterator &&it) : v(it.v), i(it.i) {}
    iterator(const iterator &it) : v(it.v), i(it.i) {}

    T& operator*();
    iterator& operator++();
    bool operator==(const iterator &it) const;
    bool operator<(const iterator &it) const;
  };

  class const_iterator
  {
    friend dvector;
  protected:
    const dvector *v;
    int i;

    const_iterator(const dvector *v, const int &i) : v(v), i(i) {}
    
  public:
    const_iterator() : v(NULL), i(0) {}
    const_iterator(const_iterator &&it) : v(it.v), i(it.i) {}
    const_iterator(const const_iterator &it) : v(it.v), i(it.i) {}

    const T& operator*();
    const_iterator& operator++();
    bool operator==(const const_iterator &it) const;
    bool operator<(const const_iterator &it) const;
  };

  class const_idx_iterator
  {
    friend dvector;
  protected:
    const dvector *v;
    int i;

    const_idx_iterator(const dvector *v, const int &j);
    
  public:
    const_idx_iterator() : v(NULL), i(0) { }
    const_idx_iterator(const_idx_iterator &&it) : v(it.v), i(it.i) {}
    const_idx_iterator(const const_idx_iterator &it) : v(it.v), i(it.i) {}

    int operator*();
    const_idx_iterator& operator++();
    bool operator==(const const_idx_iterator &it) const;
    bool operator<(const const_idx_iterator &it) const;
  };

  iterator begin();
  iterator end();

  const_iterator begin() const;
  const_iterator end() const;

  const_idx_iterator begin_idx() const;
  const_idx_iterator end_idx() const;

  
  friend std::ostream& operator<<(std::ostream &os, const dvector &v);

  std::vector<T> to_vector() const;

  // transform index of this to corresponding index in compressed version as returnd by to_vector
  int transform_idx(const int &i) const;
  
  std::string repr() const;
};

#endif
