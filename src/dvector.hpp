#ifndef DVECTOR_HPP
#define DVECTOR_HPP

#include <vector>

// dynamic vector template allowing efficient deletion of elements (size of the structure will be maximal intermediate size)
template <typename T>
class dvector
{
protected:
  std::vector<T> data;
  std::vector<bool> avail;
  std::vector<int> avail_idx;
  
public:

  dvector() : data(), avail(), avail_idx() {}
  dvector(const int &n) : data(n), avail(n, false), avail_idx() {}
  
  
  // add an element, return it's index
  int add(const T &t) {
    if (avail_idx.size()) {
      int i = avail_idx.back();
      avail_idx.pop_back();
      data[i] = t;
      avail[i] = false;
      return i;
    } else {
      int i = data.size();
      data.push_back(t);
      avail.push_back(false);
      return i;
    }
  }

  // remove element at index
  void remove(const int &i) {
    avail[i] = true;
    avail_idx.push_back(i);
  }

  // check whether index i is contained
  bool contains_idx(const int &i) const {
    return i >= 0 && i < (int)data.size() && !avail[i];
  }
  
  int size() const {
    return (int)data.size() - (int)avail_idx.size();
  }

  int capacity() const {
    return (int)data.size();
  }

  T& operator[](const int &i) {
    return data[i];
  }
  const T& operator[](const int &i) const {
    return data[i];
  }

  void reserve(const int &n) {
    data.reserve(n);
    avail.reserve(n);
  }
  
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

    T& operator*() {
      return v->data[i];
    }
    iterator& operator++() {
      while (++i < v->avail.size() && v->avail[i]);
      return *this;
    }
    bool operator==(const iterator &it) const {
      return i == it.i;
    }
    bool operator<(const iterator &it) const {
      return i < it.i;
    }
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

    const T& operator*() {
      return v->data[i];
    }
    const_iterator& operator++() {
      while (++i < (int)v->avail.size() && v->avail[i]);
      return *this;
    }
    bool operator==(const const_iterator &it) const {
      return i == it.i;
    }
    bool operator<(const const_iterator &it) const {
      return i < it.i;
    }
  };

  class const_idx_iterator
  {
    friend dvector;
  protected:
    const dvector *v;
    int i;

    const_idx_iterator(const dvector *v, const int &j) : v(v), i(j) {
      while (i < (int)v->avail.size() && v->avail[i]) i++;
    }
    
  public:
    const_idx_iterator() : v(NULL), i(0) { }
    const_idx_iterator(const_idx_iterator &&it) : v(it.v), i(it.i) {}
    const_idx_iterator(const const_idx_iterator &it) : v(it.v), i(it.i) {}

    const int& operator*() {
      return i;
    }
    const_idx_iterator& operator++() {
      while (++i < (int)v->avail.size() && v->avail[i]);
      return *this;
    }
    bool operator==(const const_idx_iterator &it) const {
      return i == it.i;
    }
    bool operator<(const const_idx_iterator &it) const {
      return i < it.i;
    }
  };

  iterator begin() {
    return iterator(this, 0);
  }
  iterator end() {
    return iterator(this, avail.size());
  }

  const_iterator begin() const {
    return const_iterator(this, 0);
  }
  const_iterator end() const {
    return const_iterator(this, avail.size());
  }

  const_idx_iterator begin_idx() const {
    return const_idx_iterator(this, 0);
  }
  const_idx_iterator end_idx() const {
    return const_idx_iterator(this, avail.size());
  }

  
  friend std::ostream& operator<<(std::ostream &os, const dvector &v) {
    os << "[";
    bool sep = false;
    for (const auto &t : v) {
      if (sep)
        os << ", ";
      os << t;
      sep = true;
    }
    os << "]";
    return os;
  }

  std::vector<T> to_vector() const {
    std::vector<T> res;
    for (const auto &x : *this)
      res.push_back(x);
    return res;
  }
  
  std::string repr() const {
    std::stringstream ss;
    ss << this;
    return ss.str();
  }
};

#endif
