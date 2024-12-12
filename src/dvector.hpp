#ifndef DVECTOR_HPP
#define DVECTOR_HPP

#include <vector>

/*! Dynamic vector template allowing efficient deletion of elements.
 * The size of the structure will be maximal intermediate size, stored elements always keep the same index.
 */
template <typename T>
class dvector
{
public:
  std::vector<T> data; //!< stored elements
  std::vector<bool> avail; //!< boolean set representation of available slots in data
  std::vector<int> avail_idx; //!< index representation of available slots in data

  /*! Construct empty dvector. */
  dvector() : data(), avail(), avail_idx() {}

  /*! Construct dvector containing n elements (default constructed as in std::vector). */
  dvector(const int &n) : data(n), avail(n, false), avail_idx() {}
  
  /*! Add an element and return it's index (in amortized O(1)).
   * \param t value to add
   * \return index such that (*this)[i] is the element that was added
   */
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

  /*! Remove element at index (in amortized O(1)).
   * \param i index to remove element from
   */
  void remove(const int &i) {
    avail[i] = true;
    avail_idx.push_back(i);
   }

  /*! Check whether index i is currently in use (in O(1)).
   * \param i index to check
   * \return whether *this currently stores a value at i
   */
  bool contains_idx(const int &i) const {
    return i >= 0 && i < (int)data.size() && !avail[i];
  }

  /*! Return number of contained elements (in O(1))
   * \return number of elements currently stored in *this
   */
  int size() const {
    return (int)data.size() - (int)avail_idx.size();
  }

  /*! Return number of elements that could be stored without needing to allocate new memory.
   * \return current size of data buffer
   */
  int capacity() const {
    return (int)data.size();
  }

  /*! Read/Write access element at given index.
   * \param i index to read / write
   * \return reference to element at index i
   */
  T& operator[](const int &i) {
    return data[i];
  }

  /*! Read access element at given index.
   * \param i index to read / write
   * \return const reference to element at index i
   */
  const T& operator[](const int &i) const {
    return data[i];
  }

  /*! Reserve n elements in underlying vectors (to reduce move operations performance). This does not actually resize the vectors, so the change will not be seen in dvector::capacity().
   * \param n number of elements to reserve
   */
  void reserve(const int &n) {
    data.reserve(n);
    avail.reserve(n);
  }

  /*! A forward iterator class over stored elements in a dvector. */
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
      while (++i < (int)v->avail.size() && v->avail[i]);
      return *this;
    }
    bool operator==(const iterator &it) const {
      return i == it.i;
    }
    bool operator<(const iterator &it) const {
      return i < it.i;
    }
  };

  /*! A const forward iterator class over stored elements in a dvector. */
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

  /*! A const forward iterator class over indices of stored elements in a dvector. */
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

    int operator*() {
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

  /*! Get forward iterator pointing to start of dvector stored elements. */
  iterator begin() {
    return iterator(this, 0);
  }
  /*! Get forward iterator pointing to end of dvector stored elements. */
  iterator end() {
    return iterator(this, avail.size());
  }

  /*! Get const forward iterator pointing to start of dvector stored elements. */
  const_iterator begin() const {
    return const_iterator(this, 0);
  }
  /*! Get const forward iterator pointing to end of dvector stored elements. */
  const_iterator end() const {
    return const_iterator(this, avail.size());
  }

  /*! Get const forward iterator pointing to start of dvector stored element indices. */
  const_idx_iterator begin_idx() const {
    return const_idx_iterator(this, 0);
  }
  /*! Get const forward iterator pointing to end of dvector stored element indices. */
  const_idx_iterator end_idx() const {
    return const_idx_iterator(this, avail.size());
  }

  /*! Output stored elements to std::ostream, assuming the elements themselves have operator<< overloaded. */
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

  /*! Return stored elements compressed into a continuous std::vector. */
  std::vector<T> to_vector() const {
    std::vector<T> res;
    for (const auto &x : *this)
      res.push_back(x);
    return res;
  }

  /*! Transform index of this to corresponding index in compressed version as returnd by to_vector.
   * \param i index into *this
   * \return index j such that this->to_vector()[j] == (*this)[i]
   */
  int transform_idx(const int &i) const {
    int sh = 0;
    for (auto k : avail_idx) {
      sh += (k <= i);
#ifdef DEBUG_CHECKS
      if (k == i)
        throw std::runtime_error("tried to transform_idx invalid index");
#endif
    }
    return i - sh;
  }

  /*! Represent contents in a string. */
  std::string repr() const {
    std::stringstream ss;
    ss << this;
    return ss.str();
  }
};

#endif
