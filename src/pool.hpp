#ifndef POOL_HPP
#define POOL_HPP

#include <vector>
#include <utility>
#include <memory>

/*!
 * A simple class for efficiently storing objects with dynamic allocation.
 * Memory lifetime is bound to the pool. It is freed when the pool is destructed.
 * Objects can only be created, not deleted. Objects can be accessed by index but are not stored
 * contiguously in memory (but rather in blocks of size BLOCKSIZE).
 */
template <typename T, int BLOCKSIZE=1024, class Allocator=std::allocator<T>>
class pool
{
protected:
  Allocator alloc; //!< allocator to use
  std::vector<T*> data; //!< vector of memory block pointers
  int lcnt; //!< number of elements currently stored in last block

public:
  /*! Construct empty pool with optional custom allocator. */
  pool(const Allocator &alloc = Allocator()) : alloc(alloc), lcnt(BLOCKSIZE) {}

  /*! Move another pool into this one. */
  pool(pool &&o) : alloc(o.alloc), data(o.data), lcnt(o.lcnt) {
    o.data.clear();
  }

  /*! Destruct pool, freeing all memory. */
  ~pool() {
    for (auto &p : data)
      alloc.deallocate(p, BLOCKSIZE);
  }

  /*! (Read/Write) access element at given index.
   * \param i index to access
   * \return reference to element at index i
   */
  T &operator[](const int &i) {
    return data[i / BLOCKSIZE][i % BLOCKSIZE];
  }
  /*! (Read only) access element at given index.
   * \param i index to access
   * \return const reference to element at index i
   */
  const T &operator[](const int &i) const {
    return data[i / BLOCKSIZE][i % BLOCKSIZE];
  }

  /*! Add new element to pool (in amortized O(1), worst case O(BLOCKSIZE)) .
   * \param t value to add
   */
  void push_back(const T &t) {
    if (lcnt == BLOCKSIZE) {
      data.push_back(alloc.allocate(BLOCKSIZE)); //new T[BLOCKSIZE]);
      lcnt = 0;
    }
    data.back()[lcnt++] = t;
  }

  /*! Move-add new element to pool (in amortized O(1), worst case O(BLOCKSIZE)) .
   * \param t value to add
   */
  void push_back(T &&t) {
    if (lcnt == BLOCKSIZE) {
      data.push_back(alloc.allocate(BLOCKSIZE)); //new T[BLOCKSIZE]);
      lcnt = 0;
    }
    data.back()[lcnt++] = std::move(t);
  }

  /*! (Read/Write) access last element added.
   * \return reference to last element added.
   */
  T &back() {
    return data.back()[lcnt - 1];
  }
  
  /*! (Read only) access last element added.
   * \return const reference to last element added.
   */
  const T &back() const {
    return data.back()[lcnt - 1];
  }

  /*! Get number of stored elements.
   * \return number of elements currently stored
   */
  int size() const {
    return ((int)data.size() - 1) * BLOCKSIZE  + lcnt;
  }

  /*! Add (partial) python bindings to m */
  static void bind(py::module_ &m) {
    py::class_<pool<T>>(m, ("Pool_" + type_name<T>::value()).c_str())
      .def(py::init<>())
      .def("__getitem__", [](const pool<T> &r, const int &i){ return r[i]; })
      .def("__len__", pool<T>::size)
      .def("push_back", pool<T>::push_back)
      .def("back", pool<T>::back);

  }
};

#endif
