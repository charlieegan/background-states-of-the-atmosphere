#ifndef POOL_HPP
#define POOL_HPP

#include <vector>
#include <utility>

// a simple class for storing elements (in blocks of size BLOCKSIZE for efficiency) in a way that their addresses do not change and they are deleted when pool goes out of scope
template <typename T, int BLOCKSIZE=1024>
class pool
{
protected:
  std::vector<T*> data;
  int lcnt; // number of elements in last block

public:
  pool() : lcnt(BLOCKSIZE) {}
  pool(pool &&o) : data(o.data), lcnt(o.lcnt) {
    o.data.clear();
  }

  ~pool() {
    for (auto &p : data)
      delete[] p;
  }

  T &operator[](const int &i) {
    return data[i / BLOCKSIZE][i % BLOCKSIZE];
  }
  const T &operator[](const int &i) const {
    return data[i / BLOCKSIZE][i % BLOCKSIZE];
  }
  void push_back(const T &t) {
    if (lcnt == BLOCKSIZE) {
      data.push_back(new T[BLOCKSIZE]);
      lcnt = 0;
    }
    data.back()[lcnt++] = t;
  }
  void push_back(T &&t) {
    if (lcnt == BLOCKSIZE) {
      data.push_back(new T[BLOCKSIZE]);
      lcnt = 0;
    }
    data.back()[lcnt++] = std::move(t);
  }
  T &back() {
    return data.back()[lcnt - 1];
  }
  const T &back() const {
    return data.back()[lcnt - 1];
  }
  int size() const {
    return ((int)data.size() - 1) * BLOCKSIZE  + lcnt;
  }
};

#define BIND_POOL(m, T, name)                                           \
  py::class_<pool<T>>(m, name)                                          \
  .def(py::init<>())                                                    \
  .def("__getitem__", [](const pool<T> &r, const int &i){ return r[i]; }) \
  .def("__len__", pool<T>::size)                                        \
  .def("push_back", pool<T>::push_back)                                 \
  .def("back", pool<T>::back);

#endif
