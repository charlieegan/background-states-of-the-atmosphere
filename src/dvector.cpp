#include "dvector.hpp"

template <typename T>
int dvector<T>::add(const T &t) {
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

template <typename T>
void dvector<T>::remove(const int &i) {
  avail[i] = true;
  avail_idx.push_back(i);
}

template <typename T>
bool dvector<T>::contains_idx(const int &i) const {
  return i >= 0 && i < (int)data.size() && !avail[i];
}

template <typename T>
int dvector<T>::size() const {
  return (int)data.size() - (int)avail_idx.size();
}

template <typename T>
int dvector<T>::capacity() const {
  return (int)data.size();
}

template <typename T>
T& dvector<T>::operator[](const int &i) {
  return data[i];
}

template <typename T>
const T& dvector<T>::operator[](const int &i) const {
  return data[i];
}

template <typename T>
void dvector<T>::reserve(const int &n) {
  data.reserve(n);
  avail.reserve(n);
}

template <typename T>
dvector<T>::iterator dvector<T>::begin() {
  return dvector<T>::iterator(this, 0);
}

template <typename T>
dvector<T>::iterator dvector<T>::end() {
  return dvector<T>::iterator(this, avail.size());
}

template <typename T>
dvector<T>::const_iterator dvector<T>::begin() const {
  return dvector<T>::const_iterator(this, 0);
}

template <typename T>
dvector<T>::const_iterator dvector<T>::end() const {
  return dvector<T>::const_iterator(this, avail.size());
}

template <typename T>
dvector<T>::const_idx_iterator dvector<T>::begin_idx() const {
  return dvector<T>::const_idx_iterator(this, 0);
}

template <typename T>
dvector<T>::const_idx_iterator dvector<T>::end_idx() const {
  return dvector<T>::const_idx_iterator(this, avail.size());
}

template <typename T>
std::vector<T> dvector<T>::to_vector() const {
  std::vector<T> res;
  for (const auto &x : *this)
    res.push_back(x);
  return res;
}

template <typename T>
int dvector<T>::transform_idx(const int &i) const {
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

template <typename T>
std::string dvector<T>::repr() const {
  std::stringstream ss;
  ss << this;
  return ss.str();
}

template <typename T>
std::ostream& operator<<(std::ostream &os, const dvector<T> &v) {
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


template <typename T>
T& dvector<T>::iterator::operator*() {
  return v->data[i];
}

template <typename T>
dvector<T>::iterator& dvector<T>::iterator::operator++() {
  while (++i < (int)v->avail.size() && v->avail[i]);
  return *this;
}

template <typename T>
bool dvector<T>::iterator::operator==(const dvector<T>::iterator &it) const {
  return i == it.i;
}

template <typename T>
bool dvector<T>::iterator::operator<(const dvector<T>::iterator &it) const {
  return i < it.i;
}


template <typename T>
const T& dvector<T>::const_iterator::operator*() {
  return v->data[i];
}

template <typename T>
dvector<T>::const_iterator& dvector<T>::const_iterator::operator++() {
  while (++i < (int)v->avail.size() && v->avail[i]);
  return *this;
}

template <typename T>
bool dvector<T>::const_iterator::operator==(const dvector<T>::const_iterator &it) const {
  return i == it.i;
}

template <typename T>
bool dvector<T>::const_iterator::operator<(const dvector<T>::const_iterator &it) const {
  return i < it.i;
}


template <typename T>
dvector<T>::const_idx_iterator::const_idx_iterator(const dvector<T> *v, const int &j) : v(v), i(j) {
  while (i < (int)v->avail.size() && v->avail[i]) i++;
}


template <typename T>
int dvector<T>::const_idx_iterator::operator*() {
  return i;
}

template <typename T>
dvector<T>::const_idx_iterator& dvector<T>::const_idx_iterator::operator++() {
  while (++i < (int)v->avail.size() && v->avail[i]);
  return *this;
}

template <typename T>
bool dvector<T>::const_idx_iterator::operator==(const dvector<T>::const_idx_iterator &it) const {
  return i == it.i;
}

template <typename T>
bool dvector<T>::const_idx_iterator::operator<(const dvector<T>::const_idx_iterator &it) const {
  return i < it.i;
}
