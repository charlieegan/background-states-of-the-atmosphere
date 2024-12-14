#ifndef REUSABLE_USED_ARRAY_HPP
#define REUSABLE_USED_ARRAY_HPP

#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <concepts>

/*! simple reusable used array that can mark, check, reset in amortized O(1)
 * a smaller type I will lead to smaller memory usage but more required hard resets (reset() will be O(n) every maxval(I) executions)
 * note that memory from this is only freed when the object is deleted
 */
template <std::signed_integral I = int>
class reusable_used_array
{
public:
  std::vector<I> a; //!< vector to store last-marked-times in
  I ctr; // current time (if a[i] == ctr, i is marked used)

  /*! Construct with initial size (default 0). Note that size does not influence behavior, only performance.
   * \param init_size how much space to initially reserve
   */
  reusable_used_array(const int &init_size = 0) : a(init_size, 0), ctr(1) {}

  /*! Reset array (i.e. mark everything unused). */
  void reset() {
    ctr++;
    if (ctr < 0) {
      std::fill(a.begin(), a.end(), 0);
      ctr = 1;
    }
  }

  /*! Check whether given elementis used.
   * \param i element to check for
   */
  bool check(const int &i) const {
    if (i < (int)a.size())
      return a[i] == ctr;
    return false;
  }

  /*! Mark given element as used.
   * \param i element to mark
   */
  void mark(const int &i) {
    while ((int)a.size() < i + 1)
      a.push_back(0);
    a[i] = ctr;
  }

  /*! Output as binary representation. */
  friend std::ostream& operator<<(std::ostream &os, const reusable_used_array &a) {
    os << "used[";
    for (int i = 0; i < (int)a.a.size(); ++i)
      os << (int)a.check(i);
    os << "]";
    return os;
  }

  /*! Generate and return binary string representation. */
  std::string repr() const {
    std::stringstream ss;
    ss << this;
    return ss.str();
  }
};

#endif
