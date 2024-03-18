#ifndef MIN_HPP_INCLUDE
#define MIN_HPP_INCLUDE

#include <vector>

/**
 * @brief A function template that returns the 
 * minimum between two objects. Note: Class K
 * must have defined comparator.
 * @tparam K data type
 * @param a first object
 * @param b second object
 * @return K minimum between a and b
 */
template <class K>
K min(const K& a, const K& b) {
    if (a < b) return a;
    return b;
}

/**
 * @brief A function template that returns the 
 * minimum among the elements in the vector x.
 * Note: Class K must have defined comparator.
 * @tparam K data type 
 * @param x vector containing elements of type K
 * @return K minimum element in x
 */
template <class K>
K min(const std::vector<K>& x) {
    K ans = x[0];
    for (size_t i = 1; i < x.size(); i++) {
        if (x[i] < ans) ans = x[i];
    }
    return ans;
}

/**
 * @brief A function template that returns the 
 * minimum among the elements in the vector x
 * and stores the index of minimum element in the idx.
 * Note: Class K must have defined comparator.
 * @tparam K data type
 * @param x vector containing elements of type K
 * @param idx index of the minimum element
 * @return K minimum element
 */
template <class K>
K min(const std::vector<K> &x, long& idx) {
    K v = x[0];
    idx = 0;
    long n = long(x.size());
    for (long j = 1; j < n; j++) {
        if (x[j] < v) {
            v = x[j];
            idx = j;
        }
    }
    return v;
}

#endif