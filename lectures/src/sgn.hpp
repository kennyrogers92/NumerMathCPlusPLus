#ifndef SGN_H
#define SGN_H

template <class K>
int sgn(const K& x) {
     if (x < K(0)) return -1;
     if (x > K(0)) return 1;
     return 0;
}

template <class K>
K abs(const K& x) {
     return sgn(x)*x;
}

#endif
