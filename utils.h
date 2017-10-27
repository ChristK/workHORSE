#ifndef KMISC_UTILS_H_
#define KMISC_UTILS_H_

inline bool IsNA( const Rcpp::internal::const_string_proxy<STRSXP>& x ) {
  return x == NA_STRING;
}

inline bool IsNA(int x) {
  return x == NA_INTEGER;
}

inline bool IsNA(SEXP x) {
  return x == NA_STRING;
}

inline bool IsNA(double x) {
  return Rcpp::internal::Rcpp_IsNA(x);
}

inline bool IsNaN(double x) {
  return Rcpp::internal::Rcpp_IsNaN(x);
}

// borrowed from data.table package;
// see https://github.com/arunsrinivasan/datatable/blob/master/pkg/src/countingcharacter.c
inline int StrCmp(SEXP x, SEXP y)
{
    if (x == NA_STRING) return (y == NA_STRING ? 0 : 1);
    else if (y == NA_STRING) return -1;
    else if (x == y) return 0;  // same string in cache
    else return strcmp(char_nocheck(x), char_nocheck(y));
}

template <typename T>
struct NACompare;

template <>
struct NACompare<int> {
  inline bool operator()(int left, int right) const {
    if (left == NA_INTEGER) return false;
    if (right == NA_INTEGER) return true;
    return left < right;
  }
};

template <>
struct NACompare<double> {
  inline bool operator()(double left, double right) const {

    bool leftNaN = (left != left);
    bool rightNaN = (right != right);

    // this branch inspired by data.table: see
    // https://github.com/arunsrinivasan/datatable/commit/1a3e476d3f746e18261662f484d2afa84ac7a146#commitcomment-4885242
    if (IsNaN(right) and IsNA(left)) return true;

    if (leftNaN != rightNaN) {
      return leftNaN < rightNaN;
    } else {
      return left < right;
    }

  }

};

template <>
struct NACompare<SEXP> {
  inline bool operator()(SEXP left, SEXP right) const {
    return StrCmp(left, right) < 0;
  }
};

#endif
