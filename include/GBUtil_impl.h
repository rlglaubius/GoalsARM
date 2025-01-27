#ifndef GBUTIL_IMPL_H
#define GBUTIL_IMPL_H

namespace GB {

// Split a delimited string into a vector of strings. For example,
// split("abc,def", ',') returns a vector containing "abc" and "def".
std::vector<std::string> split(const std::string &str, const char delimiter) {
  std::vector<std::string> rval;
  std::istringstream ss(str);
  std::string token;
  while (std::getline(ss, token, delimiter))
    rval.push_back(token);
  return rval;
}

// Return TRUE if the string str contains the substring sub, FALSE otherwise
bool contains(const std::string &str, const std::string &sub) {
  return str.find(sub) != std::string::npos;
}

// This is pretty gross but quality of floats raises an error as check will often
// fail due to machine precision. Better to check with some epsilon. But we have
// a special case here where something is initialized to 0.0, and
// we want to check for exact equality to that, so an == 0.0 is appropriate
// turn off the warnings for this single line.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
#pragma warning(push)
#pragma warning(disable : 4723)
template <typename T>
bool isZero(T x) {
    static_assert(std::is_floating_point<T>::value, "isZero can only be used with floating-point types.");
    return x == 0;
}
#pragma warning(pop)
#pragma GCC diagnostic pop

} // end namespace GB

#endif // GBUTIL_IMPL_H