#include <GBUtil.H>

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

} // end namespace GB
