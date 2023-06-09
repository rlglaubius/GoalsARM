#ifndef GBUTIL_H
#define GBUTIL_H

#include <sstream>
#include <string>
#include <vector>

namespace GB {

// Split a delimited string into a vector of strings. For example,
// split("abc,def", ',') returns a vector containing "abc" and "def".
std::vector<std::string> split(const std::string &str, const char delimiter);

// Return TRUE if the string str contains the substring sub, FALSE otherwise
bool contains(const std::string &str, const std::string &sub);

} // end namespace GB

#endif // GBUTIL_H
