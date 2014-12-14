#pragma once

#include <vector>

typedef std::_String_const_iterator<std::_String_val<std::_Simple_types<char>>> CIIterator;

class ObjReader
{
public:
  static std::vector<std::string> splitLine(const std::string& line, char sep, bool multiple = true);
};