#pragma once

#include <vector>

typedef std::string::const_iterator CSIterator;
typedef std::vector<std::string>::const_iterator CVSIterator;

class Reader
{
public:
  Reader();
  ~Reader();

  double parseDouble(const CVSIterator& iterator) const;
  int parseInteger(const CVSIterator& iterator) const;

  std::vector<std::string> split(const std::string& text, char sep, bool multiple = true);
};