#include "Reader.h"

#include <sstream>

//--------------------------------------------------------------------------------
Reader::Reader()
{
}

//--------------------------------------------------------------------------------
Reader::~Reader()
{
}

//--------------------------------------------------------------------------------
double Reader::parseDouble(const CVSIterator& it) const
{
  std::string s = it->c_str();
  std::istringstream os(s);
  double d;
  os >> d;
  return d;
}

//--------------------------------------------------------------------------------
int Reader::parseInteger(const CVSIterator& it) const
{
  return atoi((*it).c_str());
}

//--------------------------------------------------------------------------------
std::vector<std::string> Reader::split(const std::string& text, char sep, bool multiple)
{
  std::vector<std::string> result;
  std::string part;
  bool set = false;
  bool prev = false;

  CSIterator it = text.cbegin();

  while (it != text.cend())
  {
    if (*it == sep)
    {
      if (set || (!multiple && prev))
      {
        result.push_back(part);
        part = "";

        set = false;
        prev = false;
      }
      else
      {
        prev = true;
      }
    }
    else
    {
      // skip tabs
      if (*it != '\t')
      {
        set = true;
        part += *it;
      }
    }

    ++it;
  }

  // Add the remaining part
  if (part.size() > 0)
    result.push_back(part);

  return result;
}