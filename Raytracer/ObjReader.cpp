#include "ObjReader.h"

//--------------------------------------------------------------------------------
std::vector<std::string> ObjReader::splitLine(const std::string& line, char sep, bool multiple)
{
  std::vector<std::string> result;
  std::string part;
  bool set = false;
  bool prev = false;

  CIIterator it = line.cbegin();

  while (it != line.cend())
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
      set = true;
      part += *it;
    }

    ++it;
  }

  return result;
}