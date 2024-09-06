
#ifndef METADATA_H
#define METADATA_H

#include <array>
#include <cstring>
#include <iostream>
#include <string>
#include <thallium/serialization/stl/array.hpp>
#include <thallium/serialization/stl/string.hpp>
#include <thallium/serialization/stl/vector.hpp>
#include <tuple>
#include <typeinfo>
#include <vector>

// use namespace to avoid the conflicts for type defination
// comapred with system library
namespace Insitu
{

// the meta data to index the raw data and block id
struct AddrWrapper
{
  AddrWrapper(){};
  AddrWrapper(int index, std::string addr)
    : m_index(index)
    , m_addr(addr)
  {
  }
  int m_index;
  std::string m_addr;
  ~AddrWrapper(){};

  template <typename A>
  void serialize(A& ar)
  {
    ar& m_index;
    ar& m_addr;
  }
};

}
#endif
