#ifndef LOAD_DATA_HPP
#define LOAD_DATA_HPP

#include <vtkm/cont/DataSet.h>
#include <vector>
#include "AssignStrategy.hpp"

void LoadData(const std::string& visitfileName,
              const AssignStrategy& assignStrategy,
              const std::string& assignFileName,
              std::vector<vtkm::cont::DataSet> &dataSets,
              std::vector<int> &blockIDList,
              int rank,
              int nRanks);

#endif
