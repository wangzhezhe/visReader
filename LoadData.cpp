#include "LoadData.hpp"
#include "filter.hpp"
#include "AssignStrategy.hpp"
#include <vtkm/cont/ArrayCopy.h>
#include <vtkm/io/VTKDataSetReader.h>

#include <vector>

vtkm::cont::ArrayHandle<vtkm::Vec3f> CreateConstantVectorField(vtkm::Id num, const vtkm::Vec3f &vec)
{
  vtkm::cont::ArrayHandleConstant<vtkm::Vec3f> vecConst;
  vecConst = vtkm::cont::make_ArrayHandleConstant(vec, num);

  vtkm::cont::ArrayHandle<vtkm::Vec3f> vecField;
  vtkm::cont::ArrayCopy(vecConst, vecField);
  return vecField;
}

static std::vector<int> getIntList(std::string str)
{
  std::vector<int> intList;
  int leftIdx = 0;
  int rightIdx = 0;
  while (rightIdx < str.size())
  {
    while (rightIdx < str.size() && str[rightIdx] != ' ')
      rightIdx++;
    intList.push_back(std::stoi(str.substr(leftIdx, rightIdx - leftIdx)));
    leftIdx = rightIdx + 1;
    rightIdx = rightIdx + 1;
  }
  return intList;
}

static std::vector<int> assignBlocksToRank(int totalBlocks, int nRanks, int rank)
{
  std::vector<int> assignedBlocks;

  int nPer = totalBlocks / nRanks;
  int b0 = rank * nPer, b1 = (rank + 1) * nPer;

  if (rank == (nRanks - 1))
    b1 = totalBlocks;

  for (int i = 0; i < totalBlocks; i++)
  {
    if (i >= b0 && i < b1)
    {
      assignedBlocks.push_back(i);
    }
  }

  return assignedBlocks;
}

static std::vector<int> assignBlocksToRankOffset(int totalBlocks, int nRanks, int rank, int offset)
{

  std::vector<int> assignedBlocks;

  int nPer = totalBlocks / nRanks;
  int b0 = rank * nPer, b1 = (rank + 1) * nPer;

  // TODO, add better mechanism for offset setting
  // exp
  // only work for 2 ranks case
  if (rank == (nRanks - 1))
    b1 = totalBlocks;

  for (int i = 0; i < totalBlocks; i++)
  {
    if (rank == 0)
    {
      if (i >= 2 && i <= 5)
      {
        assignedBlocks.push_back(i);
      }
    }
    else
    {
      if (i >= 2 && i <= 5)
      {
        continue;
      }
      else
      {
        assignedBlocks.push_back(i);
      }
    }
  }

  return assignedBlocks;
}

// round roubin / clock wise assign
static std::vector<int> assignBlocksToRankRR(int totalBlocks, int nRanks, int rank)
{
  std::vector<int> assignedBlocks;

  for (int i = 0; i < totalBlocks; i++)
  {
    if (i % nRanks == rank)
    {
      assignedBlocks.push_back(i);
    }
  }

  return assignedBlocks;
}

// assign by file
static std::vector<int> assignByFile(const std::string& assignFileName, int totalBlocks, int nRanks, int rank)
{
  // use the default name
  std::ifstream infile(assignFileName);
  std::vector<std::vector<int>> allBlockList;
  // load the file
  std::string line;
  std::vector<int> idlist;
  // assume the input file is correct
  int blockNum = 0;
  while (std::getline(infile, line))
  {

    // split by space
    idlist = getIntList(line);
    // std::cout << idlist.size() << std::endl;
    allBlockList.push_back(idlist);
    blockNum += idlist.size();
  }

  // check the file content
  // if (allBlockList.size() != nRanks)
  //{
  //  throw std::runtime_error("allBlockList size not equals to nRanks");
  //}
  // if (blockNum != totalBlocks)
  //{
  //  throw std::runtime_error("blockNum size not equals to totalBlocks");
  //}

  return allBlockList[rank];
}


void LoadData(const std::string& visitfileName,
              const AssignStrategy& assignStrategy,
              const std::string& assignFileName,
              std::vector<vtkm::cont::DataSet> &dataSets,
              std::vector<int> &blockIDList,
              int rank,
              int nRanks)
{
  std::string buff;
  std::ifstream is;
  is.open(visitfileName);
  std::cout << "Opening: " << visitfileName << std::endl;
  if (!is)
  {
    std::cout << "File not found! : " << visitfileName << std::endl;
    throw "unknown file: " + visitfileName;
  }

  auto p0 = visitfileName.rfind(".visit");
  if (p0 == std::string::npos)
    throw "Only .visit files are supported.";
  auto tmp = visitfileName.substr(0, p0);
  auto p1 = tmp.rfind("/");
  auto dir = tmp.substr(0, p1);

  std::getline(is, buff);
  auto numBlocks = std::stoi(buff.substr(buff.find("!NBLOCKS ") + 9, buff.size()));
  if (rank == 0)
    std::cout << "numBlocks= " << numBlocks << std::endl;

  if (assignStrategy == AssignStrategy::CONTINUOUS)
  {
    blockIDList = assignBlocksToRank(numBlocks, nRanks, rank);
  }
  else if (assignStrategy == AssignStrategy::ROUNDROUBIN)
  {
    blockIDList = assignBlocksToRankRR(numBlocks, nRanks, rank);
  }
  else if (assignStrategy == AssignStrategy::ASSIGNFILE)
  {
    blockIDList = assignByFile(assignFileName, numBlocks, nRanks, rank);
  }
  else
  {
    throw std::runtime_error("unsupported assignStrategy");
  }
  //using num ranks instead of numblocks
  //since the rank may assign an empty block
  for (int i = 0; i < nRanks; i++)
  {
    std::getline(is, buff);

    auto it = std::find(blockIDList.begin(), blockIDList.end(), i);
    if (it != blockIDList.end())
    {
      std::cout << "rank " << rank << " assign blockid " << i << std::endl;
      vtkm::cont::DataSet ds;
      std::string vtkFile = dir + "/" + buff;
      vtkm::io::VTKDataSetReader reader(vtkFile);
      ds = reader.ReadDataSet();

      std::vector<vtkm::Vec3f> vecs = {{1,0,0}, {0,1,0}, {1,1,0}};
      std::vector<std::string> vecNms = {"vecX", "vecY", "vecXY"};
      for (std::size_t i = 0; i < vecs.size(); i++)
      {
        auto vecField = CreateConstantVectorField(ds.GetNumberOfPoints(), vecs[i]);
        ds.AddPointField(vecNms[i], vecField);
      }
      if (FILTER::GLOBAL_BLOCK_DUPLICATE)
      {
        // put the id that needs to assign seeds at first position
        // this can only work properly when number of blocks equals
        // to the number of ranks, and there are duplicated blocks for ranks
        if (rank == i)
        {
          FILTER::GLOBAL_BLOCKIDS.insert(FILTER::GLOBAL_BLOCKIDS.begin(), i);
          dataSets.insert(dataSets.begin(), ds);
        }
        else
        {
          FILTER::GLOBAL_BLOCKIDS.push_back(i);
          dataSets.push_back(ds);
        }
      }
      else
      {
        // normal push back
        dataSets.push_back(ds);
      }
    }
  }
}
