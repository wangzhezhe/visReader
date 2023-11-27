#include <iostream>
#include <string.h>
#include <map>
#include <list>
#include <vector>
#include <stdio.h>
#include <algorithm>
#include <vtkm/filter/flow/internal/BoundsMap.h>

#ifndef AVT_BLOCK_CHOWDER_H
#define AVT_BLOCK_CHOWDER_H


class DomainBlock
{
public:
  enum LeafBlockType {NONE, UNIFORM, NON_UNIFORM, INTERNAL, X_MIN, X_MAX, Y_MIN, Y_MAX, Z_MIN, Z_MAX};

  static std::vector<std::string> nameMap;
  static std::map<int, DomainBlock *> leafMap;
  static void CreateBlockInfo(std::vector<DomainBlock*> &v, int nDom, vtkm::filter::flow::internal::BoundsMap &it,
                              bool subdivUniform, int nX, int nY, int nZ, double pct,
                              bool skipSharedFaceStats=true);
  static int TotalNumLeaves(std::vector<DomainBlock*> &v);
  static void Dump(DomainBlock *b, std::ostream &s, int lvl, bool rawNums=false);
  static void Dump(std::vector<DomainBlock*> &v, std::ostream &s, int lvl, bool rawNums=false);
  static DomainBlock * GetBlockFromGID(std::vector<DomainBlock*> &v, int gid);
  static DomainBlock * GetLeaf_FIX_THIS(std::vector<DomainBlock*> &v, const vtkm::Vec3f &p);

  DomainBlock();
  DomainBlock(int d, vtkm::FloatDefault *b, std::string n="");
  //DomainBlock(int d, double *b, std::string n="");
  virtual ~DomainBlock();

  bool Inside(vtkm::FloatDefault *p, int &id) const;

  //bool operator< (const DomainBlock &x) const { return (dom==x.dom ? (sub < x.sub) : dom < x.dom); }

  void setBBox(vtkm::FloatDefault *b) { for(int i = 0; i < 6; i++) bbox[i]=b[i];}
  //void setBBox(double *b) { for(int i = 0; i < 6; i++) bbox[i]=b[i];}
  void setNm(const std::string &n);

  void SubdivideUniform(int nx, int ny, int nz, LeafBlockType bt=UNIFORM);
  void SubdivideFaces(int nx, int ny, int nz, vtkm::FloatDefault pct);
  void dump(std::ostream& s, int lvl, bool rawNums=false, std::string indent="") const;
  void dump(std::ostream& s, bool rawNums=false) const;
  void dumpBBox(std::ostream& s) const;

  bool InBBox(vtkm::FloatDefault *p) const;
  bool InBBox(const vtkm::Vec3f &p) const;
  bool InBBox(const vtkm::Vec3f &p, const double &tol) const;
  void GetLeaves(std::vector<DomainBlock *> &v);
  int NumLeafs() const;

  DomainBlock * GetLeafFromIndex(int idx);
  DomainBlock * GetLeaf(const vtkm::Vec3f &p);
  DomainBlock * GetRoot();
  DomainBlock * GetInternal();
  void GetBBox(vtkm::FloatDefault *bb);
  void GetExtents(vtkm::FloatDefault *ext);
  void GetExtents(vtkm::Vec3f &ext);

  bool AddBlockData(DomainBlock *dst, int numICs, int numIters, int totalNumICs);
  void UnifyData();
  int GetDataIdxFromPct(vtkm::FloatDefault p)
  {
    if (data.size() == 0) return -1;
    for (std::size_t i = 0; i < data.size(); i++)
      if (p < data[i].pct)
        return i;
      else
        p -= data[i].pct;
    return -1;
  }

  bool BlocksShareFace(DomainBlock *blk);

  class blockData
  {
  public:
    blockData(DomainBlock *b, int n, int ni, int tni):blk(b),numICs(n),numIters(ni),
                                                      totalNumICs(tni),pct(0.0),avgIt(0.0) {}
    blockData() {numICs=numIters=totalNumICs=0; blk=NULL; pct=avgIt=0.0f;}

    void computeAvg() {pct = (vtkm::FloatDefault)numICs/(vtkm::FloatDefault)totalNumICs; avgIt = (vtkm::FloatDefault)numIters/(vtkm::FloatDefault)numICs;}

    DomainBlock *blk;
    int numICs, numIters, totalNumICs;
    vtkm::FloatDefault pct, avgIt;

    void dump(std::ostream& s, bool rawNums=false) const
    {
      char str[64];
      std::string n = (blk ? blk->nm : "<NO-BLK>");

      if (rawNums)
        sprintf(str, "(%s %d %d)", n.c_str(), numICs, numIters);
      else
        sprintf(str, "(%s %.2f %4.2f)", n.c_str(), pct, avgIt);
      s<<str;
    }

    static bool cmp(const blockData &a, const blockData &b) {return a.pct < b.pct;}
    static bool rcmp(const blockData &a, const blockData &b) {return a.pct > b.pct;}
  };

  //members.
  DomainBlock *parent;
  int dom, sub, gid;
  int minIt, maxIt;
  double expIt, bbox[6];
  bool skipSharedFaceStats;
  std::string nm;
  LeafBlockType leafBlockType;
  std::vector<DomainBlock *> children;
  std::vector<blockData> data;

private:
  DomainBlock * AddChild(vtkm::FloatDefault *bb, const char *n);
  int SetIDs(int id);
  void DoSubdivideFaces(int nx, int ny, int nz, vtkm::FloatDefault pct);
  void DoSubdivideUniform(int nx, int ny, int nz);
};

#endif
