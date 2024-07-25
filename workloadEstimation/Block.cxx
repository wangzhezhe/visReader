#include <limits>
#include "Block.h"

std::vector<std::string> DomainBlock::nameMap;
std::map<int, DomainBlock *> DomainBlock::leafMap;

DomainBlock::DomainBlock()
{
    gid = -1;
    dom = -1;
    sub = -1;
    parent = NULL;
    skipSharedFaceStats = false;
    for(int i=0; i<6; i++)
        bbox[i]=0.0f;
    setNm(std::string(""));
    leafBlockType = NONE;
}

DomainBlock::DomainBlock(int d, vtkm::FloatDefault *b, std::string n)
{
    gid = -1;
    dom = d; // id of the domain
    sub = 0;
    parent = NULL;
    skipSharedFaceStats = false;
    setBBox(b); // record the domain info the current block
    setNm(n);
    leafBlockType = NONE;
}

void DomainBlock::PrintInfo(std::string indent){
std::cout << indent << " info domid " <<this->dom <<", subid "<< this->sub << ", gid " << this->gid <<", name "<< this->nm << ", LeafBlockType " << int(this->leafBlockType) << std::endl;
//when gid it -1 it is dummy node, which is supposed to be decomposed into sub-blocks
//subid is unique in each block
//gid is unique for all domain containing multiple blocks
std::cout << indent << " block info " <<this->bbox[0] <<"," << this->bbox[1] <<"," << this->bbox[2]<<"," << this->bbox[3]<<"," << this->bbox[4] <<"," << this->bbox[5] << std::endl;
std::cout << indent << " children size " << (this->children).size() << std::endl;
for(int i=0;i<(this->children).size();i++){
    std::string indentNext=indent+"    ";
    this->children[i]->PrintInfo(indentNext);
}
}
/*
DomainBlock::DomainBlock(int d, double *b, std::string n)
{
    gid = -1;
    dom = d;
    sub = 0;
    parent = NULL;
    skipSharedFaceStats = false;
    setBBox(b);
    setNm(n);
    leafBlockType = NONE;
}
*/
DomainBlock::~DomainBlock()
{
    for(int i=0; i<children.size(); i++)
        delete children[i];
    children.clear();
}

void
DomainBlock::setNm(const std::string &n)
{
    nm = "";
    if (n.empty() && dom >= 0)
    {
        char t[32];
        snprintf(t, 32 ,"%d", dom);
        nm = t;
    }
    else
        nm = n;
}

DomainBlock *
DomainBlock::AddChild(vtkm::FloatDefault *bb, const char *n)
{
    char t[32];
    snprintf(t, 32, "%s:%s", nm.c_str(), n);
    DomainBlock *blk = new DomainBlock(dom, bb, t);
    blk->parent = this;
    //add new blk into the children list attached with current node
    children.push_back(blk);
    return blk;
}

void
DomainBlock::GetLeaves(std::vector<DomainBlock *> &v)
{
    if (children.empty())
        v.push_back(this);
    else
        for (std::size_t i = 0; i < children.size(); i++)
            children[i]->GetLeaves(v);
}

void
DomainBlock::GetBBox(vtkm::FloatDefault *bb)
{
    for (std::size_t i = 0; i < 6; i++)
        bb[i] = bbox[i];
}

void
DomainBlock::GetExtents(vtkm::FloatDefault *ext)
{
    ext[0] = bbox[1]-bbox[0];
    ext[1] = bbox[3]-bbox[2];
    ext[2] = bbox[5]-bbox[4];
}

void
DomainBlock::GetExtents(vtkm::Vec3f &ext)
{
    ext[0] = bbox[1]-bbox[0];
    ext[1] = bbox[3]-bbox[2];
    ext[2] = bbox[5]-bbox[4];
}

int
DomainBlock::NumLeafs() const
{
    //the current one is leaf
    if (children.empty())
        return 1;

    int n = 0;
    //return number of leafs that the current one can hold
    for (std::size_t i = 0; i < children.size(); i++)
        n += children[i]->NumLeafs();
    return n;
}

DomainBlock *
DomainBlock::GetLeafFromIndex(int idx)
{
    std::vector<DomainBlock*> leaves;
    GetLeaves(leaves);
    return leaves[idx];
}

bool
DomainBlock::InBBox(vtkm::FloatDefault *p) const
{
    if (p[0] >= bbox[0] && p[0] <= bbox[1])
        if (p[1] >= bbox[2] && p[1] <= bbox[3])
            if (p[2] >= bbox[4] && p[2] <= bbox[5])
                return true;
    return false;
}
bool
DomainBlock::InBBox(const vtkm::Vec3f &p) const
{
    vtkm::FloatDefault pt[3] = {p[0],p[1],p[2]};
    return InBBox(pt);
}

bool
DomainBlock::InBBox(const vtkm::Vec3f &p, const double &tol) const
{
    double bb[6] = {bbox[0]-tol, bbox[1]+tol,
                    bbox[2]-tol, bbox[3]+tol,
                    bbox[4]-tol, bbox[5]+tol};
    if (p[0] >= bb[0] && p[0] <= bb[1])
        if (p[1] >= bb[2] && p[1] <= bb[3])
            if (p[2] >= bb[4] && p[2] <= bb[5])
                return true;
    return false;
}

bool
DomainBlock::Inside(vtkm::FloatDefault *p, int &id) const
{
    if (!InBBox(p))
        return false;

    //hit a node.
    int sz = children.size();
    if (sz == 0)
    {
        id = sub;
        return true;
    }

    for (std::size_t i = 0; i < sz; i++)
        if (children[i]->Inside(p, id))
            return true;

    return false;
}

DomainBlock *
DomainBlock::GetInternal()
{
  auto root = this->GetRoot();
  for (auto child : root->children)
    if (child->leafBlockType == INTERNAL)
      return child;

  return nullptr;
}

DomainBlock *
DomainBlock::GetRoot()
{
  DomainBlock* root = this;
  while (root->parent != nullptr)
    root = root->parent;

  return root;
}

DomainBlock *
DomainBlock::GetLeaf(const vtkm::Vec3f &p)
{
    std::vector<DomainBlock*> leaves;
    GetLeaves(leaves);
    for (std::size_t i = 0; i < leaves.size(); i++)
        if (leaves[i]->InBBox(p))
            return leaves[i];

    double tol = 1e-4;
    for (std::size_t i = 0; i < leaves.size(); i++)
        if (leaves[i]->InBBox(p, tol))
            return leaves[i];

    return NULL;
}

void
DomainBlock::dump(std::ostream &s, bool rawNums) const
{
    s<<"["<<nm<<"] ";
    //dumpBBox(s);
    s<<"[";
    for (std::size_t i = 0; i < data.size(); i++)
    {
        if ( i > 0)
            s<<" ";
        data[i].dump(s, rawNums);
    }
    s<<"]";
}

void
DomainBlock::dumpBBox(std::ostream& s) const
{
    s<<"{("<<bbox[0]<<","<<bbox[1]<<") ";
    s<<"("<<bbox[2]<<","<<bbox[3]<<") ";
    s<<"("<<bbox[4]<<","<<bbox[5]<<")} ";
}
//? what does the lvl and rawNums represent here??? level?
//only output first two level?
void
DomainBlock::dump(std::ostream &s, int lvl, bool rawNums, std::string indent) const
{
    s<<indent;
    dump(s, rawNums);
    s<<std::endl;
    if (lvl > 0)
    {
        for (std::size_t i=0; i<children.size(); i++)
            children[i]->dump(s, lvl-1, rawNums, indent+" ");
    }
}

static void
inflate(vtkm::FloatDefault *x, vtkm::FloatDefault xi)
{
    x[0] -= xi;
    x[1] += xi;
    x[2] -= xi;
    x[3] += xi;
    x[4] -= xi;
    x[5] += xi;
}

static void
mkCorners(std::vector<vtkm::Vec3f> &pts, vtkm::FloatDefault *bb)
{
    for (std::size_t i = 0; i < 2; i++)
        for (std::size_t j = 2; j < 4; j++)
            for (std::size_t k = 4; k < 6; k++)
                pts.push_back(vtkm::Vec3f(bb[i], bb[j], bb[k]));
}

static bool
ptInside(const vtkm::Vec3f &p, vtkm::FloatDefault *bb)
{
    return (p[0] >= bb[0] && p[0] <= bb[1] &&
            p[1] >= bb[2] && p[1] <= bb[3] &&
            p[2] >= bb[4] && p[2] <= bb[5]);
}

bool
DomainBlock::BlocksShareFace(DomainBlock *blk)
{
    if (blk == NULL)
        return false;

    if (leafBlockType == INTERNAL || blk->leafBlockType == INTERNAL)
        return false;
    if (children.size() > 0 || blk->children.size() > 0)
        return false;
    if (dom == blk->dom)
        return false;

    if (1)
    {
        vtkm::FloatDefault bb0[6], bb1[6], bb0i[6], bb1i[6];
        GetBBox(bb0);
        blk->GetBBox(bb1);
        memcpy(bb0i, bb0, 6*sizeof(vtkm::FloatDefault));
        memcpy(bb1i, bb1, 6*sizeof(vtkm::FloatDefault));
        vtkm::FloatDefault eps = std::numeric_limits<vtkm::FloatDefault>::epsilon() * 10.0f;
        inflate(bb0i, eps);
        inflate(bb1i, eps);

        std::vector<vtkm::Vec3f> pts0, pts1;
        mkCorners(pts0, bb0);
        mkCorners(pts1, bb1);

        for (std::size_t i = 0; i < pts0.size(); i++)
            if (ptInside(pts0[i], bb1i))
                return true;

        for (std::size_t i = 0; i < pts1.size(); i++)
            if (ptInside(pts1[i], bb0i))
                return true;

        return false;
    }
    else
    {
        LeafBlockType b0 = leafBlockType, b1 = blk->leafBlockType;

        //cout<<"Blocks Share face? ("<<nm<<","<<blk->nm<<")"<<std::endl;

        if (b0 == X_MIN && b1 == X_MAX) return true;
        if (b0 == Y_MIN && b1 == Y_MAX) return true;
        if (b0 == Z_MIN && b1 == Z_MAX) return true;

        if (b0 == X_MAX && b1 == X_MIN) return true;
        if (b0 == Y_MAX && b1 == Y_MIN) return true;
        if (b0 == Z_MAX && b1 == Z_MIN) return true;

        return false;
    }
}

bool
DomainBlock::AddBlockData(DomainBlock *dstBlk, int numICs, int numIters, int totalNumICs)
{
    if (dstBlk == NULL)
        return false;
    /*
    if (!skipSharedFaceStats)
    {
        if (BlocksShareFace(dstBlk))
        {
            cout<<"ERROR: "<<nm<<" --> "<<dstBlk->nm<<endl;
            cout<<"  ";
            dumpBBox(cout);
            cout<<"  ";
            dstBlk->dumpBBox(cout);
            cout<<std::endl;
        }
    }
    */

    if (skipSharedFaceStats && BlocksShareFace(dstBlk))
        return false;

    bool foundIt = false;
    for (std::size_t i = 0; i < data.size(); i++)
        if (data[i].blk == dstBlk)
        {
            data[i].numICs += numICs;
            data[i].numIters += numIters;
            data[i].totalNumICs += totalNumICs;
            foundIt = true;
        }

    if (!foundIt)
    {
        blockData b(dstBlk, numICs, numIters, totalNumICs);
        data.push_back(b);
    }

    if (parent && dstBlk->parent)
    {
        //Arg... Interior block is off by a level. No children, so skip a level.
        if (nm[nm.size()-1] == 'I')
            parent->AddBlockData(dstBlk->parent->parent, numICs, numIters, totalNumICs);
        else
            parent->AddBlockData(dstBlk->parent, numICs, numIters, totalNumICs);
    }
    return true;
}

void
DomainBlock::UnifyData()
{
    expIt = 0;
    double totNumIters = 0;
    int totICs = 0;
    for (std::size_t i = 0; i < data.size(); i++)
    {
        if(!i)
        {
            minIt = data[i].numIters;
            maxIt = data[i].numIters;
        }
        else
        {
            if(data[i].numIters < minIt)
                minIt = data[i].numIters;
            else if(data[i].numIters > maxIt)
                maxIt = data[i].numIters;
        }
        totICs += data[i].numICs;
        totNumIters += data[i].numIters;
    }
    for (std::size_t i = 0; i < data.size(); i++)
        data[i].totalNumICs = totICs;

    expIt = totNumIters/totICs;
    //cout<<minIt<<" "<<maxIt<<" "<<expIt<<std::endl;
    for (std::size_t i = 0; i < data.size(); i++)
        data[i].computeAvg();
    sort(data.begin(), data.end(), blockData::rcmp);

    for (std::size_t i = 0; i < children.size(); i++)
        children[i]->UnifyData();
}

int
DomainBlock::SetIDs(int id)
{
    sub = -1;

    int sz = children.size();
    if (sz == 0)
    {
        sub = id;
        id++;
    }
    else
    {
        for (std::size_t i = 0; i < sz; i++)
            id = children[i]->SetIDs(id);
    }

    return id;
}

//subdivide this block in an uniform way
//each dimention is nx ny and nz
void
DomainBlock::DoSubdivideUniform(int nx, int ny, int nz)
{
    //Nothing to subdivide....
    if (nx == 1 && ny == 1 && nz == 1)
        return;

    vtkm::FloatDefault b[6];
    char n[64];

    vtkm::FloatDefault dx = (bbox[1]-bbox[0]) / (vtkm::FloatDefault)nx;
    vtkm::FloatDefault dy = (bbox[3]-bbox[2]) / (vtkm::FloatDefault)ny;
    vtkm::FloatDefault dz = (bbox[5]-bbox[4]) / (vtkm::FloatDefault)nz;

    for (std::size_t i = 0; i < nx; i++)
        for (std::size_t j = 0; j < ny; j++)
            for (std::size_t k = 0; k < nz; k++)
            {
                b[0] = bbox[0] + i*dx;
                b[1] = b[0] + dx;
                b[2] = bbox[2] + j*dy;
                b[3] = b[2] + dy;
                b[4] = bbox[4] + k*dz;
                b[5] = b[4] + dz;

                snprintf(n, 64, "%ld%ld%ld", i,j,k);
                //Create a new data block meta info
                //add this meta info to the list of children within current datadomain
                DomainBlock *blk = AddChild(b, n);

                //Propagate the block type to the children.
                //the block type of children is same with parent?
                blk->leafBlockType = leafBlockType;
            }
}

void
DomainBlock::DoSubdivideFaces(int nx, int ny, int nz, vtkm::FloatDefault pct)
{
    vtkm::FloatDefault px = pct*(bbox[1]-bbox[0]);
    vtkm::FloatDefault py = pct*(bbox[3]-bbox[2]);
    vtkm::FloatDefault pz = pct*(bbox[5]-bbox[4]);
    DomainBlock *blk;
    // There are five type of subdomain, internal, x, X, y, Y, z, Z
    // In the bounudry region, we then use the nx ny nz to continue divide them into nx*ny*nz small blocks
    // The internal region does not need to do the subdivision
    //Make the interior.
    vtkm::FloatDefault b[6] = {bbox[0]+px, bbox[1]-px,
                  bbox[2]+py, bbox[3]-py,
                  bbox[4]+pz, bbox[5]-pz};
    blk = AddChild(b, "I");
    blk->leafBlockType = INTERNAL;

    //Make the Xmin:
    vtkm::FloatDefault bbx[6] = {bbox[0], bbox[0]+px,
                    bbox[2], bbox[3],
                    bbox[4], bbox[5]};
    blk = AddChild(bbx, "x");
    blk->leafBlockType = X_MIN;
    blk->SubdivideUniform(1, ny, nz, blk->leafBlockType);

    //Make the Xmax:
    vtkm::FloatDefault bbX[6] = {bbox[1]-px, bbox[1],
                    bbox[2], bbox[3],
                    bbox[4], bbox[5]};
    blk = AddChild(bbX, "X");
    blk->leafBlockType = X_MAX;
    blk->SubdivideUniform(1, ny, nz, blk->leafBlockType);

    //Make the Ymin:
    vtkm::FloatDefault bby[6] = {bbox[0]+px, bbox[1]-px,
                    bbox[2], bbox[2]+py,
                    bbox[4], bbox[5]};
    blk = AddChild(bby, "y");
    blk->leafBlockType = Y_MIN;
    blk->SubdivideUniform(nx, 1, nz, blk->leafBlockType);

    //Make the Ymax:
    vtkm::FloatDefault bbY[6] = {bbox[0]+px, bbox[1]-px,
                    bbox[3]-py, bbox[3],
                    bbox[4], bbox[5]};
    blk = AddChild(bbY, "Y");
    blk->leafBlockType = Y_MAX;
    blk->SubdivideUniform(nx, 1, nz, blk->leafBlockType);

    //Make the Zmin:
    vtkm::FloatDefault bbz[6] = {bbox[0]+px, bbox[1]-px,
                    bbox[2]+py, bbox[3]-py,
                    bbox[4], bbox[4]+pz};
    blk = AddChild(bbz, "z");
    blk->leafBlockType = Z_MIN;
    blk->SubdivideUniform(nx, ny, 1, blk->leafBlockType);

    //Make the Zmax:
    vtkm::FloatDefault bbZ[6] = {bbox[0]+px, bbox[1]-px,
                                 bbox[2]+py, bbox[3]-py,
                                 bbox[5]-pz, bbox[5]};
    blk = AddChild(bbZ, "Z");
    blk->leafBlockType = Z_MAX;
    blk->SubdivideUniform(nx, ny, 1, blk->leafBlockType);
}

void
DomainBlock::SubdivideUniform(int nx, int ny, int nz, LeafBlockType bt)
{
    leafBlockType = bt;
    DoSubdivideUniform(nx, ny, nz);
    SetIDs(0);
}

void
DomainBlock::SubdivideFaces(int nx, int ny, int nz, vtkm::FloatDefault pct)
{
    leafBlockType = NON_UNIFORM;
    DoSubdivideFaces(nx, ny, nz, pct);
    SetIDs(0);
}

/*
inline std::ostream&
operator<<(std::ostream &out, const DomainBlock &b)
{
    out<<"["<<b.dom;
    if (b.sub != -1)
        out<<"."<<b.sub;
    if (b.nm != "")
        out<<":"<<b.nm;
    out<<"] ";

    out<<"{("<<b.bbox[0]<<","<<b.bbox[1]<<") ";
    out<<"("<<b.bbox[2]<<","<<b.bbox[3]<<") ";
    out<<"("<<b.bbox[4]<<","<<b.bbox[5]<<")} ";
    //out<<b.data;
    return out;
}
*/

void
DomainBlock::Dump(DomainBlock *b, std::ostream &s, int lvl, bool rawNums)
{
    b->dump(s,lvl, rawNums);
}

void
DomainBlock::Dump(std::vector<DomainBlock*> &v, std::ostream &s, int lvl, bool rawNums)
{
    for (std::size_t i = 0; i < v.size(); i++)
        DomainBlock::Dump(v[i], s, lvl, rawNums);
}

int
DomainBlock::TotalNumLeaves(const std::vector<DomainBlock*> &v)
{
    int N = 0;
    for (std::size_t i = 0; i < v.size(); i++)
        N += v[i]->NumLeafs();
    return N;
}

void
DomainBlock::CreateBlockInfo(std::vector<DomainBlock*> &v, int nDom, vtkm::filter::flow::internal::BoundsMap &it,
                             bool subdivUniform, int nX, int nY, int nZ, double pct,
                             bool _skipSharedFaceStats)
{
    double bbox[6];
    v.resize(nDom);

    int totLeafs = 0;
    for (int i = 0; i < nDom; i++)
    {
      //this api is removed?
      //get spatial bounds for each domain
      auto bounds = it.GetBlockBounds(i);
      bbox[0] = bounds.X.Min;
      bbox[1] = bounds.X.Max;
      bbox[2] = bounds.Y.Min;
      bbox[3] = bounds.Y.Max;
      bbox[4] = bounds.Z.Min;
      bbox[5] = bounds.Z.Max;

      //it->GetElementExtents(i, bbox);
      //init a struct for new data domain
      //ser the id and bbx for the current domain
      v[i] = new DomainBlock(i,bbox);
      v[i]->skipSharedFaceStats = _skipSharedFaceStats;
      if (subdivUniform)
        // when need to sub divide the uniform
        // subdivide curernt domain in distributed way
        // nX, nY nZ represents number of subdomains in each dim
        // this condition is actually set as false in StreamlinMPI
        v[i]->SubdivideUniform(nX, nY, nZ);
      else
        // this will subdivide the domain into the internal region and the bounudary region
        // pct represents the ratio of the border width to the entire width
        v[i]->SubdivideFaces(nX, nY, nZ, pct);
      
      totLeafs += v[i]->NumLeafs(); // compute how many leaf for the current block

    }

    //Set the GIDs, and the name map.
    nameMap.resize(totLeafs);
    int id = 0;
    for (std::size_t i = 0; i < v.size(); i++)
    {
        std::vector<DomainBlock*> leaves;
        //this is a recursive call
        //all leafs of current nodes are stored into the vector
        v[i]->GetLeaves(leaves);
        v[i]->skipSharedFaceStats = _skipSharedFaceStats;
        for (std::size_t j = 0; j < leaves.size(); j++)
        {
            //this is the global id including all leafes
            //assign global id to each leaf
            leaves[j]->gid = id;
            leaves[j]->skipSharedFaceStats = _skipSharedFaceStats;
            id++;
            //map leaf id to its name
            //map leaf id to its meatdata
            //during the process of dividing faces
            //the name of the leaf can be I x X y Y z Z etc
            nameMap[leaves[j]->gid] = leaves[j]->nm;
            leafMap[leaves[j]->gid] = leaves[j];
        }
    }

    //check the all block info as needed
    // for(int i=0;i<nDom;i++){
    //     std::cout << "------info for domain " << i << std::endl;
    //     v[i]->PrintInfo("");
    // }
      
}

DomainBlock *
DomainBlock::GetBlockFromGID(const std::vector<DomainBlock*>& vtkmNotUsed(v), int gid)
{
    std::map<int,DomainBlock*>::iterator it = leafMap.find(gid);
    if (it == leafMap.end())
        return NULL;

    return it->second;

    /*
    for (std::size_t i = 0; i < v.size(); i++)
    {
        std::vector<DomainBlock*> leaves;
        v[i]->GetLeaves(leaves);
        int n = leaves.size();
        if ( n > 0 && gid >= leaves[0]->gid && gid <= leaves[n-1]->gid)
            for (std::size_t j = 0; j < n; j++)
                if (leaves[j]->gid == gid)
                    return leaves[j];
    }
    return NULL;
    */
}

DomainBlock *
DomainBlock::GetLeaf_FIX_THIS(std::vector<DomainBlock*> &v, const vtkm::Vec3f &p)
{
    DomainBlock *blk = NULL;
    for (std::size_t i = 0; !blk && i < v.size(); i++)
        blk = v[i]->GetLeaf(p);

    return blk;
}
