#include <iostream>
#include <fstream>
#include <thallium/serialization/stl/string.hpp>
#include <thallium/serialization/stl/vector.hpp>
#include <thallium.hpp>
#include <mpi.h>

namespace tl = thallium;

int globalRank = 0;
int globalProc = 0;
std::vector<std::string> globalAddrList;

std::string loadMasterAddr(std::string masterConfigFile)
{

    std::ifstream infile(masterConfigFile);
    std::string content = "";
    std::getline(infile, content);
    // spdlog::debug("load master server conf {}, content -{}-", masterConfigFile,content);
    if (content.compare("") == 0)
    {
        std::getline(infile, content);
        if (content.compare("") == 0)
        {
            throw std::runtime_error("failed to load the master server\n");
        }
    }
    return content;
}

// The tightly part should be like a client and send data into the server(loosely part)
int main(int argc, char **argv)
{

    MPI_Init(NULL, NULL);

    MPI_Comm_rank(MPI_COMM_WORLD, &globalRank);
    MPI_Comm_size(MPI_COMM_WORLD, &globalProc);

    if (argc != 3)
    {
        std::cerr << "Usage: " << argv[0] << " <protocol> <masterconf>" << std::endl;
        exit(0);
    }

    // let rank 0 to detect the address of vis server
    std::string protocol = argv[1];
    std::string masterconf = argv[2];

    std::string masterAddr = loadMasterAddr(masterconf);

    std::cout << "detect server addr: " << masterAddr << std::endl;

    // broad cast store all server in list
    tl::engine myEngine(protocol, THALLIUM_CLIENT_MODE);
    // tl::remote_procedure sum   = myEngine.define("sum");
    // tl::remote_procedure prod  = myEngine.define("prod");
    tl::remote_procedure hello = myEngine.define("hello").disable_response();
    tl::remote_procedure print = myEngine.define("print").disable_response();
    tl::remote_procedure getaddrlist = myEngine.define("getaddrlist");
    tl::remote_procedure stage = myEngine.define("stage");
    tl::remote_procedure runfilter = myEngine.define("runfilter");

    // load the master server addr

    tl::endpoint masterep = myEngine.lookup(masterAddr);
    // just fix it
    uint16_t provider_id = 22;

    // create provider handler
    tl::provider_handle ph(masterep, provider_id);
    std::cout << "identity: " << ph.get_identity() << std::endl;

    if (globalRank == 0)
    {
        std::string testStr = "tightly in situ client";
        hello.on(ph)(testStr);
        std::cout << "Done sending hello RPC, no response expected" << std::endl;
    }

    // get addrlist
    std::vector<std::string> globalAddrList = getaddrlist.on(ph)();

    if (globalRank == 0)
    {
        std::cout << "addr list len:" << globalAddrList.size() << std::endl;
        for (auto addr : globalAddrList)
        {
            std::cout << "get server addr: " << addr << std::endl;
        }
    }

    // TODO call run filter for testing
    // maybe we should call the master, and master will call run for all associated ranks
    if (globalRank == 0)
    {
        // for controller rank
        // send run api for all staging service
        std::vector<tl::async_response> reqlist;
        for (auto addr : globalAddrList)
        {
            std::cout << "sent stage api to " << addr << std::endl;
            tl::endpoint addrEndPoint = myEngine.lookup(addr);
            tl::provider_handle phVisServer(addrEndPoint, provider_id);
            //use async call
            auto response=stage.on(phVisServer).async();
            reqlist.push_back(std::move(response));
        }

        //make sure all response is ok
        for (int i=0;i<reqlist.size();i++){
            int status = reqlist[i].wait();
            if(status!=0){
                std::cout << "failed to stage the data for some request " << i << std::endl;
            }
            //TODO, hangs at the reqlist here
        }
        std::cout << "Controller ok to call the stage" << std::endl;

        //clean the reqlist
        reqlist.clear();
        for (auto addr : globalAddrList)
        {
            std::cout << "sent runfilter api to " << addr << std::endl;
            tl::endpoint addrEndPoint = myEngine.lookup(addr);
            tl::provider_handle phVisServer(addrEndPoint, provider_id);
            //use async call
            int cycle=0;
            std::string field="velocity";
            auto response=runfilter.on(phVisServer).async(0,field);
            reqlist.push_back(std::move(response));
        }

        //make sure all response is ok
        for (int i=0;i<reqlist.size();i++){
            int status = reqlist[i].wait();
            if(status!=0){
                std::cout << "failed to call runfilter for some request " << i << std::endl;
            }
            //TODO, hangs at the reqlist here
        }

    }

    std::cout << "client " << globalRank << " close" << std::endl;
    MPI_Finalize();
    return 0;

    return 0;
}