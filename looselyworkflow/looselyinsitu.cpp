#include <iostream>
#include <thallium.hpp>
#include <thallium/serialization/stl/string.hpp>
#include <thallium/serialization/stl/vector.hpp>
#include <spdlog/spdlog.h>
#include <mpi.h>
#include <fstream>

// The loosely part can be a visualization service
// it need to be started when all the tightly insitu (client) is well started
namespace tl = thallium;

int globalRank = 0;
int globalProc = 0;
std::vector<std::string> globalAddrList;

class my_sum_provider : public tl::provider<my_sum_provider>
{

private:
    // for testing
    tl::auto_remote_procedure m_prod_rpc;
    tl::auto_remote_procedure m_sum_rpc;
    tl::auto_remote_procedure m_hello_rpc;
    tl::auto_remote_procedure m_print_rpc;
    tl::auto_remote_procedure m_getaddrlist_rpc;

    // for rank 0, get all addresses

    void prod(const tl::request &req, int x, int y)
    {
        std::cout << "Computing " << x << "*" << y << std::endl;
        req.respond(x * y);
    }

    int sum(int x, int y) const
    {
        std::cout << "Computing " << x << "+" << y << std::endl;
        return x + y;
    }

    void hello(const std::string &name)
    {
        std::cout << "Hello, " << name << ", from " << identity() << std::endl;
    }

    int print(const std::string &word)
    {
        std::cout << "Printing " << word << std::endl;
        return word.size();
    }

    void getaddrlist(const tl::request &req, int x, int y)
    {
        std::cout << "Computing " << x << "*" << y << std::endl;
        req.respond(globalAddrList);
    }


public:
    my_sum_provider(const tl::engine &e, uint16_t provider_id)
        : tl::provider<my_sum_provider>(e, provider_id, "myprovider"), 
        m_prod_rpc{define("prod", &my_sum_provider::prod)}, 
        m_sum_rpc{define("sum", &my_sum_provider::sum)}, 
        m_hello_rpc{define("hello", &my_sum_provider::hello).disable_response()}, 
        m_print_rpc{define("print", &my_sum_provider::print, tl::ignore_return_value())},
        m_getaddrlist_rpc{define("getaddrlist", &my_sum_provider::getaddrlist)}
    {
    }
};

std::vector<std::string> AddrSplit(const char *s, int size, char seperatorH, char seperatorE)
{
    std::vector<std::string> result;
    typedef std::string::size_type string_size;
    string_size i = 0;
    int flag = 0;

    while (i != size)
    {
        // if flag =0 , estra str
        // if flag =1 , real content
        // int flag = 0;
        while (i != size && flag == 0)
        {
            // caculate start position

            if (s[i] == seperatorH)
            {
                flag = 1;
                ++i;
                break;
            }
            else
            {
                ++i;
            }
        }

        // caculate end position
        string_size j = i;

        while (j != size && flag == 1)
        {

            if (s[j] == seperatorE)
            {
                flag = 0;
                break;
            }

            if (flag == 1)
                ++j;
        }

        if (i != j)
        {
            char substr[100];
            std::memcpy(substr, s + i, j - i);
            substr[j - i] = '\0';
            result.push_back(std::string(substr));
            i = j;
        }
    }
    return result;
}

void gatherAddr(std::string endpoint)
{
    int msgPaddingLen = 200;
    if (endpoint.size() > msgPaddingLen)
    {
        throw std::runtime_error("current addr is longer than msgPaddingLen, reset the addr buffer to "
                                 "make it larger");
        return;
    }

    int sendLen = msgPaddingLen;
    int sendSize = sendLen * sizeof(char);
    char *sendipStr = (char *)malloc(sendSize);
    sprintf(sendipStr, "H%sE", endpoint.c_str());

    std::cout << "check send ip: " << std::string(sendipStr) << std::endl;
    int rcvLen = sendLen;
    char *rcvString = NULL;
    int rcvSize = msgPaddingLen * globalProc * sizeof(char);
    rcvString = (char *)malloc(rcvSize);
    {
        if (rcvString == NULL)
        {
            MPI_Abort(MPI_COMM_WORLD, 1);
            return;
        }
    }

    int error_code =
        MPI_Gather(sendipStr, sendLen, MPI_CHAR, rcvString, rcvLen, MPI_CHAR, 0, MPI_COMM_WORLD);
    if (error_code != MPI_SUCCESS)
    {
        std::cout << "error for rank " << globalRank << " get MPI_GatherError: " << error_code
                  << std::endl;
    }

    if (globalRank == 0)
    {
        std::vector<std::string> ipList = AddrSplit(rcvString, msgPaddingLen * globalProc, 'H', 'E');

        for (int i = 0; i < ipList.size(); i++)
        {
            // store all server addrs
            spdlog::debug("rank {} add raw data server {}", globalRank, ipList[i]);
            globalAddrList.push_back(ipList[i]);
        }
    }

    free(sendipStr);
    free(rcvString);
}

int main(int argc, char **argv)
{

    MPI_Init(NULL, NULL);

    MPI_Comm_rank(MPI_COMM_WORLD, &globalRank);
    MPI_Comm_size(MPI_COMM_WORLD, &globalProc);

    if (argc != 3)
    {
        std::cerr << "looselyinsitu <protocol>" << std::endl;
        exit(0);
    }

    uint16_t provider_id = 22;
    std::string protocol = argv[1];
    std::string logLevel = argv[2];

    if (logLevel == "info")
    {
        spdlog::set_level(spdlog::level::info);
    }
    else
    {
        spdlog::set_level(spdlog::level::debug);
    }

    tl::engine myEngine(protocol, THALLIUM_SERVER_MODE);
    std::string selfAddr = myEngine.self();
    std::cout << "Server running at address " << selfAddr
              << " with provider id " << provider_id << std::endl;

    tl::abt scope;

    if (globalRank == 0)
    {
        std::cout << "total process num: " << globalProc << std::endl;
    }



    // start provider
    auto provider = new my_sum_provider(myEngine, provider_id);
    myEngine.push_finalize_callback(provider, [provider]()
                                    { delete provider; });

    // collecting all addresses
    gatherAddr(selfAddr);

  if (globalRank == 0)
  {
    std::string confile="masterinfo.conf";
    std::ofstream confFile;
    spdlog::debug("master info file: {}", confile);
    confFile.open(confile);

    if (!confFile.is_open())
    {
      spdlog::info("Could not open file: {}", confile);
      exit(-1);
    }
    confFile << selfAddr << "\n";
    confFile.close();
  }

    myEngine.wait_for_finalize();

    std::cout << "server close" << std::endl;
    MPI_Finalize();
    return 0;
}