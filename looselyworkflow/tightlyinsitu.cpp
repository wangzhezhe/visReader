#include <iostream>
#include <thallium/serialization/stl/string.hpp>
#include <thallium.hpp>

namespace tl = thallium;

// The tightly part should be like a client and send data into the server(loosely part)
int main(int argc, char** argv) {
    if(argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <address> <provider_id>" << std::endl;
        exit(0);
    }

    //let rank 0 to detect the address of vis server

    //broad cast store all server in list


    tl::engine myEngine("tcp", THALLIUM_CLIENT_MODE);
    tl::remote_procedure sum   = myEngine.define("sum");
    tl::remote_procedure prod  = myEngine.define("prod");
    tl::remote_procedure hello = myEngine.define("hello").disable_response();
    tl::remote_procedure print = myEngine.define("print").disable_response();
    tl::endpoint server = myEngine.lookup(argv[1]);
    uint16_t provider_id = atoi(argv[2]);
    tl::provider_handle ph(server, provider_id);
    std::cout << "identity: " << ph.get_identity() << std::endl;
    int ret = sum.on(ph)(42,63);
    std::cout << "(sum) Server answered " << ret << std::endl;
    ret = prod.on(ph)(42,63);
    std::cout << "(prod) Server answered " << ret << std::endl;
    std::string name("Matthieu");
    hello.on(ph)(name);
    std::cout << "Done sending hello RPC, no response expected" << std::endl;
    print.on(ph)(name);
    std::cout << "Done sending print RPC, no response expected" << std::endl;

    return 0;
}