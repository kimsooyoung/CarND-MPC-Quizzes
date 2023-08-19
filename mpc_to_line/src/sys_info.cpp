#include <iostream>
#include <sys/utsname.h>

int main() {
    struct utsname systemInfo;

    uname(&systemInfo);

    // if (uname(&systemInfo) != 0) {
    //     std::cerr << "Error getting system information." << std::endl;
    //     return 1;
    // }

    std::cout << "System name: " << systemInfo.sysname << std::endl;
    std::cout << "Node name: " << systemInfo.nodename << std::endl;
    std::cout << "Release: " << systemInfo.release << std::endl;
    std::cout << "Version: " << systemInfo.version << std::endl;
    std::cout << "Machine: " << systemInfo.machine << std::endl;

    return 0;
}