#include "apa_hic.h"

#include <exception>
#include <iostream>

int main(int argc, char* argv[]) {
    try {
        Options options = parseOptions(argc, argv);
        return runAPA(options);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}
