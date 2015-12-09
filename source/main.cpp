
#include <iostream>
#include <cstdlib>
#include "controller.h"


using namespace std;

int main(int Narguments, char* arguments[])
{
    if (Narguments != 2)
    {
        cerr << "Usage: platosim <inputfile>" << endl;
        exit(EXIT_FAILURE);
    }

    string inputFileName(arguments[1]);

    Controller controller(inputFileName);
    controller.run();

    return EXIT_SUCCESS;
}