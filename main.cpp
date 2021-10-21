#pragma GCC optimize(2)
#include "MSA_worker.h"
using namespace std;

int main(int argc, char** argv)
{
    MSA_worker msa;
    msa.data_load(argv[1]);
    msa.make_worker();
    msa.query(argv[2]);
    return 0;
}