#include "CutFEMTestApp.h"
#include "MooseMain.h"

int
main(int argc, char * argv[])
{
  return Moose::main<CutFEMTestApp>(argc, argv);
}
