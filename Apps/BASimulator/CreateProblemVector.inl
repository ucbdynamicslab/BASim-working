#include "Problems/SewingMachine.hh"
#include "Problems/HangingRod.hh"

void CreateProblemVector()
{
  problems.push_back(NULL);
  problems.push_back(new SewingMachine());
  problems.push_back(new HangingRod());
}
