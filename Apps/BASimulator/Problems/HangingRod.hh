#include "ProblemBase.hh"
#include "stdio.h"
#include <fstream>


class HangingRod : public Problem
{
public:
    HangingRod();
    virtual ~HangingRod();

protected:
  
    void Setup();
    void AtEachTimestep();
    void AfterEachTimestep();
    void Render();
    
    ElasticRod* rod;
    RodTimeStepper* stepper;
    
    Scalar m_curvature; // natural curvature
    Scalar rps; // rotation per second
    Scalar dtheta; // angle every dt time
    Scalar l; // length
};

