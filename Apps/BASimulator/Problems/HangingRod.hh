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
    Scalar switchTime; // time at which switch happens
    Scalar l; // length
    
    Scalar m_curvatureSwitch; // natural curvature after switch
    std::vector<Vec3d> undeformedSwitch; // nodes of the natural state of the rod (after switch)
    void computeNaturalCurvatureSwitch(); // compute the natural curvature vectors at each node
};

