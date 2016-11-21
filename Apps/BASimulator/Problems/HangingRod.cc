#include "HangingRod.hh"
/* for use with input_shape.txt'*/
#include <iostream>

using namespace std;

HangingRod::HangingRod() : Problem("Hanging Rod in air", "A soft filament in air")
, rod(NULL)
, stepper(NULL)
{
    addDynamicsProps();
    addRodOptions();
    addRodTimeStepperOptions();
    
    AddOption("length", "length of rod", 0.1);
    AddOption("natural-curvature", "natural curvature", m_curvature);
    AddOption("rps", "rotation per second", rps);
    
    // default to no velocity-dependent forces
    GetScalarOpt("viscosity") = 0.0;
    GetScalarOpt("mass-damping") = 0.0;
}

HangingRod::~HangingRod()
{
    if (rod != NULL) delete rod;
    if (stepper != NULL) delete stepper;
}

void HangingRod::Setup()
{
    loadDynamicsProps();
    
    RodOptions opts;
    getRodOptions(opts);
    
    m_curvature = GetScalarOpt("natural-curvature");
    
    l = GetScalarOpt("length");
    rps = GetScalarOpt("rps");
    
    std::vector<Vec3d> vertices, undeformed;
    
    if (m_curvature == 0.0) //straight rod
    {
        for (int i = 0; i < opts.numVertices; i++)
        {
              vertices.push_back(Vec3d(l * i / double(opts.numVertices - 1), 0.0, 0.0));
            undeformed.push_back(Vec3d(l * i / double(opts.numVertices - 1), 0.0, 0.0));

        }
    }
    else if (m_curvature > 0.0)
    {
        Scalar R0  = 1.0 / m_curvature;
        Scalar x, y, z; /*Utility variables*/
        
        double phi;
        for (int i = 0; i < opts.numVertices; i++)
        {
            phi = ( ((double)i) * l /
                   ( (double) ( opts.numVertices - 1.0 ) ) ) /R0;
            x = R0 * sin(phi);
            y = R0 * ( 1.0 - cos(phi) );
            z = 0.0;
            
            vertices.push_back(Vec3d(x, y, z));
            undeformed.push_back(Vec3d(x, y, z));
        }
    }
    else
    {
        std::cout << "Negative curvature not supported: exiting\n";
        exit(1);
    }
    
    
    rod = setupRod(opts, vertices, undeformed);
    int nv = rod->nv();
    int ne = rod->ne();
    
    stepper = getRodTimeStepper(*rod);

    m_world->addObject(rod);
    m_world->addController(stepper);
}


void HangingRod::AtEachTimestep()
{
    // Clamping: set the position of all the nodes
    for (int i = 0; i < 2; i++)
    {
        Vec3d pos = rod->getVertex(i);
        stepper->getBoundaryCondition()->setDesiredVertexPosition(i, Vec3d(pos(0),
                                            pos(1), pos(2)));
    }
}

void HangingRod::AfterEachTimestep()
{
    ;
}


void HangingRod::Render()
{
    Scalar l = 0.2;
    
    Scalar scaleX = 0.328;  // conversion from dm to feet
    Scalar scaleYZ = 0.328; // 3.937; // conversion from dm to inches
    glScalef(scaleX, scaleYZ, scaleYZ);
    
    glDisable(GL_LIGHTING);
    
    glBegin(GL_LINES);
    glColor3f(1.0f, 0.0f, 0.0f);
    glVertex3f(-l, 0, 0);
    glVertex3f(l, 0, 0);     
    glEnd();

    if (true)
    {
        glDepthMask(GL_TRUE);
        //glDisable(GL_DEPTH_TEST);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        //glBegin(GL_QUADS);
        glBegin(GL_LINES);
        glEnd();
        glDisable(GL_BLEND);
    }

}
