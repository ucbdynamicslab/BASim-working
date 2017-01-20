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
    AddOption("natural-curvature-switch", "natural curvature after switch", m_curvatureSwitch);
    AddOption("switch-time", "Time at which natural curvature changes", switchTime);
    
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
    
    l = GetScalarOpt("length");
    m_curvature = GetScalarOpt("natural-curvature");
    m_curvatureSwitch = GetScalarOpt("natural-curvature-switch");
    switchTime = GetScalarOpt("switch-time");
    
    cout << "Initial curvature: " << m_curvature << endl;
    cout << "Secondary curvature: " << m_curvatureSwitch << endl;
    cout << "Switch time: " << switchTime << endl;
        
    std::vector<Vec3d> vertices; // initial configuration
    std::vector<Vec3d> undeformed; // undeformed configuration
    
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

    // For dual curvature simulation ( mkjawed@andrew.cmu.edu )
    // Figure out the configuration after switch
    if (m_curvatureSwitch == 0.0) //straight rod
    {
        for (int i = 0; i < opts.numVertices; i++)
        {
            undeformedSwitch.push_back(Vec3d(l * i / double(opts.numVertices - 1), 0.0, 0.0));
        }
    }
    else if (m_curvatureSwitch > 0.0)
    {
        Scalar R0  = 1.0 / m_curvatureSwitch;
        Scalar x, y, z; /*Utility variables*/
        
        double phi;
        for (int i = 0; i < opts.numVertices; i++)
        {
            phi = ( ((double)i) * l /
                   ( (double) ( opts.numVertices - 1.0 ) ) ) /R0;
            x = R0 * sin(phi);
            y = R0 * ( 1.0 - cos(phi) );
            z = 0.0;
            undeformedSwitch.push_back(Vec3d(x, y, z));
        }
    }
    else
    {
        std::cout << "Negative switch-curvature not supported: exiting\n";
        exit(1);
    }
    rod->setswitchTime(switchTime);
    computeNaturalCurvatureSwitch();
    
    int nv = rod->nv();
    int ne = rod->ne();
    
    stepper = getRodTimeStepper(*rod);

    m_world->addObject(rod);
    m_world->addController(stepper);
}


void HangingRod::AtEachTimestep()
{
    // I am setting the time in the rod object so that bending force can access it
    rod->setcurrentTime(getTime());
    
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

void HangingRod::computeNaturalCurvatureSwitch()
{
    // Step 1: compute reference director
    // Copied from ElasticRod::computeSpaceParallel()
    // transport first edge in time
    Vec3d t0 = rod->getTangent(0);
    Vec3d tSwitch = undeformedSwitch[1] - undeformedSwitch[0];
    tSwitch = tSwitch / tSwitch.norm();
    Vec3d d10 = rod->getReferenceDirector1(0);
    Vec3d u = parallel_transport(d10, t0, tSwitch);
    u = (u - u.dot(tSwitch) * tSwitch).normalized();
    
    std::vector<Vec3d> d1Switch, d2Switch; // reference directors
    d1Switch.push_back(u);
    d2Switch.push_back(tSwitch.cross(u));
    
    t0 = tSwitch;
    // transport along centerline (Bishop frame)
    for (int i=1; i < rod->ne(); i++) {
        tSwitch = undeformedSwitch[i+1] - undeformedSwitch[i];
        tSwitch = tSwitch / tSwitch.norm();
        u = parallel_transport(u, t0, tSwitch);
        u = (u - u.dot(tSwitch) * tSwitch).normalized();
        d1Switch.push_back(u);
        d2Switch.push_back(tSwitch.cross(u));
        t0 = tSwitch;
    }
    
    // Step 2: compute material director
    // We could just set material directors equal to reference directors
    // Copied from ElasticRod::computeMaterialDirectors()
    std::vector<Vec3d> m1Switch, m2Switch; // material directors
    for (int j = 0; j < rod->ne(); ++j)
    {
        Scalar c = cos(rod->getTheta(j));
        Scalar s = sin(rod->getTheta(j));
        Vec3d u = d1Switch[j];
        Vec3d v = d2Switch[j];
        Vec3d Material1 =  c * u + s * v;
        Vec3d Material2 = -s * u + c * v;
        m1Switch.push_back(Material1);
        m2Switch.push_back(Material2);
    }
    
    //Step 3: compute curvature after switch
    for (int i = 1; i < rod->ne(); ++i)
    {
        Vec3d t1 = undeformedSwitch[i] - undeformedSwitch[i-1];
        Vec3d t2 = undeformedSwitch[i+1] - undeformedSwitch[i];
        t1 = t1 / t1.norm();
        t2 = t2 / t2.norm();
        Vec3d kb;
        Vec3d& kbT = kb;
        computeCurvatureBinormal(kb, t1, t2);
        Vec3d m1e = m1Switch[i - 1];
        Vec3d m2e = m2Switch[i - 1];
        Vec3d m1f = m1Switch[i];
        Vec3d m2f = m2Switch[i];
        Scalar curvSwitch1 = 0.5 * kb.dot(m2e + m2f);
        Scalar curvSwitch2 = - 0.5 * kb.dot(m1e + m1f);
        rod->setkappaBarSwitch1(i, curvSwitch1);
        rod->setkappaBarSwitch2(i, curvSwitch2);
    }
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
