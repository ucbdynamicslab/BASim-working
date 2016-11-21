// rods/src/problems/SewingMachine.cc
// author: Miklos Bergou (miklos@cs.columbia.edu)
// December 19, 2008
//----------------------------------------------------------------------------
#include "SewingMachine.hh"
#include "BASim/src/Core/Util.hh"
//#include "BASim/src/Util/Utils.hh"
#include "JitterFunction.hh"

//#include "objects/SimpleAdaptivity.hh"
#include <cstdlib>
#include <ctime>
#include <limits>
#include <iomanip>
#include <queue>

#include <iostream>
#include <sstream>

#define REPORT_PHI              (false)
#define REPORT_COILING_RADIUS   (false)
#define REPORT_PD_DATA          (false)

// the max simulation time, used to determine when the program ends
extern Scalar max_time;

using namespace std;

//----------------------------------------------------------------------------
SewingMachine::SewingMachine()
  : Problem("SewingMachine", "SewingMachine")
  , m_rod(NULL)
  , m_br_stepper(NULL)  
  , m_rodLength(1.0)
  , m_fallHeight(30.0)
  , m_beltSpeed(0.0)
  , m_nozzleSpeed(3.0)
  , m_headVertices(3)
  , m_tailVertices(3)
  , m_numRemovedVertices(0)
  , m_jitter(2)
	, m_contact_method(0)
	, m_friction_method(0)
  , fft_n(0)
  , fft_in(0)
  , fft_out(0)
  , fft_plan()
  , contact_history()
  , m_step(0)
  , m_base_frequency(1)
{

  addDynamicsProps();
  addRodOptions();
  addRodTimeStepperOptions();

  GetScalarOpt("atol") = 1e-8;
 
  AddOption("rod-length", "rod length", m_rodLength);
	AddOption("fall-height", "fall height", m_fallHeight);
  AddOption("belt-speed", "belt speed", m_beltSpeed);
  AddOption("nozzle-speed", "nozzle speed", m_nozzleSpeed);
  AddOption("head-vertices", "head vertices", m_headVertices);
  AddOption("tail-vertices", "tail vertices", m_tailVertices);

  AddOption("jitter-amplitude", "amplitude of jitter", 0.0001);
  AddOption("jitter-frequency", "frequency of jitter", 1.0);
  
  AddOption("penalty-collision", "penalty-collision", true);
  AddOption("original-freeze", "original-freeze", false);
//  AddOption("friction-ground", "friction-ground", 0.3);
//  AddOption("friction-coulomb-coeff", "friction-coulomb-coeff", 0.3);


  AddOption("penalty-collision-stiff", "penalty-collision-stiff", 0.0);
  AddOption("penalty-collision-damp", "penalty-collision-damp", 0.0);
  AddOption("penalty-collision-extrathick", "penalty-collision-extrathick", 0.0);

  AddOption("penalty-collision-edge-edge", "penalty-collision-edge-edge", 0.0);

  AddOption("cutting-point", "cutting-point", 0.0);

  AddOption("fix-twist-after-hit", "fix-twist-after-hit", false);

  AddOption("contact-method", "contact-method", 0);
  AddOption("friction-method", "friction-method", 0);

  AddOption("coulomb-static", "coulomb-static", 0.36);
  AddOption("coulomb-kinetic", "coulomb-kinetic", 0.31);
  
  AddOption("natural-curvature", "natural-curvature", 0.0);
  AddOption("natural-rad-curvature", "natural-rad-curvature", 0.0);
  AddOption("use-curved-nozzle", "whether or not to use the curved nozzle implementation", false);
  AddOption("floor-bc-fix-twist", "whether or not to fix the material frame on the floor", false);
  
  AddOption("inversion-test", "inversion-test", false);
  AddOption("hopf-test", "hopf-test", false);
  AddOption("hopf-epsilon", "hopf-epsilon", 0.0);
  AddOption("hopf-end", "hopf-end", 0.0);
  AddOption("hopf-cut", "hopf-cut", false);
  
  AddOption("draw-recorded-rod", "draw-recorded-rod", false);
  AddOption("draw-recorded-rod-file", "draw-recorded-rod-file", "");

  AddOption("hollow", "Use a hollow rod with inner-radius", false);
  AddOption("inner-radius", "Radius of empty space at the center of hollow rod", 0.0);

  AddOption("hopf-test2", "hopf-test2", false);
  AddOption("hopf-epsilon-additional-time", "hopf-epsilon-additional-time", 5000.0);

  AddOption("hopf-critical-curv-test", "hopf-critical-curv-test", false);
//  AddOption("hopf-critical-curv-start", "hopf-critical-curv-start", 0.005);

  AddOption("hopf-critical-curv-jump", "hopf-critical-curv-jump", 0.01);
  AddOption("hopf-critical-nonzero-threshold", "hopf-critical-nonzero-threshold", 0.1);
  
  AddOption("hopf-offset", "hopf-offset", false);
  AddOption("hopf-offset-length", "hopf-offset-length", 300.0);
  
  AddOption("fft-window-size", "size of the sliding window (in frames) used in FFT of contact point coordinates", 2000);  // default comes from Brun 2012
  AddOption("fft-interval", "interval between two FFT computations", 500);// default comes from Brun 2012
  
  AddOption("stepping-variable", "the name of the variable being stepped in the multi-step simulation ('empty' or 'none' = no stepping)", "belt-speed");
  AddOption("stepping-increment-method", "method for incrementing the stepping variable ('add' or 'multiply')", "add");
  AddOption("stepping-increment", "the value added/multiplied to the stepping variable at every step", 0.0);
  AddOption("stepping-number-of-steps", "the number of steps of the multi-step simulation", 1);
  
  AddOption("floor-bc-fix-twist-alpha", "alpha for the fixing the twist on belt. the twist will be fixed after a distance of alpha * Lgb away from the latest contact point", 0.0);
  
  is_epsilon_switch = false;
  epsilon_switch_step = 0;
  epsilon_switch_time[0] = 300;
  epsilon_switch_time[1] = 1600;
  epsilon_switch_time[2] = 2000;
  epsilon_switch_time[3] = 1600;
  epsilon_switch_time[4] = 2000;
  
  hopf_msg_time = 0;
  time_tick = 0;
  jitter_done = false;
  
  m_inv_times.clear();
  m_cut_pos.clear();
  m_cut_time.clear();
    
  hopf_eps_switch_varying_curv = false;
  
  hopf_crt_curv = 0.0;
  hopf_max_zero_curv = 0.0;
  hopf_min_nonzero_curv = 0.0;
  hopf_varying_step = 0.0;
  hopf_varying_curv_jump = 0.0;
  
  cut_length = 0.0;
  cut_length_udf = 0.0;
      
  length_info.clear();
  length_info_size = 0;
  
  contactp_info.clear();
      
  multi_n = 0;
  
  srand( time(0) );
}
//----------------------------------------------------------------------------
SewingMachine::~SewingMachine()
{
  if (m_rod != 0) delete m_rod;
}

/*
# contact-method for rod - belt
# 0 - do nothing
# 1 - fixing y dimension
# 2 - applying a force
# 3 - inelastic impulse after dynamic (original)
contact-method 1

# friction-method
# 0 - no friction
# 1 - implicit damping friction, -cv
friction-method 1
*/

void SewingMachine::loadAuxOptions()
{
  m_contact_method = GetIntOpt("contact-method");
  m_friction_method = GetIntOpt("friction-method");
  
  m_rodLength = GetScalarOpt("rod-length");
	m_fallHeight = GetScalarOpt("fall-height");
  m_beltSpeed = GetScalarOpt("belt-speed");

	if (GetBoolOpt("hopf-test") && GetScalarOpt("hopf-epsilon")>0.0)
	{
		Scalar epsilon = GetScalarOpt("hopf-epsilon");
		m_nozzleSpeed =  m_beltSpeed / (1.0 - epsilon);
	
	} else {
	  m_nozzleSpeed = GetScalarOpt("nozzle-speed");
	}

  m_nozzleSpeedInit = m_nozzleSpeed;

  m_headVertices = GetIntOpt("head-vertices");
  m_tailVertices = GetIntOpt("tail-vertices");
  m_jitter.setAmplitude(GetScalarOpt("major-radius") * GetScalarOpt("jitter-amplitude"));
  m_jitter.setFrequency(GetScalarOpt("jitter-frequency"));

  GetIntOpt("nv") = max(GetIntOpt("nv"), m_headVertices + 1);

  is_epsilon_switch = GetBoolOpt("hopf-test2");
  if (is_epsilon_switch) {
    cout << "Injection speed will be set to the belt speed\n";
    
    // start with epsilon = 0
    m_nozzleSpeed = m_beltSpeed;
  }
  
  hopf_crt_curv = GetScalarOpt("natural-curvature");
  if (GetScalarOpt("natural-rad-curvature") > 0.0) {
    hopf_crt_curv = 1.0 / GetScalarOpt("natural-rad-curvature");
  }
  
  hopf_eps_switch_varying_curv = GetBoolOpt("hopf-critical-curv-test");
  if (hopf_eps_switch_varying_curv) {
    hopf_varying_curv_jump = GetScalarOpt("hopf-critical-curv-jump");
  }
  
  
  youngM = GetScalarOpt("youngs-modulus");
  shearM = youngM / (2.0 * (1.0 + 0.5));
  
  cout << "youngM : " << youngM << "\n";
  cout << "shearM : " << shearM << "\n";
  
  
  hopf_offset_test = GetBoolOpt("hopf-offset");
  if (hopf_offset_test) {
    cout << "hopf_offset_test " << hopf_crt_curv << "\n";
    m_nozzleSpeedInit = m_nozzleSpeed = m_beltSpeed;
    hopf_offset_length = GetScalarOpt("hopf-offset-length");
    hopf_offset_endtime = (m_fallHeight + hopf_offset_length) / m_nozzleSpeed;
    cout << "hopf_offset_endtime " << hopf_offset_endtime << "\n";
    
    epsilon_switch_time[0] = hopf_offset_endtime + 1000.0;
    
  } else {
    hopf_offset_length = 0.0;
    hopf_offset_endtime = 0.0;
  }
  
  

}

void SewingMachine::setupTimeStepper()
{
  RodTimeStepper::Method integrator;
  if( GetStringOpt("integrator") == "implicit" ) integrator = RodTimeStepper::IMPL_EULER;
  else assert( false );
  
  m_br_stepper = new BridsonStepper( RodTimeStepper::IMPL_EULER, GetIntOpt("iterations"), GetScalarOpt("dt"), GetScalarOpt("mass-damping"), GetVecOpt("gravity") );

//  m_br_stepper->setGravity(GetVecOpt("gravity"));

	m_br_stepper->setContactAndFrictionMethod( m_contact_method, m_friction_method );
	if (m_contact_method == 1 || m_contact_method == 2 ) {
		m_br_stepper->setGroundFix( true );
	} else {
		m_br_stepper->setGroundFix( false );
	}
		
  // damping friction
	m_br_stepper->setFriction( true ); //GetScalarOpt("friction-ground"), GetScalarOpt("friction-coulomb-coeff"));
	
	// Coulomb's friction
	m_br_stepper->setCoulombFriction(GetScalarOpt("coulomb-static"), GetScalarOpt("coulomb-kinetic"));

	m_br_stepper->setTwistFixAfterHit( GetBoolOpt("fix-twist-after-hit") );
  
  if (m_contact_method == 3) {
		m_br_stepper->disablePenaltyImpulses();
	  m_br_stepper->enableIterativeInelasticImpulses();
//  m_br_stepper->setEdgeEdgePenalty(GetScalarOpt("penalty-collision-edge-edge"));
//  m_br_stepper->setVertexFacePenalty(00.0); 
	  m_br_stepper->setEdgeEdgePenalty(0.0); 
  }
  
  if (m_contact_method == 4) {
		m_br_stepper->disablePenaltyImpulses();
	  m_br_stepper->disableIterativeInelasticImpulses();
//  m_br_stepper->setEdgeEdgePenalty(GetScalarOpt("penalty-collision-edge-edge"));
//  m_br_stepper->setVertexFacePenalty(00.0); 
	  m_br_stepper->setEdgeEdgePenalty(0.0); 
  }

}

void SewingMachine::makeRod()
{
  RodOptions opts;
  getRodOptions(opts);

  opts.YoungsModulus = youngM;
  opts.ShearModulus = shearM;
  
  opts.viscosity = 0.0;
  opts.viscosity = 1e-10;
  
  // setting for hollow rod
  opts.hollow = GetBoolOpt("hollow");
  opts.radiusInner = GetScalarOpt("inner-radius");
  
  if (opts.hollow) {
    cout << "hollow rod : " << opts.radiusInner << " cm\n";
  }

  Scalar radius = GetScalarOpt("major-radius");
  Scalar rho = GetScalarOpt("density");
  Scalar edgeLen = m_rodLength; 

  std::vector<Scalar> vals(2);
  m_jitter.getValues(vals, getTime());
  std::vector<Vec3d> vertices, undeformed;
  
  for (int i = 0; i < opts.numVertices; ++i) {
    vertices.push_back(Vec3d(vals[0], m_fallHeight + (Scalar)i * edgeLen, vals[1]));
	}

  m_rod = setupRod(opts, vertices, vertices);
  
  int nv = m_rod->nv();
  int ne = m_rod->ne();
  
  std::cout << "Injection Height : " << m_fallHeight << "\n";
  std::cout << "Injection Speed : " << m_nozzleSpeed << "\n";
  std::cout << "Belt Speed : " << m_beltSpeed << "\n";

  for (int i = 0; i < nv; ++i) {
  	m_rod->setVelocity(i, Vec3d(0, -m_nozzleSpeed, 0));
  }
 	 
  //////////////////
  // set undeformed curvature
//  if (GetScalarOpt("natural-curvature") > 0.0) {
//    std::cout << "Natural-curvature   " << GetScalarOpt("natural-curvature") << "\n";
  if (hopf_crt_curv > 0.0) {
    std::cout << "Natural-curvature   " << hopf_crt_curv << "\n";
        
    ElasticRod *rod = m_rod;
    std::vector<RodForce*>& forces = rod->getForces();
    std::vector<RodForce*>::iterator fIt;
    for (fIt = forces.begin(); fIt != forces.end(); ++fIt) {
      if (dynamic_cast<RodBendingForceSym*> (*fIt) != NULL ) {
        RodBendingForceSym *f = dynamic_cast<RodBendingForceSym*> (*fIt);
          
        if (!f->viscous()) {
          ElasticRod::vertex_iter vit = rod->vertices_begin();
          int i=0;
          ++vit;
      
          for (; vit != rod->vertices_end() && i < rod->nv() - 1; ++vit, ++i) {
            ElasticRod::vertex_handle vh = *vit;
            Vec2d kap = f->getKappaBar(vh);
//            kap(0) = hopf_crt_curv;
//            kap(1) = hopf_crt_curv;
            
//            kap(0) = hopf_crt_curv / sqrt(2.0);
//            kap(1) = hopf_crt_curv / sqrt(2.0);
            
//            kap(0) = -hopf_crt_curv;
//            kap(1) = 0.0;

//            kap(0) = 0.0;
//            kap(1) = hopf_crt_curv;

            kap(0) = hopf_crt_curv;
            kap(1) = 0.0;

            f->setKappaBar(vh, kap);
          }
        }
      }
    }
  }
  
  m_br_stepper->addRod(m_rod);
  m_world->addObject(m_rod);

  RodBoundaryCondition* boundary = m_rod->getBoundaryCondition();
  
  topVert = -1;

  for (int i = 1; i < nv; ++i) {
    if (m_rod->getVertex(i)[1] > m_fallHeight) {
      if (topVert < 0) topVert = i;
    
      boundary->setDesiredVertexPosition(i, m_rod->getVertex(i));
    
      if (boundary->isVertexScripted(i-1)) {
        boundary->setDesiredEdgeAngle(i-1, m_rod->getTheta(i-1));
      }
    }
  }

}

void SewingMachine::makeBelt() {
  tri_mesh = new TriangleMesh();
	Vec3d v = Vec3d(m_beltSpeed, 0, 0);
	
	{	  
	  TriangleMesh::vertex_handle vhx = tri_mesh->addVertex();
	  tri_mesh->getVertex(vhx) = Vec3d(-1000.0, 0.0, -1000.0);	
	  tri_mesh->getVelocity(vhx) = v;
	  
	  TriangleMesh::vertex_handle vhy = tri_mesh->addVertex();
	  tri_mesh->getVertex(vhy) = Vec3d(-1000.0, 0.0, +1000.0);
	  tri_mesh->getVelocity(vhy) = v;

	  TriangleMesh::vertex_handle vhz = tri_mesh->addVertex();
	  tri_mesh->getVertex(vhz) = Vec3d(+1000.0, 0.0, 0.0);
	  tri_mesh->getVelocity(vhz) = v;

	  tri_mesh->addEdge(vhx,vhy);
	  tri_mesh->addEdge(vhy,vhz);
	  tri_mesh->addEdge(vhz,vhx);
	  
	  // Generate a triangular face
	  tri_mesh->addFace(vhx,vhy,vhz);
	}
	
	tri_mesh->fixed = true;
	
	m_world->addObject(tri_mesh);
	m_br_stepper->addTriangleMesh(tri_mesh);
	
}

//----------------------------------------------------------------------------
void SewingMachine::Setup()
{
	std::cout << "SewingMachine : Start setup\n";
	
  loadDynamicsProps();

	loadAuxOptions();
	
	setupTimeStepper();
	
//	if (GetBoolOpt("draw-recorded-rod")) {
//		makeRodFromFile();
//		makeBelt();
//		m_world->addController(m_br_stepper);
//	} else {
  {	
		makeRod();
	
		makeBelt();

		m_br_stepper->prepareForExecution();
		m_world->addController(m_br_stepper);
	}  
  
  // initialize fft
  int N = GetIntOpt("fft-window-size");
  assert(N % 2 == 0);
  fft_n = N;
  fft_in =  (double *)      fftw_malloc(sizeof (double) * N);
  fft_out = (fftw_complex *)fftw_malloc(sizeof (fftw_complex) * (N / 2 + 1));
  fft_plan = fftw_plan_dft_r2c_1d(N, fft_in, fft_out, FFTW_MEASURE);
  
}

//----------------------------------------------------------------------------
void SewingMachine::updateFloorBoundary()
{
  int nv = m_rod->nv();
  int i = nv - 1;
  
  if (GetBoolOpt("original-freeze")) {
		RodBoundaryCondition* boundary = m_rod->getBoundaryCondition();
		
		while (i >= 0 && boundary->isVertexScripted(i)) --i;
		while (i >= 0 && !boundary->isVertexScripted(i)) --i;
		while (i >= 0) {
		  if (boundary->isVertexScripted(i)) {
		    m_rod->setVelocity(i, Vec3d(m_beltSpeed, 0, 0));
		  }
		  --i;
		}
	} else {
		// for method 7
		RodBoundaryCondition* boundary = m_rod->getBoundaryCondition();
		
		i = nv - 1;
		
		while (i >= 0 && boundary->isVertexScripted(i)) --i;
		while (i >= 0) {
		  if (boundary->isVertexScripted(i)) {
		    m_rod->setVelocity(i, Vec3d(m_beltSpeed, 0, 0));
		  }		
		
			i--;
		}
	
	}
  
}
//----------------------------------------------------------------------------
void SewingMachine::updateContainerBoundary()
{
  int nv = m_rod->nv();
  int ne = m_rod->ne();

  RodBoundaryCondition* boundary = m_rod->getBoundaryCondition(); 

  // release particles/edges that have fallen out of the container
  for (int i = nv - 1; i >= 0; --i) {
    if (!boundary->isVertexScripted(i)) break;
    if (m_rod->getVertex(i)[1] <= m_fallHeight) {

    	if (topVert >= i) {
    		topVert = i + 1;
    	}
    	
      boundary->releaseVertex(i);
      if (i < ne && boundary->isEdgeScripted(i)) {
        boundary->releaseEdge(i);
      }
    }
  }
	
  Scalar curvature = GetScalarOpt("natural-curvature");
  Scalar curve_radius = 1 / curvature;
  Vec3d curve_center;
  curve_center.x() = 0 + curve_radius;
  curve_center.y() = m_fallHeight + m_rodLength / 2;
  curve_center.z() = 0;
  
	bool is_vertex_added = false;
	
  // create new particles if too few remain in container
  std::vector<Scalar> vals(2);
  m_jitter.getValues(vals, getTime());
  Scalar edgeLen = m_rodLength; //min(m_fallHeight.getValue() / 20.0, m_rodLength);
  while (!boundary->isVertexScripted(nv - m_headVertices)
         || m_rod->getVertex(nv - 1)[1] <= m_fallHeight)
  {
    Vec3d v, vel;
    if (curvature == 0 || !GetBoolOpt("use-curved-nozzle"))
    {
      v = Vec3d(vals[0], m_rod->getVertex(nv - 1)[1] + edgeLen, vals[1]);
      vel = Vec3d(0, -m_nozzleSpeed, 0);
    } else
    {
      Scalar theta = ((m_rod->nv() - topVert) * 2 + 1) * asin(edgeLen / 2 / curve_radius);
      v = curve_center + curve_radius * Vec3d(-cos(theta), sin(theta), 0);

      Scalar theta2 = theta - 2 * asin(m_nozzleSpeed * getDt() / 2 / curve_radius);
      Vec3d v2 = curve_center + curve_radius * Vec3d(-cos(theta2), sin(theta2), 0);
//      vel = (v2 - v) / getDt(); // should be equivalent to the line below
      vel = (v2 - v).normalized() * m_nozzleSpeed;
    }
		
		m_rod->addVertexAtEnd(v, vel);
	  m_br_stepper->resizeRodVertex();

    is_vertex_added = true;
    nv = m_rod->nv();
    
    boundary->setDesiredVertexPosition(nv-1, m_rod->getVertex(nv-1));
    if (boundary->isVertexScripted(nv-2)) {
    	boundary->setDesiredEdgeAngle(m_rod->ne()-1, m_rod->getTheta(m_rod->ne()-1));
    }
  }
  
  // modify the vertex scripting within the container to create a curved shape due to rod's natural curvature
  if (topVert > 0 && curvature != 0 && GetBoolOpt("use-curved-nozzle"))
  {
    Scalar theta0 = asin((m_rod->getVertex(topVert)[1] - curve_center[1]) / curve_radius);
    for (int i = topVert; i < m_rod->nv(); i++)
    {
      Scalar theta = theta0 + 2 * (i - topVert) * asin(edgeLen / 2 / curve_radius);
      Vec3d v = curve_center + curve_radius * Vec3d(-cos(theta), sin(theta), 0);
      
      Scalar theta2 = theta - 2 * asin(m_nozzleSpeed * getDt() / 2 / curve_radius);
      Vec3d v2 = curve_center + curve_radius * Vec3d(-cos(theta2), sin(theta2), 0);
//      Vec3d vel = (v2 - v) / getDt(); // should be equivalent to the line below
      Vec3d vel = (v2 - v).normalized() * m_nozzleSpeed;
      
      m_rod->setVertex(i, v);
      boundary->setDesiredVertexPosition(i, v);
      m_rod->setVelocity(i, vel);
    }
    
    // alter the bending stiffness of the first (lowest) vertex in the nozzle (all the other vertices have regular stiffness)
    for (int i = 0; i < m_rod->nv(); i++)
      m_rod->setBendingStiffnessMultiplier(i, 1.0);

    Scalar multiplier = pow(1000.0, (theta0 / asin(edgeLen / 2 / curve_radius) + 1) / 2);
    multiplier = 1;
    m_rod->setBendingStiffnessMultiplier(topVert, multiplier); // topVert is the first vertex in the nozzle
    
    for (int i = 0; i < m_rod->nv(); i++)
      m_rod->setViscousBendingStiffnessOverride(i, 0.0);
    
    Scalar viscosity = 1e+8;
    m_rod->setViscousBendingStiffnessOverride(topVert, viscosity);
    
  }
  
  if (hopf_crt_curv > 0.0) {
    ElasticRod *rod = m_rod;
    std::vector<RodForce*>& forces = rod->getForces();
    std::vector<RodForce*>::iterator fIt;
    for (fIt = forces.begin(); fIt != forces.end(); ++fIt) {
      if (dynamic_cast<RodBendingForceSym*> (*fIt) != NULL ) {
        RodBendingForceSym *f = dynamic_cast<RodBendingForceSym*> (*fIt);
          
        if (!f->viscous()) {
          ElasticRod::vertex_iter vit = rod->vertices_begin();
          int i=0;
            
          ++vit;
      
          for (; vit != rod->vertices_end() && i < rod->nv() - 1; ++vit, ++i) {
            ElasticRod::vertex_handle vh = *vit;
            Vec2d kap = f->getKappaBar(vh);
//            kap(0) = hopf_crt_curv;
//            kap(1) = hopf_crt_curv;

//            kap(0) = hopf_crt_curv / sqrt(2.0);
//            kap(1) = hopf_crt_curv / sqrt(2.0);
            
//            kap(0) = -hopf_crt_curv;
//            kap(1) = 0.0;

//            kap(0) = 0.0;
//            kap(1) = hopf_crt_curv;

            kap(0) = hopf_crt_curv;
            kap(1) = 0.0;
            
            f->setKappaBar(vh, kap);
            //std::cout << f->getKappaBar(vh) << "\n";
          }
        }
      }
    }
    m_rod->updateProperties();
  }
}

//----------------------------------------------------------------------------
void SewingMachine::freezeCollidedVertices(bool saveToFile)
{
  int nv = m_rod->nv();
  RodBoundaryCondition* boundary = m_rod->getBoundaryCondition();
	
  if (!GetBoolOpt("original-freeze")) return;
  
  for (int i = 0; i < nv; ++i) {
    const Vec3d& v = m_rod->getVertex(i);
    Scalar rad = m_rod->radiusA(min(i, nv-2));
    if (v[1] <= rad) {
      if (!boundary->isVertexScripted(i)) {
        boundary->setDesiredVertexPosition(i, m_rod->getVertex(i));
        m_rod->setVelocity(i, Vec3d(m_beltSpeed, 0, 0));
        
        if (!GetBoolOpt("fix-twist-after-hit")) continue;

        if (i > 0 && boundary->isVertexScripted(i-1) && !boundary->isEdgeScripted(i-1)) {
        	boundary->setDesiredEdgeAngle(i-1, m_rod->getTheta(i-1));
          m_rod->setThetaDot(i-1, 0);
        }
        if (i < nv - 1 && boundary->isVertexScripted(i+1) && !boundary->isEdgeScripted(i)) {
					boundary->setDesiredEdgeAngle(i, m_rod->getTheta(i));
          m_rod->setThetaDot(i, 0);
        }
      }
    }
  }
}
//----------------------------------------------------------------------------
void SewingMachine::removeVertex()
{
  ++m_numRemovedVertices;

	m_rod->removeVertexFromFront();
  m_br_stepper->resizeRodVertex();
	
}
      
void SewingMachine::cutRod()
{
  if (GetBoolOpt("hopf-cut") && GetBoolOpt("hopf-test"))
  {
    const int k = 150;
    RodBoundaryCondition* boundary = m_rod->getBoundaryCondition();    
    
    int nfix = 0;
    for(int i=0; i<m_rod->nv() && i<k; i++) {
      if (boundary->isVertexScripted(i)) {    
        nfix++;
      }
    }
    
    if (nfix > k/2) {
        m_cut_pos.push_back(m_rod->getVertex(0));
        m_cut_time.push_back(getTime());
        
        if (epsilon_switch_step == 4) {
          Vec3d e = m_rod->getVertex(0) - m_rod->getVertex(1);
          cut_length += e.norm();
          cut_length_udf += m_rodLength;
        }
        removeVertex();
    }

/*    
    for(int i=0; i<m_rod->nv(); i++) {
      if (boundary->isVertexScripted(i)) {
        if (i == k) {
          m_cut_pos.push_back(m_rod->getVertex(0));
          m_cut_time.push_back(getTime());
          removeVertex();
          break;
        }
      } else {
        break;
      }
    }  
*/  
    return;
  }
  
	Scalar cutDist = GetScalarOpt("cutting-point");
	
	if (cutDist > 0.0) 
	{
		int cutId = -1;
		
//		for(int i=0; i<m_rod->nv(); i++) {
//			if (fabs(m_rod->getVertex(i)[0]) > cutDist) {
//				cutId = i; 
//			}
//		}
    
    cutId = (int)m_rod->nv() - (int)(cutDist / m_rodLength);
    if (cutId < 0)
      cutId = -1;
		
		if (cutId == -1) return;

		while(cutId > -1) {	// delete first vertex and edge
			removeVertex();
		
			cutId--;
		}		
	}

}

// Returns the id of vertex which was fixed last.
// Most operations only use vertices fixed on the belt, ignoring vertices still moving.
// If a vertex has an ID smaller than this return value, it has been already fixed.
// If -1, nothing is fixed except vertices around injector.
int SewingMachine::GetVertexIdFixedLast()
{
  int last_vertex_fixed = -1;
  for(int i=m_rod->nv() - 5; i>= 0; i--) {
    if (m_rod->getBoundaryCondition()->isVertexScripted(i)) {          
      last_vertex_fixed = i; break;
    }
  }	

  return last_vertex_fixed;
}

// [n_points] : number of extreme points that we want to get. ex) 10 means getting the last 10 extreme points
// if n_points == 0, extract all extreme points.
// This returns the number of vertices actually extracted.
int SewingMachine::GetExtremePoints(vector<Vec2d> &ext_points, const int n_points)
{
  ext_points.clear();
  
  int last_vertex_fixed = GetVertexIdFixedLast();
  if (last_vertex_fixed < 3) return 0;

  vector<double> x_vals;
  vector<double> z_vals;
    
  for(int i=last_vertex_fixed; i >= 0; i--) {
    x_vals.push_back(m_rod->getVertex(i)(0));
    z_vals.push_back(m_rod->getVertex(i)(2));
  }
  
	for(int i=(int)m_cut_pos.size()-1; i>=0; i--) {
    x_vals.push_back(m_cut_pos[i](0));
    z_vals.push_back(m_cut_pos[i](2));
	}
	
	int n = 0;
	for(int i=1; i<z_vals.size()-1; i++) {
		if ((z_vals[i+1]>z_vals[i] && z_vals[i-1]>z_vals[i]) || (z_vals[i+1]<z_vals[i] && z_vals[i-1]<z_vals[i])) {
			ext_points.push_back( Vec2d(x_vals[i], z_vals[i]) );
      n++;		
      if (n == n_points) break;
		}	
  }
  
  return n;
}


// This computes amplitudes and periods of recent [n_peaks] peaks of oscillations.
// if [n_peaks] == 0, compute all of them
// This returns the number of oscillations actually computed.
// [x_positions] x in (x,y,z)
int SewingMachine::ComputeAmplitudeAndWaveLength(vector<double> &x_positions, double &mean_amp, vector<double> &amplitudes, double &mean_wavelength, vector <double>& wavelengths, const int n_peaks)
{
  mean_amp = 0;
  mean_wavelength = 0;
  amplitudes.clear();
  wavelengths.clear();
  x_positions.clear();
  
  vector<Vec2d> ext_points;
  
  int npts;
  if (n_peaks == 0) {
    npts = GetExtremePoints(ext_points, 0);
    if (npts < 2) return 0;
  } else {
    npts = GetExtremePoints(ext_points, n_peaks + 1);
    if (npts < n_peaks+1) return 0;
  }
  

  int n = 0;
  double amp_sum = 0.0;
  double wvlen_sum = 0.0;
  for(int i=0; i<npts-1; i++) {
    double amp = abs(ext_points[i+1](1) - ext_points[i](1));
    double wvlen = abs(ext_points[i+1](0) - ext_points[i](0)) * 2.0;
    
    amplitudes.push_back(amp);
    wavelengths.push_back(wvlen);
    x_positions.push_back((ext_points[i+1](0) + ext_points[i](0)) * 0.5);
    
    amp_sum += amp;
    wvlen_sum += wvlen;
    n++;
  }

  mean_amp = amp_sum / (double)n;
  mean_wavelength = wvlen_sum / (double)n;
  
  return n;
}

std::vector<SewingMachine::FFTPeak> SewingMachine::fftFindPeaks(fftw_complex * out)
{
  std::vector<FFTPeak> peaks;
  for (int i = 0; i < fft_n / 2 + 1; i++)
  {
    Scalar i1 = (i > 0 ? sqrt(square(fft_out[i - 1][0]) + square(fft_out[i - 1][1])) : 0);
    Scalar i2 = sqrt(square(fft_out[i][0]) + square(fft_out[i][1]));
    Scalar i3 = (i < fft_n / 2 ? sqrt(square(fft_out[i + 1][0]) + square(fft_out[i + 1][1])) : 0);
    
    if (i2 > i1 && i2 > i3)
    {
      FFTPeak p;
      p.intensity = i2;

      Scalar omega = i;
      Scalar phi = atan2(fft_out[i][1], fft_out[i][0]);
      
      // estimate peak position with the peak sample and two neighboring sample
      Scalar i_nb = (i3 > i1 ? i3 : i1);
      Scalar omega_nb = (i3 > i1 ? (i + 1) : (i - 1));
      Scalar phi_nb = (i3 > i1 ? atan2(fft_out[i + 1][1], fft_out[i + 1][0]) : atan2(fft_out[i - 1][1], fft_out[i - 1][0]));  // the side with stronger signal cannot be the side with boundary
      
      p.omega = omega + (omega_nb - omega) * i_nb / (i2 + i_nb);
      p.phi = phi + (phi_nb - phi) * i_nb / (i2 + i_nb);
      
      peaks.push_back(p);
    } 
  }
  
//  std::cout << "peak# = " << peaks.size() << std::endl;
  
//  std::sort(peaks.begin(), peaks.end());
  
  Scalar maxpeak = 0;
  for (size_t i = 0; i < peaks.size(); i++)
    if (peaks[i].intensity > maxpeak)
      maxpeak = peaks[i].intensity;
  
  std::vector<FFTPeak> significant_peaks;
  for (size_t i = 0; i < peaks.size(); i++)
    if (peaks[i].intensity / maxpeak > .15)
      significant_peaks.push_back(peaks[i]);
  
  std::sort(significant_peaks.begin(), significant_peaks.end());
  
  return significant_peaks;
}

bool SewingMachine::ApproxMatch(Scalar omega, Scalar omegabase, int multiplier, Scalar tolerance)
{
  return fabs(omega / omegabase - multiplier) < tolerance;
}

bool SewingMachine::steadyStateReached()
{
  int nsample = (int)(1.5 * fft_n / m_base_frequency);
  int ncycle = 10;    // required number of stable cycles (< 5% variation in amplitude)
  int npattern = 10;  // required number of stable pattern (same patterns detected for the last npattern ffts)
  
  // reject steady state if the detected patterns haven't been staying constant, or not enough of them exist, or we're just seeing unknown patterns
  if ((int)pattern_detected.size() < npattern)
    return false;
  
  if (pattern_detected.back() == 0)
    return false;
      
  int last_pattern = pattern_detected.back();
  for (int i = 1; i < npattern; i++)
    if (last_pattern != pattern_detected[pattern_detected.size() - i - 1])
      return false;
  
  // not enough samples in contact_history for steady state determination
  if ((int)contact_history.size() < nsample + (int)((ncycle - 1) * fft_n / m_base_frequency))
    return false;
  
  // scan ncycle periods from the end of contact_history, finding the amplitude of each cycle
  Scalar ymin = std::numeric_limits<Scalar>::infinity();
  Scalar ymax = -std::numeric_limits<Scalar>::infinity();
  std::vector<Scalar> amplitudes(nsample);
  for (int i = 0; i < ncycle; i++)
  {
    int start = (int)(i * fft_n / m_base_frequency);
    for (int j = 0; j < nsample; j++)
    {
      int sampleid = contact_history.size() - start - j - 1;
      Vec2d & sample = contact_history[sampleid];
      if (sample.y() > ymax)
        ymax = sample.y();
      if (sample.y() < ymin)
        ymin = sample.y();
    }
    
    Scalar amplitude = ymax - ymin;
    assert(amplitude > 0);
    amplitudes[i] = amplitude;
  }
  
  // check the deviation of each individual amplitude against mean
  Scalar amp_mean = 0;
  for (int i = 0; i < ncycle; i++)
    amp_mean += amplitudes[i];
  amp_mean /= ncycle;
  
  for (int i = 0; i < ncycle; i++)
    if (fabs(amplitudes[i] - amp_mean) / amp_mean > 0.05)
      return false;

  // if all amplitudes are within range of +/- 5% of mean, and all recently detected patterns match. then steady state is reached
  return true;
}

void SewingMachine::updateMultiStep()
{
  bool steady_state = steadyStateReached();
  Scalar time = getTime();
  Scalar step_time = time - m_step_start_time;
  
  if ((time * m_nozzleSpeed > 500 && steady_state) || step_time * m_nozzleSpeed > 10000) // total > 500cm of rod + steady state, or time out after 10000cm has been deposited for this step (at injection speed 1.0cm/s, 10000cm is roughly two orders of magnitude higher than the wavelength of some of the long patterns
  {
    if (steady_state && REPORT_PD_DATA)
      std::cerr << "Steady state reached." << std::endl;
    else
      std::cerr << "Current step time out." << std::endl;
    
    // time up or steady state is reached. next step.
    m_step++;
    m_step_start_time = time;
    if (m_step >= GetIntOpt("stepping-number-of-steps"))
    {
      max_time = time; // this will terminate the program
      return;
    }
    
    std::string variablename = GetStringOpt("stepping-variable");
    Scalar * variable;
    if (variablename == "belt-speed")
    {
      variable = &m_beltSpeed;
    } else if (variablename == "nozzle-speed")
    {
      variable = &m_nozzleSpeed;
    } else if (variablename == "fall-height")
    {
      variable = &m_fallHeight;
    } else if (variablename == "natural-curvature")
    {
      variable = NULL;
    }
    
    if (variable)
    {
      if (GetStringOpt("stepping-increment-method") == "add")
        (*variable) += GetScalarOpt("stepping-increment");
      else
        (*variable) *= GetScalarOpt("stepping-increment");
    }
    
    std::cerr << "Next step at T = " << time << ": " << variablename << " = " << *variable << std::endl;
    std::cout << "Next step at T = " << time << ": " << variablename << " = " << *variable << std::endl;
    
    contact_history.clear();
    pattern_detected.clear();
    
    m_base_frequency = 1;
  }
  
}

//----------------------------------------------------------------------------
void SewingMachine::AtEachTimestep()
{
    // print the current twist at every vertex, and sum them up to get Phi
    if (REPORT_PHI)
    {
        Scalar Phi = 0;
        static int s_vifl = -1;
        bool new_fixed = false;
        if (GetVertexIdFixedLast() > 0)
        {
            if (s_vifl != GetVertexIdFixedLast())
            {
                s_vifl = GetVertexIdFixedLast();
                new_fixed = true;
            }

            for (int i = GetVertexIdFixedLast(); i < m_rod->nv() - 1; i++)
            {
                Scalar twist = m_rod->getReferenceTwist(i) + m_rod->getTheta(i) - m_rod->getTheta(i - 1);
                Phi += twist;
            }
//            std::cout << "Phi = " << Phi << std::endl;
       
            // the rotation is fixed after a segment of free rod of length alpha * Lgb
            Scalar alpha = GetScalarOpt("floor-bc-fix-twist-alpha");
            Scalar r = GetScalarOpt("major-radius");
            Scalar E = GetScalarOpt("youngs-modulus");
            Scalar rho = GetScalarOpt("density");
            Scalar g = GetVecOpt("gravity").norm();
            Scalar Lgb = pow(r * r * E / (8 * rho * g), 1.0 / 3.0);
            int free_edges = Lgb * alpha / GetScalarOpt("rod-length");
            int hook_edges = Lgb * M_PI / GetScalarOpt("rod-length");  // this is the length of the hook, between the contact point and the relatively straight and homogeneous "heel" up to the nozzle
            
            Scalar fixed_max_twist = 0;
            Scalar fixed_total_twist = 0;
            for (int i = 1; i < GetVertexIdFixedLast() - free_edges; i++)
            {
                Scalar twist = m_rod->getReferenceTwist(i) + m_rod->getTheta(i) - m_rod->getTheta(i - 1);
                fixed_total_twist += twist;
                fixed_max_twist = std::max(fixed_max_twist, fabs(twist));
            }
           
            Scalar ground_max_twist = 0;
            Scalar ground_total_twist = 0;
            for (int i = 1; i < GetVertexIdFixedLast(); i++)
            {
                Scalar twist = m_rod->getReferenceTwist(i) + m_rod->getTheta(i) - m_rod->getTheta(i - 1);
                ground_total_twist += twist;
                ground_max_twist = std::max(ground_max_twist, fabs(twist));
            }
            
            Scalar free_max_twist = 0;
            Scalar free_min_twist = 0;
            Scalar free_total_twist = 0;
            Scalar free_twist_ss = 0;
            for (int i = GetVertexIdFixedLast() + hook_edges; i < m_rod->nv() - 2; i++)
            {
                Scalar twist = m_rod->getReferenceTwist(i) + m_rod->getTheta(i) - m_rod->getTheta(i - 1);
                free_total_twist += twist;
                free_twist_ss += twist * twist;
                
                if (i == GetVertexIdFixedLast() + hook_edges)
                {
                    free_max_twist = twist;
                    free_min_twist = twist;
                } else
                {
                    free_max_twist = std::max(free_max_twist, twist);
                    free_min_twist = std::min(free_min_twist, twist);
                }
            }
            int N = (m_rod->nv() - 2 - GetVertexIdFixedLast() - hook_edges);
            Scalar free_mean_twist = free_total_twist / N;
            Scalar free_twist_std = sqrt((free_twist_ss * N - free_total_twist * free_total_twist) / (N * (N - 1)));
            
            static Scalar free_twist_alltime_max = 0;
            if (fabs(free_min_twist) > free_twist_alltime_max)
                free_twist_alltime_max = fabs(free_min_twist);
            if (fabs(free_max_twist) > free_twist_alltime_max)
                free_twist_alltime_max = fabs(free_max_twist);
            
            // turning direction in contact point
            Scalar turning = 0;
            if (GetVertexIdFixedLast() >= 2)
            {
                Vec3d contact0 = m_rod->getVertex(GetVertexIdFixedLast());
                Vec3d contact1 = m_rod->getVertex(GetVertexIdFixedLast() - 1);
                Vec3d contact2 = m_rod->getVertex(GetVertexIdFixedLast() - 2);
                turning = (contact0 - contact1).cross(contact1 - contact2).y();
            }
            
//            std::cout << fixed_max_twist << " " << fixed_total_twist << " " <<    ground_max_twist << " " << ground_total_twist << " " << free_max_twist << " " << free_total_twist << std::endl;
            
            Scalar recently_fixed_twist = 0;
            int recently_fixed_vertex = GetVertexIdFixedLast() - free_edges - 3;
            if (recently_fixed_vertex > 0)
            {
                recently_fixed_twist = m_rod->getReferenceTwist(recently_fixed_vertex) + m_rod->getTheta(recently_fixed_vertex) - m_rod->getTheta(recently_fixed_vertex - 1);
            }
            
            static std::deque<Scalar> recently_fixed_twists;
            if (new_fixed)
                recently_fixed_twists.push_back(recently_fixed_twist);
            while (recently_fixed_twists.size() > 214)
                recently_fixed_twists.pop_front();
            Scalar recently_fixed_twists_mean = 0;
            for (size_t i = 0; i < recently_fixed_twists.size(); i++)
                recently_fixed_twists_mean += recently_fixed_twists[i];
            recently_fixed_twists_mean /= recently_fixed_twists.size();
            
            if (new_fixed)
            {
                std::cout << getTime() << " " << free_mean_twist << " " << free_min_twist << " " << free_max_twist << " " << free_twist_std << " " << free_twist_alltime_max << " " << recently_fixed_twist << " " << recently_fixed_twists_mean << std::endl;
                
//                std::cout << getTime() << " ";
//                for (int i = GetVertexIdFixedLast() + + hook_edges; i < m_rod->nv() - 2; i++)
//                {
//                    Scalar twist = m_rod->getReferenceTwist(i) + m_rod->getTheta(i) - m_rod->getTheta(i - 1);
//                    std::cout << twist << " ";
//                }
//                std::cout << std::endl;
            }
        }
    }
    
    
  // multi-step simulation, update the step
  updateMultiStep();
  
  // compute coiling radius for static coiling
  Scalar coiling_radius = 0;
  if (contact_vertices.size() > 3)
  {
    Vec2d dir_last1 = ((*(contact_vertices.rbegin() + 0) - *(contact_vertices.rbegin() + 1))).normalized();
    Vec2d dir_last2 = ((*(contact_vertices.rbegin() + 1) - *(contact_vertices.rbegin() + 2))).normalized();
    Scalar turn = dir_last1.x() * dir_last2.y() - dir_last1.y() * dir_last2.x();
    Vec2d dir_opposite;
    size_t i = 0;
    for (i = 3; i < contact_vertices.size(); i++)
    {
      dir_opposite = (*(contact_vertices.rbegin() + i - 1) - *(contact_vertices.rbegin() + i)).normalized();
      if ((dir_last1.x() * dir_opposite.y() - dir_last1.y() * dir_opposite.x()) * turn < 0)
        break;
    }
    
    if (i < contact_vertices.size())
    {
      coiling_radius = (*(contact_vertices.rbegin() + i) - *(contact_vertices.rbegin() + 1)).norm() / 2;
    }
  }  
  
  time_tick++;
  
  // fourier analysis
  int fft_interval = GetIntOpt("fft-interval");
  static int last_contactid = -1;
  int step_since_last_contactid = 0;
  static Scalar last_time = getTime();
  int contactvid = GetVertexIdFixedLast();
  int contactid = contactvid + m_numRemovedVertices;
  static int fft_count_down = 0;
  
  if (contactid != last_contactid)
  {
    Vec3d contact3d = m_rod->getVertex(contactvid);
    Vec2d contact = Vec2d(contact3d.x(), -contact3d.z());
    
    Vec3d contact_mat = m_rod->getMaterial1(contactvid > 0 ? contactvid - 1 : 0);

    if (REPORT_PD_DATA)
//    std::cerr << "new contact point: simulation time = " << getTime() << " position = " << contact3d.x() << " " << contact3d.y() << " " << contact3d.z() << std::endl;
      std::cerr << "new contact point: simulation time = " << getTime() << " position = " << contact3d.x() << " " << contact3d.y() << " " << contact3d.z() << " material frame 1 = " << contact_mat.x() << " " << contact_mat.y() << " " << contact_mat.z() << std::endl;
    if (REPORT_COILING_RADIUS)
      std::cout << coiling_radius << std::endl;

    Scalar current_time = getTime();
    step_since_last_contactid = (current_time - last_time) * 20 * m_nozzleSpeed;  // since rod resolution is 1 vertex/cm, this corresponds to ~20 samples per new contact point, enough resolution to reflect the small shift of samples earlier/later around their expected position
    
    if (contact_history.size() > 0)
    {
      Vec2d last_contact = contact_history.back();  // the vertex identified by last_contactid
      
      // push interpolation beween the last contact and the new contact into the history
      for (int i = 0; i < step_since_last_contactid; i++)
        contact_history.push_back(last_contact + (contact - last_contact) * ((Scalar)(i + 1) / step_since_last_contactid));
      
      fft_count_down -= step_since_last_contactid;
      
    } else
    {
      contact_history.push_back(contact);
    }
    
    contact_vertices.push_back(contact);
    if (contact_vertices.size() > 1000)
      contact_vertices.pop_front();
    
    step_since_last_contactid = 0;
    last_contactid = contactid;
    last_time = current_time;
    
    // perform FFT as soon as new data become available
    if ((int)contact_history.size() > fft_n)
    {
      // chop off samples that are too old
      while ((int)contact_history.size() > fft_n)
        contact_history.pop_front();
      
      if (fft_count_down < 0) // since this check doesn't happen every time step, the interval can be larger; but usually step_since_last_contactid << fft_interval so it's fine.
      {
        fft_count_down = fft_interval;
        
        // perform FFT
        // x axis
        Scalar xmean = 0;
        Scalar ymean = 0;
        for (int i = 0; i < fft_n; i++)
        {
          xmean += contact_history[i].x();
          ymean += contact_history[i].y();
        }
        xmean /= fft_n;
        ymean /= fft_n;
        
        for (int i = 0; i < fft_n; i++)
          fft_in[i] = contact_history[i].x() - xmean;
        memset(fft_out, 0, sizeof (fftw_complex) * (fft_n / 2 + 1));
        
        fftw_execute(fft_plan);

        std::vector<FFTPeak> xpeaks = fftFindPeaks(fft_out);
        
        // y axis
        for (int i = 0; i < fft_n; i++)
          fft_in[i] = contact_history[i].y() - ymean;
        memset(fft_out, 0, sizeof (fftw_complex) * (fft_n / 2 + 1));
        
        fftw_execute(fft_plan);
        
        std::vector<FFTPeak> ypeaks = fftFindPeaks(fft_out);

        // pattern matching
//        FFTPeak & xtoppeak = xpeaks[0];
//        FFTPeak & ytoppeak = ypeaks[0];
        Scalar baseomega = 0;
        Scalar baseintensity = 0;
//        if (xtoppeak.intensity > ytoppeak.intensity)
//          baseomega = xtoppeak.omega;
//        else
//          baseomega = ytoppeak.omega;
        
        // find the highest peak in both x and y peak list
        for (size_t i = 0; i < xpeaks.size(); i++)
          if (xpeaks[i].intensity > baseintensity)
            baseintensity = xpeaks[i].intensity, baseomega = xpeaks[i].omega;
        for (size_t i = 0; i < ypeaks.size(); i++)
          if (ypeaks[i].intensity > baseintensity)
            baseintensity = ypeaks[i].intensity, baseomega = ypeaks[i].omega;
        
        // normalize all the peak intensities using this top peak
        for (size_t i = 0; i < xpeaks.size(); i++)
          xpeaks[i].intensity /= baseintensity;
        for (size_t i = 0; i < ypeaks.size(); i++)
          ypeaks[i].intensity /= baseintensity;
        
        if (REPORT_PD_DATA)
        {
          for (int i = 0; i < (int)xpeaks.size(); i++)
            std::cout << "X peak intensity = " << xpeaks[i].intensity << " omega = " << xpeaks[i].omega << " phi = " << xpeaks[i].phi / pi << " pi"<< std::endl;
          for (int i = 0; i < (int)ypeaks.size(); i++)
            std::cout << "Y peak intensity = " << ypeaks[i].intensity << " omega = " << ypeaks[i].omega << " phi = " << ypeaks[i].phi / pi << " pi"<< std::endl;
        }
          
        Scalar x0 = xpeaks[0].omega;
        Scalar x1 = xpeaks[1].omega;
        Scalar y0 = ypeaks[0].omega;
        Scalar y1 = ypeaks[1].omega;
        Scalar base = baseomega;
//        std::cout << "x0 = " << x0 << " y0 = " << y0 << " base = " << base << " amx0 = " << ApproxMatch(x0, base, 2) << " amy0 = " << ApproxMatch(y0, base, 1) << std::endl;
//        if (ApproxMatch(x0, base, 2) && ApproxMatch(y0, base, 1))
//        {
//          // match: meanders
//          std::cout << "Match: meaders" << std::endl;
//        } else if (ApproxMatch(x0, base, 2) && ApproxMatch(y0, base, 1))
//        {
//          // match: translated coiling
//          std::cout << "Match: translated coiling" << std::endl;
//        } else if (ApproxMatch(x0, base, 1) && ((ApproxMatch(y0, base, 0.5) && ApproxMatch(y1, base, 1.5)) || (ApproxMatch(y0, base, 1.5) && ApproxMatch(y1, base, 0.5))))
//        {
//          // match: alternating loops
//          std::cout << "Match: alternating loops" << std::endl;
//        } else if (ApproxMatch(x0, base, 1) && ApproxMatch(x1, base, 0.5) && ApproxMatch(y0, base, 1) && ApproxMatch(y1, base, 0.5))
//        {
//          // match: double coiling
//          std::cout << "Match: double coiling" << std::endl;
//        } else if (ApproxMatch(x0, base, 0.5) && ApproxMatch(y1, base, 1))
//        {
//          // match: double meanders
//          std::cout << "Match: double meanders" << std::endl;
//        } else if (ApproxMatch(x0, base, 1) && ApproxMatch(x1, base, 2) && ApproxMatch(y0, base, 1) && ApproxMatch(y1, base, 2))
//        {
//          // match: stretched coiling
//          std::cout << "Match: stretched coiling" << std::endl;
//        } else if (ApproxMatch(x0, base, 1) && ApproxMatch(x1, base, 2) && ApproxMatch(y0, base, 2) && ApproxMatch(y1, base, 1))
//        {
//          // match: w-pattern
//          std::cout << "Match: w-pattern" << std::endl;
//        } else if (xpeaks[0].intensity < 10 && ypeaks[0].intensity < 10)  // empirical noise threshold
//        {
//          // match: catenary
//          std::cout << "Match: catenary" << std::endl;
//        } else
//        {
//          // match: disorder
//          std::cout << "Match: disorder" << std::endl;
//        }
        
        int nx = (int)xpeaks.size();
        int ny = (int)ypeaks.size();
        
        Scalar xbarycenter = 0;
        Scalar ybarycenter = 0;
        Scalar xtotalpower = 0;
        Scalar ytotalpower = 0;
        for (size_t i = 0; i < xpeaks.size(); i++)
          xbarycenter += xpeaks[i].intensity * xpeaks[i].omega,
          xtotalpower += xpeaks[i].intensity;
        xbarycenter /= xtotalpower;
        for (size_t i = 0; i < ypeaks.size(); i++)
          ybarycenter += ypeaks[i].intensity * ypeaks[i].omega,
          ytotalpower += ypeaks[i].intensity;
        ybarycenter /= ytotalpower;
        
        // this pattern recognition code is adapted from Pierre-Thomas's Mathematica code
        // whenever P-T's code doesn't agree with the paper, we stick to his code (e.g. the detection of alternating loops)
        // a lot of the checks in P-T's code are redundant (e.g. nx and ny are at least 1 by definition but there're still checks
        //  for that). however we decide to keep all these redundant checks here to ensure best agreement in the source code.
        int match = 0;
        if (nx == 1 && ny == 1 && ApproxMatch(x0, y0, 1, 0.05) && xpeaks[0].intensity * ypeaks[0].intensity > 0.01)
        {
          match = 1;  // translated coiling (P-T: Coiling)
        } else if (nx >= 2 && nx <= 4 && ny >= 2 && ny <= 4 && y0 < x1 && ApproxMatch(x1, x0, 2, 0.2) && ApproxMatch(y1, y0, 2, 0.2))
        {
          match = 2;  // double coiling (P-T: MODCoilingDouble)
        } else if (nx >= 1 && nx <= 6 && ny >= 1 && ny <= 6 && (ApproxMatch(x0, y0, 1.5, 0.2) || (ny >= 2 && ApproxMatch(y0, y1, 3, 0.2))) && xpeaks[0].intensity > 0.1 && ypeaks[0].intensity > 0.1)
        {
          match = 3;  // alternating loops (P-T: figOf8)
        } else if (nx >= 2 && nx <= 4 && ny >= 2 && ny <= 4 && x0 > x1 && y0 < y1 && ApproxMatch(x0, x1, 2, 0.2) && ApproxMatch(y1, y0, 2, 0.2) && ypeaks[1].intensity / ypeaks[0].intensity > 0.15)
        {
          match = 4;  // W-pattern (P-T: MODCoilingW)
        } else if (nx >= 1 && ny >= 1 && ApproxMatch(x0, y0, 2, 0.2) && ApproxMatch(xbarycenter, x0, 1, 0.3) && ApproxMatch(ybarycenter, y0, 1, 0.3))
        {
          match = 5;  // meandering (P-T: DOWNHeight)
        } else if (nx >= 2 && nx <= 4 && ny >= 2 && ny <= 4 && x0 > x1 && y0 > y1 && ApproxMatch(x0, x1, 2, 0.2) && ApproxMatch(y0, y1, 2, 0.2) && ypeaks[1].intensity / ypeaks[0].intensity > 0.15)
        {
          match = 6;  // stretched coiling (P-T: MODCoilingStretched)
        } else if (nx >= 1 && ny >= 1 && ApproxMatch(y0, x0, 2, 0.2))
        {
          match = 7;  // double meandering (P-T: UPHeight)
        } else if (ny >= 1 && ypeaks[0].intensity < 0.05)
        {
          match = 8;  // catenary (P-T: Catenary)
        } else if (nx >= 4 && ny >= 4)
        {
          match = 9;  // disorder (P-T: Disorder)
        } else
        {
          match = 0;  // unknown (P-T: Don't Know)
        }
          
        if (REPORT_PD_DATA)
          std::cout << "match = " << match << std::endl;
        pattern_detected.push_back(match);
        m_base_frequency = baseomega;

      }
    }
  }
  
  // cutting
//  std::cout << "v0 = " << m_rod->getVertex(0) << std::endl;
  
//	if (getTime() > 80) exit(1);

//  cout << getTime() << "\n";
//  if (getTime() >= 5) {
//    StartNewSimulation(true);
//    return;
//  }
  
	if (getTime() > 60) {
//		static bool done = false;
		
		if (!jitter_done) {
		  m_jitter.setAmplitude(GetScalarOpt("major-radius") * 0.1);
		  jitter_done = true;
		}	
	}

//  // "EPSILON SWITCHING" test. 0 -> e -> 0
//  if (is_epsilon_switch) CheckEpsilonSwitch();
//	
//	if (GetBoolOpt("inversion-test")) CheckInversion();
//	if (!is_epsilon_switch && GetBoolOpt("hopf-test")) CheckHopf();
//	
//	if (hopf_offset_test) CheckHopfOffset();
	
	
	// for drawing..
	if (1 || GetBoolOpt("inversion-test"))
	{
  
		double tw = 0;
		double len = 0;

		ElasticRod::RodForces& fcs = m_rod->getForces();
		int nv = m_rod->nv();
		RodBoundaryCondition* boundary = m_rod->getBoundaryCondition();
	
		int vfrom = 0;
		int vend = nv - 1;
 
		for (int i = vend-1; i >= vfrom; --i) {
	    if (boundary->isEdgeScripted(i)) {		
	    	vend--;
	    } else {
	    	break;
	    }
	  }
	  
		vfrom = vend - 1;
		for (int i = vend - 1; i >= 0; --i) {
	    if (boundary->isVertexScripted(i)) {
	    	break;
	    } else {
				vfrom = i;
	    }
	  }
	  
	  int nev = vend - vfrom  ;

	  if (nev <= 0) {
	  	tw = 0;
	  	len = 1.0;
	  
	  } else {
			
				  
			for(int i=0; i<fcs.size(); i++) {
				RodForce* fc = fcs[i];
        if (!fc->viscous() && dynamic_cast< RodTwistingForceSym* > (fc) != NULL) {
					RodTwistingForceSym *tf = dynamic_cast< RodTwistingForceSym* > (fc);
				
					int j=0;
					for (ElasticRod::vertex_iter vit = m_rod->vertices_begin(); vit != m_rod->vertices_end(); ++vit, ++j) {
						if (j >= vfrom && j < vend) {
							ElasticRod::vertex_handle vh = (*vit);
							double e = tf->localEnergy(vh);
							tw+=e;
							
							len += m_rod->getEdgeLength(j);
						}
					}
					
					break;
							}
			}
		}
		
    m_rod->m_twisting_energy = VecXd(m_rod->nv()+1);
    m_rod->m_twisting_energy.setZero();
    
    for(int i=0; i<fcs.size(); i++) {
      RodForce* fc = fcs[i];
      if (!fc->viscous() && dynamic_cast< RodTwistingForceSym* > (fc) != NULL) {
        RodTwistingForceSym *tf = dynamic_cast< RodTwistingForceSym* > (fc);
        
        int j=0;
        for (ElasticRod::vertex_iter vit = m_rod->vertices_begin(); vit != m_rod->vertices_end(); ++vit, ++j) {
          ElasticRod::vertex_handle vh = (*vit);
          //double e = tf->localEnergy(vh);
          double e = tf->getTwist(vh);
          
          m_rod->m_twisting_energy(j) = e;
        }
      }
    }
		msg_val = tw/len;
	}

//	if (getTime() > 130.0) exit(0);

  freezeCollidedVertices(false);
  updateFloorBoundary();
  updateContainerBoundary();

  RodBoundaryCondition* boundary = m_rod->getBoundaryCondition(); 
  int i = 0; 
  for( i = (int)m_rod->nv() - 1; i >= 0; i-- )
  {
  	if (boundary->isVertexScripted(i)) { 		
//	  	m_rod->setVelocity(i, Vec3d(0, -m_nozzleSpeed, 0));

  		Vec3d displacement = m_rod->getVelocity(i) * getDt();
  		boundary->setDesiredVertexPosition(i, (boundary->getDesiredVertexPosition(i) + displacement));
	  } else break;
  }

	int end_head = i;
	
  for( int i = 0; i < end_head; i++ )
  {
  	if (boundary->isVertexScripted(i)) { 		
	  	m_rod->setVelocity(i, Vec3d(m_beltSpeed, 0, 0));

  		Vec3d displacement = m_rod->getVelocity(i) * getDt();
  		boundary->setDesiredVertexPosition(i, (boundary->getDesiredVertexPosition(i) + displacement));

	  } 
  }    
  
  if (GetBoolOpt("floor-bc-fix-twist"))
  {
    // set edge twist constraints on the belt, after a segment of free rod of length alpha * Lgb
    Scalar alpha = GetScalarOpt("floor-bc-fix-twist-alpha");
    Scalar r = GetScalarOpt("major-radius");
    Scalar E = GetScalarOpt("youngs-modulus");
    Scalar rho = GetScalarOpt("density");
    Scalar g = GetVecOpt("gravity").norm();
    Scalar Lgb = pow(r * r * E / (8 * rho * g), 1.0 / 3.0);
    int free_edges = Lgb * alpha / GetScalarOpt("rod-length");
//    std::cout << "Lgb = " << Lgb << " free_edges = " << free_edges << std::endl;
    for (int i = 0; i < GetVertexIdFixedLast() - free_edges; i++)
    {
      if (boundary->isVertexScripted(i) && boundary->isVertexScripted(i + 1))
      {
        Scalar angle = m_rod->getTheta(i);
        
        boundary->setDesiredEdgeAngle(i, angle);
      }
    }
  }
}
//----------------------------------------------------------------------------
void SewingMachine::AfterStep()
{

  cutRod();

}
//----------------------------------------------------------------------------

void SewingMachine::Render()
{
  
}

void SewingMachine::HUDRender(int w, int h)
{  
  if (RodRenderer::lastInstance() && RodRenderer::lastInstance()->getMode() == RodRenderer::PHOTOREALISTIC)
    return;
    
  // background occlusion
  glDisable(GL_DEPTH_TEST);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glColor4f(1.0, 1.0, 1.0, 0.7);
  glBegin(GL_QUADS);
  glVertex2i(w - 105, h - 105);
  glVertex2i(w -   5, h - 105);
  glVertex2i(w -   5, h -   5);
  glVertex2i(w - 105, h -   5);
  glVertex2i(      5, h - 105);
  glVertex2i(w - 110, h - 105);
  glVertex2i(w - 110, h -   5);
  glVertex2i(      5, h -   5);
  glEnd();
  
  // axises
  glColor4f(0.0, 0.0, 0.0, 0.05);
  glBegin(GL_LINES);
  glVertex2i(w - 105, h -  55);
  glVertex2i(w -   5, h -  55);
  glVertex2i(w -  55, h - 105);
  glVertex2i(w -  55, h -   5);
  glVertex2i(      5, h -  55);
  glVertex2i(w - 110, h -  55);
  glVertex2i(w - 160, h - 105);
  glVertex2i(w - 160, h -   5);
  glEnd();
  
  // render "Motion in nozzle's frame" and "Stitch pattern on the belt" as in Brun 2012, figure 9(a)
  Scalar xrange = 17; // these numbers are taken from the belt tract renderer found in TriangleMeshRenderer::render()
  Scalar yrange = 17;
  Mat2d contact_scale;
  contact_scale << 100.0 / xrange, 0, 0, 100.0 / yrange;
  Vec2d contact_center(w - 55, h - 55);
  Vec2d belt_start(w - 160, h - 55);
  Scalar belt_speed = m_beltSpeed;
  Scalar dt = 1 / m_nozzleSpeed / 20; // contact history grows at a rate of 20 samples per centimeter of rod
  
//  glLineWidth(3);
  glBegin(GL_LINES);
  for (size_t i = 0; i < contact_history.size(); i++)
  {
    Vec2d previous = contact_scale * contact_history[i > 0 ? i - 1 : i];
    Vec2d current =  contact_scale * contact_history[i];
    Scalar belt_offset_previous = belt_speed * ((contact_history.size() - i)     * dt) * contact_scale.x();
    Scalar belt_offset_current  = belt_speed * ((contact_history.size() - i - 1) * dt) * contact_scale.x();
    
    // motion in nozzle's frame
    glColor4f(0.8, 0, 0, std::min(1.0f, (float)i * 10 / contact_history.size()));
    glVertex2f(previous.x() + contact_center.x(), previous.y() + contact_center.y());
    glVertex2f(current.x() + contact_center.x(), current.y() + contact_center.y());
    
    // stitch pattern on the belt
    glColor4f(0, 0.6, 0, 1);
    glVertex2f(previous.x() + belt_start.x() + belt_offset_previous, previous.y() + belt_start.y());
    glVertex2f(current.x()  + belt_start.x() + belt_offset_current,  current.y()  + belt_start.y());
  }
  glEnd();
  glDisable(GL_BLEND);
  
  glEnable(GL_DEPTH_TEST);
  
}
