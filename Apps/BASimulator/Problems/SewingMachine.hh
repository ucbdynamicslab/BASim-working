// rods/src/problems/SewingMachine.hh
// author: Miklos Bergou (miklos@cs.columbia.edu
// December 19, 2008

#ifndef SEWING_MACHINE_HH
#define SEWING_MACHINE_HH

#include "ProblemBase.hh"

#include "BASim/src/Util/KeyedValue.hh"
#include "BASim/src/Physics/ElasticRods/BridsonStepper.hh"
#include "JitterFunction.hh"

#include <fftw3.h>
#include <queue>

class SewingMachine : public Problem
{
public:
  SewingMachine();
  virtual ~SewingMachine();

protected:

  void Setup();
  void AtEachTimestep();
	void AfterStep();
  void Render();
  void HUDRender(int w, int h);
  
  ElasticRod* m_rod;
//  RodTimeStepper* stepper;

  Scalar m_maxTwist;
  Scalar m_twistRate;
  Scalar m_currentTwist;
  
  
private:

	void loadAuxOptions();
	void setupTimeStepper();	
	void makeRod();
	void makeBelt();
	
	void CheckInversion();	
	void CheckHopf();
//	void CheckEpsilonSwitch();
	
  bool CheckSteadyState(int n_verts_skip);
	
//	void makeRodFromFile();
	
	int GetVertexIdFixedLast();
	int GetExtremePoints(std::vector<Vec2d> &ext_points, const int n_points);
  int ComputeAmplitudeAndWaveLength(std::vector<double> &x_positions, double &mean_amp, std::vector<double> &amplitudes, double &mean_wavelength, std::vector <double>& wavelengths, const int n_peaks);
  
//  double RecordAmplitudes();
//  void RecordAllPoints();
	
	
//	void CheckHopfOffset();
	
  Scalar m_rodLength;
	Scalar m_fallHeight;
  Scalar m_beltSpeed;
  Scalar m_nozzleSpeed;
  Scalar m_nozzleSpeedInit;

  int m_headVertices;
  int m_tailVertices;
  int m_numRemovedVertices;
  
  void updateFloorBoundary();
  void updateContainerBoundary();

  void freezeCollidedVertices(bool saveToFile);

  void removeVertex();
  void cutRod();

//  void StartNewSimulation(bool is_zero_Ao , double final_Ao);

  void measureLength(Scalar &l_whole, Scalar &l_hanging, Scalar &l_belt, Scalar &ol_whole, Scalar &ol_hanging, Scalar &ol_belt);
  
  void RecordLength();
  
//  void CheckStaticCoils(bool inversion);
  
//  void GetTextFileName(std::string &file_name, const std::string &header, const std::string &ext);
  
  // FFT utilities
  struct FFTPeak
  {
    Scalar intensity; // strength
    Scalar omega;     // frequency (in cycles per fft_n samples, i.e. in the window) estimated by interpolating the neighboring sample
    Scalar phi;       // phase
    
    bool operator < (const FFTPeak & p) const { return intensity > p.intensity; }
//    bool operator < (const FFTPeak & p) const { return omega < p.omega; }
  };
  
  std::vector<FFTPeak> fftFindPeaks(fftw_complex * out);
  bool ApproxMatch(Scalar omega, Scalar omegabase, int multiplier, Scalar tolerance = 0.1);
  
  // multistep simulation
  void updateMultiStep();
  bool steadyStateReached();
    
  // after switching epsilon back to zero
  Scalar cut_length;
  Scalar cut_length_udf;
  std::vector <double> length_info;
  int length_info_size;

  std::vector <double> contactp_info;

  JitterFunction m_jitter;
  int topVert;
  
  BridsonStepper* m_br_stepper;
	TriangleMesh* tri_mesh;
	
	int m_contact_method;
	int m_friction_method;
  
  std::vector<double> m_inv_times;
  
  std::vector<Vec3d> m_cut_pos;
  std::vector<double> m_cut_time;

  // List of positions of vertices fixed on the belt. To avoid updating positions at each time step, it records position w.r.t belt at t = 0.
  // ex) If a vertex is fixed at (1, 2) at t = 3.0 and the belt speed = 2.0, recorded position would be (1 - 3.0*2.0, 2)
//  std::vector<Vec3d> m_fixed_position;
  
  // The number of fixed vertices excluding the part which has been deleted. It only cares about vertices being simulated.
  // To get the number of all fixed vertices including removed part, use m_fixed_position.size()
//  int m_nFixed;
  
  
  bool is_epsilon_switch;
  int epsilon_switch_step;
  double epsilon_switch_time[6];
  int epsilon_switch_early_verts;
  int epsilon_switch_early2_verts;
  int epsilon_switch_main_verts;
  double epsilon_switch_t[5];
//  std::vector<Vec3d> m_all_fixed_pos;

  // move some static values...
  int hopf_msg_time;
  bool jitter_done;
  int time_tick;

  
  // Test for finding the critical curvature value at the boundary of zero/nonzero A_o (final amplitude)
  bool hopf_eps_switch_varying_curv;
  
  // The current curvature
  double hopf_crt_curv;

  // Maximum curvature value which has been confirmed to yield zero A_o
  double hopf_max_zero_curv;

  // Minimum curvature value which has been confirmed to yield non-zero A_o
  double hopf_min_nonzero_curv;
  
  // Stage. 
  // 0 : increasing curvature phase by "hopf_varying_curv_jump". The next value will be set to current value + "hopf_varying_curv_jump"
  // 1 : After hitting any non-zero A_o. The next value will be set to the mid-point of current value & (max_zero or min_nonzero)
  // ex) Initial curvature = 0.02 & hopf_varying_curv_jump = 0.01. Then it starts with 0.02 and keep increasing 0.03, 0.04, ... until it meets a non-zero A_o.
  //     Then 0.07, 0.065, 0.0625, .... to refine the value.
  int hopf_varying_step;
  
  double hopf_varying_curv_jump;
  
  double youngM;
  double shearM;
  

  // records for multiple simulations
  int multi_n;
  std::vector<double> multi_c;  
  std::vector<bool> multi_zeroA;
  std::vector<double> multi_finalTime;
  std::vector<double> multi_finalA;
  
  bool hopf_offset_test;
  double hopf_offset_endtime;
  double hopf_offset_length;

  // for FFT analysis of the contact point
  int fft_n;
  double * fft_in;
  fftw_complex * fft_out;
  fftw_plan fft_plan;
  
  std::deque<Vec2d> contact_history;
  std::vector<int> pattern_detected;
  std::deque<Vec2d> contact_vertices;
  
  // for multistep simulation
  int m_step; // goes from 0 to nstep
  Scalar m_step_start_time;
  
  Scalar m_base_frequency;
  
};

#endif  // SEWING_MACHINE_HH
