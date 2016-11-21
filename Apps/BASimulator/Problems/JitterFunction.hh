/**
 * \file JitterFunction.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 01/16/2010
 */

#ifndef JITTERFUNCTION_HH
#define JITTERFUNCTION_HH

#include "Option.hh"
//#include <libnoise/noise.h>
#include <vector>


/// Map from R to R^n, where the output will be have slowly varying
/// values.
class JitterFunction
{
public:

  JitterFunction(int n, double amplitude = 1, double frequency = 1);
  ~JitterFunction();

  void setAmplitude(double a) { m_a = a; }
  double getAmplitude() const { return m_a; }

  void setFrequency(double w) { m_w = w; }
  double getFrequency() const { return m_w; }

  void getValues(std::vector<double>& values, double t);

private:

//  std::vector<noise::module::Perlin*> m_noise;
  double m_a;
  double m_w;

  std::vector<Vec3d> m_directions;
};

#endif // JITTERFUNCTION_HH
