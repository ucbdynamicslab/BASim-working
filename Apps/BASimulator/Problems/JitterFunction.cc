#include "JitterFunction.hh"
#include <cstdlib>

JitterFunction::JitterFunction(int n, double amplitude, double frequency)
  : m_a(amplitude)
  , m_w(frequency)
{
  for (int i = 0; i < n; ++i) {
//    m_noise.push_back(new noise::module::Perlin());
    m_directions.push_back(Vec3d::Random().normalized());
  }
}

JitterFunction::~JitterFunction()
{
  //for (size_t i = 0; i < m_noise.size(); ++i) {
//    delete m_noise[i];
  //}
}

void JitterFunction::getValues(std::vector<double>& values, double t)
{
	for(int i=0; i< values.size(); i++) {
		values[i] = m_a * (double)rand() / (double)RAND_MAX * m_directions[0](i);
//    values[i] = 0;
	}
	
//  for (size_t i = 0; i < m_noise.size(); ++i) {
//    Vec3d x = t * m_w * m_directions[i];
//    values[i] = m_a * m_noise[i]->GetValue(x[0], x[1], x[2]);
//  }
}
