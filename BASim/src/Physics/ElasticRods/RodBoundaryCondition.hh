/**
 * \file RodBoundaryCondition.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 12/09/2009
 */

#ifndef RODBOUNDARYCONDITION_HH
#define RODBOUNDARYCONDITION_HH

#include "BASim/src/Core/ObjectHandles.hh"
#include "BASim/src/Core/TopologicalObject/TopObjHandles.hh"

namespace BASim {

class ElasticRod;

/** Class for managing fixed/scripted vertices and edges of a rod */
// TODO: For consistency, methods in RodBoundaryCondition should acccept
//       edge and/or vertex iterators to be consistent with the rest of
//       BASim.
class RodBoundaryCondition
{
public:

  typedef std::vector<int> BCList;

  RodBoundaryCondition(ElasticRod& rod);

  const BCList& scriptedVertices() const;
  bool isVertexScripted(int vertIdx) const;
  void setDesiredVertexPosition(int vertIdx, const Vec3d& position);
  const Vec3d& getDesiredVertexPosition(int vertIdx);
  void releaseVertex(int vertIdx);
  
  void adjustListAfterDeleteVertex(int vertIdx);

  const BCList& scriptedEdges() const;
  bool isEdgeScripted(int edgeIdx) const;
  void setDesiredEdgeAngle(int edgeIdx, const Scalar& theta);
  const Scalar& getDesiredEdgeAngle(int edgeIdx);
  void releaseEdge(int edgeIdx);

  void adjustListAfterDeleteEdge(int edgeIdx);
  
  void clearVerticalConstraints() 
  {
  	scriptedVerticalVertices.clear();
  	desiredVerticalPositions.clear();
  	desiredVerticalVelocities.clear();
  	unconstrainedVerticalPositions.clear();
  	isStaticMode.clear();
  }
  
  void setVerticalConstraint(int vertex_id, Scalar desired_pos, Scalar desired_vel, Scalar uncons_pos, bool is_static)
  {
  	scriptedVerticalVertices.push_back(vertex_id);
  	desiredVerticalPositions.push_back(desired_pos);
  	desiredVerticalVelocities.push_back(desired_vel);
  	unconstrainedVerticalPositions.push_back(uncons_pos);
  	isStaticMode.push_back(is_static);
  }
  
  int getNumberOfVerticallyScriptedVertices() const {
  	return (int)scriptedVerticalVertices.size();
  }
  
  // 'id' means the order in the vector of scripted vertices
  int getVerticalScriptedVertId(int id) const
  {
  	return scriptedVerticalVertices[id];
  }

  Scalar getDesiredVerticalPosition(int id) const
  {
  	return desiredVerticalPositions[id];
  }

  Scalar getDesiredVerticalVelocity(int id) const
  {
  	return desiredVerticalVelocities[id];
  }
  
  Scalar getUnconstrainedVerticalPosition(int id) const
  {
  	return unconstrainedVerticalPositions[id];
  }
  
  bool isStatic(int id) const
  {
  	return isStaticMode[id];
  }  

  void setStatic(int id, bool is_static) 
  {
  	isStaticMode[id] = is_static;
  }  
  

protected:

  ElasticRod& m_rod;

  ObjPropHandle<BCList> m_scriptedVerts;
  VPropHandle<Vec3d> m_desiredPositions;
  VPropHandle<bool> m_isVertexScripted;

  ObjPropHandle<BCList> m_scriptedEdges;
  EPropHandle<Scalar> m_desiredTheta;
  EPropHandle<bool> m_isMaterialScripted;

  std::vector<int> scriptedVerticalVertices;
  std::vector<Scalar> desiredVerticalPositions;
  std::vector<Scalar> desiredVerticalVelocities;
  std::vector<Scalar> unconstrainedVerticalPositions;
  std::vector<bool> isStaticMode;
  
};

} // namespace BASim

#endif // RODBOUNDARYCONDITION_HH
