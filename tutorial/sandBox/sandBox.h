#pragma once
#include "igl/opengl/glfw/Viewer.h"
#include "igl/AABB.h"
#include "simplifier.h"

typedef std::set<std::pair<double, int>> PriorityQueue;

class SandBox : public igl::opengl::glfw::Viewer
{
public:
	SandBox();
	~SandBox();
	void Init(const std::string& config);
	void InitObjectData(ObjectData& od, Eigen::MatrixXd& V, Eigen::MatrixXi& F);
	void ReInitObjectData(ObjectData& od);
	void ClearObjectData(ObjectData& od);
	void Simplify();
	void Simplify(int num_to_collapse, ObjectData& od);
	void MoveTo(double x, double y);
	void InitVelocity(ObjectData& od);
	void AddBoundingBox(Eigen::AlignedBox<double, 3>& m_box, Eigen::RowVector3d& color);


private:
	// Prepare array-based edge data structures and priority queue
	std::vector<ObjectData*>* objectsData;
	
	void Animate();
};

