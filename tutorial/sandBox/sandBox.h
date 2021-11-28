#pragma once
#include "igl/opengl/glfw/Viewer.h"
#include "simplifier.h"

typedef std::set<std::pair<double, int>> PriorityQueue;

class SandBox : public igl::opengl::glfw::Viewer
{
public:
	SandBox();
	~SandBox();
	void Init(const std::string& config);
	void InitObjectData(ObjectData& od, Eigen::MatrixXd& V, Eigen::MatrixXi& F);
	void Simplify();
	void Simplify(int num_to_collapse, ObjectData& od);

private:
	// Prepare array-based edge data structures and priority queue
	std::vector<ObjectData*>* objectsData;
	void Animate();
};

