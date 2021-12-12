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
	void OnNewMeshLoad();
	void Init(const std::string& config);
	void InitObjectData(ObjectData& od, Eigen::MatrixXd& V, Eigen::MatrixXi& F);
	void ReInitObjectData(ObjectData& od);
	void ClearObjectData(ObjectData& od);
	void Simplify();
	void Simplify(int num_to_collapse, ObjectData& od);
	void MoveTo(double x, double y);
	bool ObjectsCollide(igl::AABB<Eigen::MatrixXd, 3>* firstTree, igl::AABB<Eigen::MatrixXd, 3>* secondTree);
	bool BoxesIntersect(Eigen::AlignedBox <double, 3>& firstBox, Eigen::AlignedBox <double, 3>& secondBox);
	void SetVelocity(double x, double y);

private:
	// Prepare array-based edge data structures and priority queue
	std::vector<ObjectData*>* objectsData;
	double xVelocity;
	double yVelocity;
	void Animate();
};

