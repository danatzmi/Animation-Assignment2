#pragma once

#include "igl/opengl/glfw/Display.h"
#include "igl/opengl/ViewerData.h"
#include "igl/edge_flaps.h"
#include "igl/collapse_edge.h"
#include "set"

//#define SEMP

typedef std::set<std::pair<double, int>> PriorityQueue;

struct _ObjectData {
	Eigen::MatrixXd* V;
	Eigen::MatrixXi* F;
	Eigen::VectorXi* EMAP;
	Eigen::MatrixXi* E;
	Eigen::MatrixXi* EF;
	Eigen::MatrixXi* EI;
	PriorityQueue* Q;
	std::vector<PriorityQueue::iterator> Qit;
	Eigen::MatrixXd* C;
	Eigen::MatrixXd* F_NORMALS;
	std::vector<Eigen::Matrix4d*> QMATRICES;
	int num_collapsed;
} typedef ObjectData;

void ComputePriorityQueue(ObjectData& od);
void ComputeNormals(ObjectData& od, igl::opengl::ViewerData& viewerData);
std::vector<int> ComputeVertexFaces(const Eigen::MatrixXi& F, const int v);
void ComputeQMatrices(ObjectData& od);
Eigen::Matrix4d ComputeQMatrix(ObjectData& od, int v);
double ComputeCost(Eigen::RowVectorXd& v, Eigen::Matrix4d& Q_Matrix);
Eigen::RowVectorXd ComputePlace(ObjectData& od, Eigen::Matrix4d& _Q, int v1, int v2);
Eigen::Matrix4d ComputeKp(ObjectData& od, Eigen::RowVectorXd& plane, int v);
bool collapse_edge(ObjectData& od);
void update_priority_queue(ObjectData& od);
