#pragma once

#include "igl/opengl/glfw/Display.h"
#include "igl/opengl/ViewerData.h"
#include "igl/edge_flaps.h"
#include "igl/collapse_edge.h"
#include "set"
#include "DataObject.h"

//#define SEMP

void ComputePriorityQueue(ObjectData& od);
void ComputeNormals(ObjectData& od, igl::opengl::ViewerData& viewerData);
std::vector<int> ComputeVertexFaces(const Eigen::MatrixXi& F, const int v);
void ComputeQMatrices(ObjectData& od);
Eigen::Matrix4d ComputeQMatrix(ObjectData& od, int v);
double ComputeCost(Eigen::RowVectorXd& v, Eigen::Matrix4d& Q_Matrix);
Eigen::RowVectorXd ComputePlace(ObjectData& od, Eigen::Matrix4d& _Q, int v1, int v2);
Eigen::Matrix4d ComputeKp(ObjectData& od, Eigen::RowVectorXd& plane, int v);
bool collapse_edge(ObjectData& od);
