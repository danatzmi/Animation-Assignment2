#pragma once

#include "set"
#include "igl/AABB.h"

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
	igl::AABB<Eigen::MatrixXd, 3>* tree;
} typedef ObjectData;
