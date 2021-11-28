#include "simplifier.h"
#include <igl/shortest_edge_and_midpoint.h>

void ComputePriorityQueue(ObjectData& od)
{

#ifdef SEMP
    const auto& cost_and_placement = igl::shortest_edge_and_midpoint;
#else
    const auto& cost_and_placement = [&](const int e,
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& F,
        const Eigen::MatrixXi& E,
        const Eigen::VectorXi& EMAP,
        const Eigen::MatrixXi& EF,
        const Eigen::MatrixXi& EI,
        double& cost,
        Eigen::RowVectorXd& p
        )
    {
        int v1 = E(e, 0);
        int v2 = E(e, 1);

        Eigen::Matrix4d* Q1 = od.QMATRICES[v1];
        Eigen::Matrix4d* Q2 = od.QMATRICES[v2];

        Eigen::Matrix4d* _Q = new Eigen::Matrix4d(*Q1 + *Q2);

        p = ComputePlace(od, *_Q, v1, v2); //0.5 * (V.row(E(e, 0)) + V.row(E(e, 1))); //

        cost = ComputeCost(p, *_Q); //(V.row(E(e, 0)) - V.row(E(e, 1))).norm(); //
    };
#endif

    od.Qit.clear();
    od.Qit.resize(od.E->rows());
    od.Q->clear();
    for (int e = 0; e < od.E->rows(); e++)
    {
        if ((*(od.E))(e, 0) == IGL_COLLAPSE_EDGE_NULL && (*(od.E))(e, 1) == IGL_COLLAPSE_EDGE_NULL)
        {
            continue;
        }
        double cost = e;
        Eigen::RowVectorXd p(1, 3);
        cost_and_placement(e, *od.V, *od.F, *od.E, *od.EMAP, *od.EF, *od.EI, cost, p);
        od.C->row(e) = p;
        od.Qit[e] = od.Q->insert(std::pair<double, int>(cost, e)).first;
    }
}

void ComputeNormals(ObjectData& od, igl::opengl::ViewerData& viewerData)
{
    od.F_NORMALS = new Eigen::MatrixXd();
    od.F_NORMALS->resize(od.F->rows(), 3);
    for (int i = 0; i < od.F_NORMALS->rows(); i++)
    {
        od.F_NORMALS->row(i) = viewerData.F_normals.row(i);
    }
}

std::vector<int> ComputeVertexFaces(const Eigen::MatrixXi& F, const int v)
{
    std::vector<int> faces{};
    for (int f = 0; f < F.rows(); f++) {

        if (F(f, 0) == IGL_COLLAPSE_EDGE_NULL && F(f, 1) == IGL_COLLAPSE_EDGE_NULL && F(f, 2) == IGL_COLLAPSE_EDGE_NULL) {
            continue;
        }

        if (F(f, 0) == v || F(f, 1) == v || F(f, 2) == v) {
            faces.push_back(f);
        }
    }

    return faces;
}

void ComputeQMatrices(ObjectData& od)
{
    for (int v = 0; v < od.V->rows(); v++)
    {
        od.QMATRICES.push_back(new Eigen::Matrix4d(ComputeQMatrix(od, v)));
    }
}

Eigen::Matrix4d ComputeQMatrix(ObjectData& od, int v)
{
    Eigen::Matrix4d Kp = Eigen::Matrix4d::Zero();
    std::vector<int> faces = ComputeVertexFaces(*od.F, v);
    for (int f : faces)
    {
        Eigen::RowVectorXd plane1 = od.F_NORMALS->row(f).normalized();
        Eigen::RowVectorXd plane(3);
        plane << plane1(0), plane1(1), plane1(2); 

        Kp += ComputeKp(od, plane, v);
    }
    return Kp;
}

double ComputeCost(Eigen::RowVectorXd& v, Eigen::Matrix4d& Q_Matrix)
{
    Eigen::RowVector4d vertex_expaneded;
    vertex_expaneded << v(0), v(1), v(2), 1;
    Eigen::Matrix<double, 1, 1> res = vertex_expaneded * Q_Matrix * vertex_expaneded.transpose();
    return res(0);
}

Eigen::RowVectorXd ComputePlace(ObjectData& od, Eigen::Matrix4d& _Q, int v1, int v2)
{
    Eigen::RowVectorXd place(3);
    Eigen::RowVector3d _v1 = od.V->row(v1);
    Eigen::RowVector3d _v2 = od.V->row(v2);
    //return (_v1 + _v2) / 2;
    Eigen::Matrix4d Qtmp(_Q);
    Qtmp(3, 0) = 0;
    Qtmp(3, 1) = 0;
    Qtmp(3, 2) = 0;
    Qtmp(3, 3) = 1;
    if (Qtmp.determinant() != 0) // Matrix is invertible
    {
        Eigen::Matrix4d _place = Qtmp.inverse(); // *Eigen::Vector4d{ 0, 0, 0, 1 };
        place << _place(0, 3), _place(1, 3), _place(2, 3);
    }
    else // Matrix is not invertible - we take midpoint
    {
        place = (_v1 + _v2) / 2;
    }
    return place;
}

Eigen::Matrix4d ComputeKp(ObjectData& od, Eigen::RowVectorXd& plane, int v)
{
    double a = plane(0);
    double b = plane(1);
    double c = plane(2);
	
    Eigen::RowVector3d _v = od.V->row(v);
    double x = _v(0), y = _v(1), z = _v(2);
    double d = 0 - (a * x) - (b * y) - (c * z);

    double a2 = a * a;
    double ab = a * b;
    double ac = a * c;
    double ad = a * d;
    double b2 = b * b;
    double bc = b * c;
    double bd = b * d;
    double c2 = c * c;
    double cd = c * d;
    double d2 = d * d;

    Eigen::Matrix4d res;
    res <<  a2, ab, ac, ad,
            ab, b2, bc, bd,
            ac, bc, c2, cd,
            ad, bd, cd, d2;

    return res;
}

bool collapse_edge(ObjectData& od)
{

#ifdef SEMP
    const auto& cost_and_placement = igl::shortest_edge_and_midpoint;
#else
    const auto& cost_and_placement = [&](const int e,
        const Eigen::MatrixXd & V,
        const Eigen::MatrixXi & F,
        const Eigen::MatrixXi & E,
        const Eigen::VectorXi & EMAP,
        const Eigen::MatrixXi & EF,
        const Eigen::MatrixXi & EI,
        double& cost,
        Eigen::RowVectorXd & p
        )
    {
        int v1 = E(e, 0);
        int v2 = E(e, 1);

        Eigen::Matrix4d* Q1 = od.QMATRICES[v1];
        Eigen::Matrix4d* Q2 = od.QMATRICES[v2];

        Eigen::Matrix4d* _Q = new Eigen::Matrix4d(*Q1 + *Q2);

        p = ComputePlace(od, *_Q, v1, v2); //0.5 * (V.row(E(e, 0)) + V.row(E(e, 1))); //

        cost = ComputeCost(p, *_Q); //(V.row(E(e, 0)) - V.row(E(e, 1))).norm(); //
    };
#endif

    int e = od.Q->begin()->second;
    double cost = od.Q->begin()->first;
    Eigen::RowVector3d p = od.C->row(e);

#if 1
    if (igl::collapse_edge(cost_and_placement,
        *od.V,
        *od.F,
        *od.E,
        *od.EMAP,
        *od.EF,
        *od.EI,
        *od.Q,
        od.Qit,
        *od.C))
    {
        std::cout << "edge " << e 
                << ", cost = " << cost 
                << ", new v position (" << p(0) << ", " << p(1) << ", " << p(2) << ")" 
                << std::endl;
        return true;
    }
    return false;
#endif

}

void update_priority_queue(ObjectData& simplifyDataObject)
{
    PriorityQueue::iterator it;
    PriorityQueue tmp;
    PriorityQueue new_priority_queue;
    for (it = simplifyDataObject.Q->begin(); it != simplifyDataObject.Q->end(); ++it) {
        int e = it->second;
        if (it->first == std::numeric_limits<double>::infinity()) {
            tmp.insert(std::pair<double, int>(it->first, e));
            continue;
        }
        int v1 = (*simplifyDataObject.E)(e, 0);
        int v2 = (*simplifyDataObject.E)(e, 1);
        Eigen::Matrix4d* q_matrix_v1 = simplifyDataObject.QMATRICES[v1];
        Eigen::Matrix4d* q_matrix_v2 = simplifyDataObject.QMATRICES[v2];
        Eigen::Matrix4d* qq = new Eigen::Matrix4d(*q_matrix_v1 + *q_matrix_v2);
        Eigen::RowVectorXd p = ComputePlace(simplifyDataObject, *qq, v1, v2);
        simplifyDataObject.C->row(e) = p;

        double cost = ComputeCost(p, *qq);
        tmp.insert(std::pair<double, int>(cost, e));
    }
    delete simplifyDataObject.Q;
    simplifyDataObject.Q = new PriorityQueue();
    for (it = tmp.begin(); it != tmp.end(); ++it) {
        simplifyDataObject.Q->insert(std::pair<double, int>(it->first, it->second));
    }
}