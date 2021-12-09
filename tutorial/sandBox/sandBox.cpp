#include "sandBox.h"
#include "igl/edge_flaps.h"
#include "igl/collapse_edge.h"
#include "Eigen/dense"
#include <functional>
#include <Eigen/Core>
#include "igl/opengl/ViewerCore.h"
#include "igl/opengl/glfw/renderer.h"
#include "igl/decimate.h"
#include "igl/writeOBJ.h"
#include <igl/circulation.h>
#include <igl/shortest_edge_and_midpoint.h>
#include <igl/parallel_for.h>
#include <igl/opengl/glfw/Viewer.h>
#include <set>
#include "simplifier.h"
#include "igl/swept_volume_bounding_box.h"

SandBox::SandBox() : objectsData{ new std::vector<ObjectData*>{} }
{

}

SandBox::~SandBox()
{
    delete objectsData;
}

void SandBox::Init(const std::string &config)
{
	
	std::string item_name;
	std::ifstream nameFileout;
	nameFileout.open(config);
	if (!nameFileout.is_open())
	{
		std::cout << "Can't open file " << config << std::endl;
	}
	else
	{
        double n = 0;
		while (nameFileout >> item_name)
		{
			std::cout << "openning " << item_name << std::endl;
			load_mesh_from_file(item_name);

			parents.push_back(-1);
			data().add_points(Eigen::RowVector3d(0, 0, 0), Eigen::RowVector3d(0, 0, 1));
			data().show_overlay_depth = false;
			data().point_size = 10;
			data().line_width = 2;
			data().set_visible(false, 1);

            ObjectData* od = new ObjectData();
            InitObjectData(*od, data().V, data().F);
            objectsData->push_back(od);

            double x = -1 + (n * 2);
            double y = -1 + (n * 2);
            MoveTo(x, y);
            InitVelocity(*od);
            Eigen::RowVector3d green(0, 1, 0);
            AddBoundingBox(od->tree->m_box, green);
            data().dirty = 157; //this line prevents texture coordinates
            n++;
		}
		nameFileout.close();
	}
	MyTranslate(Eigen::Vector3d(0, 0, -1), true);
	
	data().set_colors(Eigen::RowVector3d(0.9, 0.1, 0.1));
    isActive = true;
}

void SandBox::InitObjectData(ObjectData& od, Eigen::MatrixXd& V, Eigen::MatrixXi& F)
{
    
    od.V = new Eigen::MatrixXd(V);
    od.F = new Eigen::MatrixXi(F);
    od.E = new Eigen::MatrixXi();
    od.EF = new Eigen::MatrixXi();
    od.EI = new Eigen::MatrixXi();
    od.EMAP = new Eigen::VectorXi();
    od.Q = new PriorityQueue();
    od.num_collapsed = 0;
    od.tree = new igl::AABB<Eigen::MatrixXd, 3>{};
    od.subTree = new igl::AABB<Eigen::MatrixXd, 3>{};
    od.velocity << 0, 0, 0;

    igl::edge_flaps(*od.F, *od.E, *od.EMAP, *od.EF, *od.EI);

    od.tree->init(*od.V, *od.F);

    od.C = new Eigen::MatrixXd(od.E->rows(), od.V->cols());

    ComputeNormals(od, data());

    ComputeQMatrices(od);

    ComputePriorityQueue(od);
}

void SandBox::ReInitObjectData(ObjectData& od)
{
    Eigen::MatrixXd V = *od.V; 
    Eigen::MatrixXi F = *od.F;

    // Remove old object data
    ObjectData* odToRemove = objectsData->at(selected_data_index);
    ClearObjectData(*odToRemove);
    delete odToRemove;

    // Add new object data after collapsing edges 
    ObjectData* odToAdd = new ObjectData();
    InitObjectData(*odToAdd, V, F);
    objectsData->at(selected_data_index) = odToAdd;
}

void SandBox::ClearObjectData(ObjectData& od)
{
    delete od.V;
    delete od.F;
    delete od.EMAP;
    delete od.E;
    delete od.EF;
    delete od.EI;
    delete od.Q;
    delete od.C;
    delete od.F_NORMALS;

    std::for_each(od.QMATRICES.begin(), od.QMATRICES.end(), [](Eigen::Matrix4d* m) -> void { delete m; });
    od.QMATRICES.clear();
}

void SandBox::Simplify()
{
    ObjectData* objectToSimplify = objectsData->at(selected_data_index);
    int num_to_collapse = std::ceil(0.05 * objectToSimplify->Q->size());
    Simplify(num_to_collapse, *objectToSimplify);
}

void SandBox::Simplify(int num_to_collapse, ObjectData& od)
{

    if (!od.Q->empty())
    {
        bool something_collapsed = false;
        for (int j = 0; j < num_to_collapse; j++)
        {
              if (!collapse_edge(od))
            {
                break;
            }
            something_collapsed = true;
            od.num_collapsed++;
        }
         
        if (something_collapsed)
        {
            data().clear();
            data().set_mesh(*od.V, *od.F);
            data().set_face_based(true);
            data().dirty = 157; //this line prevents texture coordinates
            ReInitObjectData(od);
        }
    }
}

void SandBox::MoveTo(double x, double y)
{
    data().TranslateInSystem(GetRotation(), Eigen::Vector3d(x, 0, 0));
    data().TranslateInSystem(GetRotation(), Eigen::Vector3d(0, y, 0));
    WhenTranslate();
}

void SandBox::InitVelocity(ObjectData& od)
{

}

void SandBox::AddBoundingBox(Eigen::AlignedBox<double, 3>& m_box, Eigen::RowVector3d& color) {
    // Corners of the bounding 
    Eigen::MatrixXd V_box(8, 3);
    V_box << m_box.corner(m_box.BottomLeftCeil).transpose(),
        m_box.corner(m_box.BottomLeftFloor).transpose(),
        m_box.corner(m_box.BottomRightCeil).transpose(),
        m_box.corner(m_box.BottomRightFloor).transpose(),
        m_box.corner(m_box.TopLeftCeil).transpose(),
        m_box.corner(m_box.TopLeftFloor).transpose(),
        m_box.corner(m_box.TopRightCeil).transpose(),
        m_box.corner(m_box.TopRightFloor).transpose();
    // Edges of the bounding box
    Eigen::MatrixXi E_box(12, 2);
    E_box <<
        0, 1,
        1, 3,
        2, 3,
        2, 0,
        4, 5,
        5, 7,
        6, 7,
        6, 4,
        0, 4,
        1, 5,
        2, 6,
        7, 3;
    // Plot the corners of the bounding box as points
    data().add_points(V_box, color);
    // Plot the edges of the bounding box
    for (unsigned i = 0; i < E_box.rows(); ++i)
    {
        data().add_edges(V_box.row(E_box(i, 0)), V_box.row(E_box(i, 1)), color);
    }
}

void SandBox::Animate()
{
	if (isActive)
	{
        data().MyTranslate(Eigen::Vector3d(-0.005, 0, 0), true);
	}
}
