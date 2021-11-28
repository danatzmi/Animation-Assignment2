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
		}
		nameFileout.close();
	}
	MyTranslate(Eigen::Vector3d(0, 0, -1), true);
	
	data().set_colors(Eigen::RowVector3d(0.9, 0.1, 0.1));
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

    igl::edge_flaps(*od.F, *od.E, *od.EMAP, *od.EF, *od.EI);

    od.C = new Eigen::MatrixXd(od.E->rows(), od.V->cols());

    ComputeNormals(od, data());

    ComputeQMatrices(od);

    ComputePriorityQueue(od);
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
            ComputeNormals(od, data());
            ComputeQMatrices(od);
            ComputePriorityQueue(od);
        }
         
        if (something_collapsed)
        {
            data().clear();
            data().set_mesh(*od.V, *od.F);
            data().set_face_based(true);
            data().dirty = 157; //this line prevents texture coordinates
            InitObjectData(od, *od.V, *od.F);
        }
    }
}

void SandBox::Animate()
{
	if (isActive)
	{
		
		
		
	}
}
