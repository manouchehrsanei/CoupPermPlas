
#include <cmath>

#include <iostream>
#include <fstream>
#include <string>
#include "pzgmesh.h"
#include "TPZVTKGeoMesh.h"
#include "pzanalysis.h"
#include "pzbndcond.h"

#include "TPZGmshReader.h"


#include "TPZCPPDarcyMat.h"
#include "TPZCPPDarcyWithMem.h"

#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"
#include <pzgeoel.h>
#include "pzgeoelbc.h"


using namespace std;

//---------------------------------Model-------------------------------------

// Define Geometry
const int  dim            =  2;
int       pOrder          =  2;
int      h_level          =  8;
double     lx             =  1.0;
double     ly             =  1.0;
int       nelx            =  h_level;
int       nely            =  h_level;
int        nx             =  nelx + 1;
int        ny             =  nely + 1;

// Define id for a material (weak formulation)
int     matid             =  1;
int     m_matBCbott       = -1;
int     m_matBCtop        = -2;
int     m_matBCleft       = -3;
int     m_matBCright      = -4;
int     m_matPoint        = -5;

// Define Boundary condition
const int dirichlet       =  0;
//const int neumann         =  1;

// Define post processing resolution
int postProcessResolution = 0;

// @brief Function to create the geometric mesh
TPZGeoMesh * ReadGeometry();
void PrintGeometry(TPZGeoMesh * gmesh);

// @brief Function to create the geometric mesh
TPZGeoMesh *CreateGMesh(int nelx, int nely, double hx, double hy);

// @brief Function to create  the computational mesh
TPZCompMesh *CMesh(TPZGeoMesh *gmesh, int pOrder);

// @brief Function of analysis
TPZAnalysis * Analysis(TPZCompMesh * cmesh);

// Post-process utilities
void PostProcess(TPZAnalysis *an);



// Main function of the program:
int main()
{
    
    TPZGeoMesh *gmesh = CreateGMesh(nx, ny, lx, ly); // Function to create geometry
    
    TPZCompMesh *cmesh = CMesh(gmesh, pOrder); // Function to create polynomial mesh

    TPZAnalysis * analysis = Analysis(cmesh);
    
    analysis->Assemble();
    analysis->Solve();

    PostProcess(analysis);
    std::cout << "Execution complete." << std::endl;
    return 0;
    
    
}


// Create geometry
TPZGeoMesh *CreateGMesh(int nx, int ny, double lx, double ly)
{
    
    int64_t id, index;
    int dim = 2;
    
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    gmesh->SetDimension(dim);
    
    TPZVec <REAL> coord (3,0.);
    
    // Initialization of nodes:
    for(int i = 0; i < ny; i++)
    {
        for(int j = 0; j < nx; j++)
        {
            id = i*nx + j;
            coord[0] = (j)*lx/(nx - 1);
            coord[1] = (i)*ly/(ny - 1);
            coord[2] = 0.;
            //Get the index of the node
            index = gmesh->NodeVec().AllocateNewElement();
            //Set the value of the node in the mesh nodes vector
            gmesh->NodeVec()[index] = TPZGeoNode(id,coord,*gmesh);
        }
    }
    
    // Point 1
    TPZVec<int64_t> pointtopology(1);
    pointtopology[0] = 0;
    
    gmesh->CreateGeoElement(EPoint,pointtopology,m_matPoint,id);
    
    // Auxiliary vector to store the connections between elements:
    TPZVec <int64_t> connect(4,0);
    
    // Connectivity of elements:
    for(int i = 0; i < (ny - 1); i++)
    {
        for(int j = 0; j < (nx - 1); j++)
        {
            index      = (i)*(nx - 1)+ (j);
            connect[0] = (i)*ny + (j);
            connect[1] = connect[0]+1;
            connect[2] = connect[1]+(nx);
            connect[3] = connect[0]+(nx);
            gmesh->CreateGeoElement(EQuadrilateral,connect,matid,id);
        }
    }
    
    
    // Generating neighborhood information:
    gmesh->BuildConnectivity();
    
    {
        TPZCheckGeom check(gmesh);
        check.CheckUniqueId();
    }
    
    int64_t el, numelements = gmesh->NElements();
    TPZManVector <int64_t> TopolPlate(4);
    
    for (el=0; el<numelements; el++)
    {
        int64_t totalnodes = gmesh->ElementVec()[el]->NNodes();
        TPZGeoEl *plate = gmesh->ElementVec()[el];
        for (int i=0; i<4; i++)
        {
            TopolPlate[i] = plate->NodeIndex(i);
        }
        
        // The boundary conditions:
        TPZManVector <TPZGeoNode> Nodefinder(totalnodes);
        TPZManVector <REAL,3> nodecoord(3);
        
        // In face x = 1
        TPZVec<int64_t> ncoordzbottVec(0);
        TPZVec<int64_t> ncoordztopVec(0);
        TPZVec<int64_t> ncoordzleftVec(0);
        TPZVec<int64_t> ncoordzrightVec(0);
        
        int64_t sizeOfbottVec = 0;
        int64_t sizeOftopVec = 0;
        int64_t sizeOfleftVec = 0;
        int64_t sizeOfrightVec = 0;
        
        for (int64_t i = 0; i < totalnodes; i++)
        {
            Nodefinder[i] = gmesh->NodeVec()[TopolPlate[i]];
            Nodefinder[i].GetCoordinates(nodecoord);
            if (nodecoord[2] == 0. & nodecoord[1] == 0.)
            {
                sizeOfbottVec++;
                ncoordzbottVec.Resize(sizeOfbottVec);
                ncoordzbottVec[sizeOfbottVec-1] = TopolPlate[i];
            }
            if (nodecoord[2] == 0. & nodecoord[1] == ly)
            {
                sizeOftopVec++;
                ncoordztopVec.Resize(sizeOftopVec);
                ncoordztopVec[sizeOftopVec-1] = TopolPlate[i];
            }
            if (nodecoord[2] == 0. & nodecoord[0] == 0.)
            {
                sizeOfleftVec++;
                ncoordzleftVec.Resize(sizeOfleftVec);
                ncoordzleftVec[sizeOfleftVec-1] = TopolPlate[i];
            }
            if (nodecoord[2] == 0. & nodecoord[0] == lx)
            {
                sizeOfrightVec++;
                ncoordzrightVec.Resize(sizeOfrightVec);
                ncoordzrightVec[sizeOfrightVec-1] = TopolPlate[i];
            }
        }
        
        if (sizeOfbottVec == 2)
        {
            int sidesbott = plate->WhichSide(ncoordzbottVec);
            TPZGeoElSide platesidebott(plate, sidesbott);
            TPZGeoElBC(platesidebott,m_matBCbott);
        }
        
        if (sizeOftopVec == 2)
        {
            int sidestop = plate->WhichSide(ncoordztopVec);
            TPZGeoElSide platesidetop(plate, sidestop);
            TPZGeoElBC(platesidetop,m_matBCtop);
        }
        
        if (sizeOfleftVec == 2)
        {
            int sidesleft = plate->WhichSide(ncoordzleftVec);
            TPZGeoElSide platesideleft(plate, sidesleft);
            TPZGeoElBC(platesideleft,m_matBCleft);
        }
        
        if (sizeOfrightVec == 2)
        {
            int sidesright = plate->WhichSide(ncoordzrightVec);
            TPZGeoElSide platesideright(plate, sidesright);
            TPZGeoElBC(platesideright,m_matBCright);
        }
        
        ncoordzbottVec.Resize(0);
        sizeOfbottVec = 0;
        ncoordztopVec.Resize(0);
        sizeOftopVec = 0;
        ncoordzleftVec.Resize(0);
        sizeOfleftVec = 0;
        ncoordzrightVec.Resize(0);
        sizeOfrightVec = 0;
        
    }
    
    
#ifdef PZDEBUG
    std::ofstream out("geomesh.vtk"), outtxt("gmesh.txt");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out,true);// Printing the geometric mesh in vtk
    gmesh->Print(outtxt);
#endif
    
    return gmesh;
}


// Create read geometry
TPZGeoMesh * ReadGeometry()
{
    TPZGeoMesh * gmesh = new TPZGeoMesh;
    std::string grid("sample_test.msh");
    
    TPZGmshReader Geometry;
    REAL s = 1.0;
    Geometry.SetfDimensionlessL(s);
    gmesh = Geometry.GeometricGmshMesh(grid);
    const std::string name("Sample of test");
    gmesh->SetName(name);
    
    return gmesh;
}

// Create print geometry
void PrintGeometry(TPZGeoMesh * gmesh)
{
    std::stringstream text_name;
    std::stringstream vtk_name;
    text_name   << "geometry" << ".txt";
    vtk_name    << "geometry" << ".vtk";
    std::ofstream textfile(text_name.str().c_str());
    gmesh->Print(textfile);
    std::ofstream vtkfile(vtk_name.str().c_str());
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtkfile, true);
}


// Create Computational mesh
TPZCompMesh *CMesh(TPZGeoMesh *gmesh, int pOrder)
{
    
    // Create computational mesh
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder); // Set the order of approximation polynomial
    cmesh->SetDimModel(dim); // Set the dimension of the model
    cmesh->SetAllCreateFunctionsContinuous(); // Create approximation H1 space

    /** @{
     * @brief of creating material that implements the weak formulation of the model problem
     */
    
    // *************** Begin of checking the type of material ******************************************************
    // *************** TPZCPPDarcyMat ******************************************************************************
    
    
    TPZCPPDarcyMat * material = new TPZCPPDarcyMat(matid);
    
    // *************** TPZCPPDarcyWithMem **************************************************************************

//    TPZCPPDarcyWithMem * material = new TPZCPPDarcyWithMem(matid, dim);
    
    // *************** End of checking the type of material ********************************************************

    
    // Inserting material into the mesh
    cmesh->InsertMaterialObject(material);
        
    // Insert left contour condition
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    
    val2(0,0) = 0.0;
    TPZMaterial * BCond0 = material->CreateBC(material, m_matBCbott, dirichlet, val1, val2); //Creates material that implements the bottom contour condition
    cmesh->InsertMaterialObject(BCond0); // Insert material into the mesh
    
    val2(0,0) = 10.0;
    TPZMaterial * BCond1 = material->CreateBC(material, m_matBCtop, dirichlet, val1, val2);//creates material that implements the top contour condition
    cmesh->InsertMaterialObject(BCond1); // Insert material into the mesh
    
    val2(0,0) = 0.0;
    TPZMaterial * BCond2 = material->CreateBC(material, m_matBCleft, dirichlet, val1, val2); // Creates material that implements the left contour condition
    cmesh->InsertMaterialObject(BCond2); // Insert material into the mesh
    
    val2(0,0) = 0.0;
    TPZMaterial * BCond3 = material->CreateBC(material, m_matBCright, dirichlet, val1, val2); // Creates material that implements the right contour condition
    cmesh->InsertMaterialObject(BCond3); // Insert material into the mesh
    
    //  Creates computational elements that will manage the approach space of the mesh
    cmesh->AutoBuild();
    
    return cmesh;
    
}


TPZAnalysis * Analysis(TPZCompMesh * cmesh)
{
    int numthreads = 0;
    TPZAnalysis * analysis = new TPZAnalysis(cmesh,true);
    TPZSkylineStructMatrix matrix(cmesh);
    //    TPZSymetricSpStructMatrix matrix(cmesh);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    matrix.SetNumThreads(numthreads);
    analysis->SetStructuralMatrix(matrix);
    analysis->SetSolver(step);
    return analysis;
}

void PostProcess(TPZAnalysis *an)
{
    const int dim = an->Mesh()->Dimension();
    TPZStack<std::string> scalnames, vecnames;
    std::string plotfile = "DarcyModel.vtk";
    
    scalnames.Push("p");
    scalnames.Push("k");
//    vecnames.Push("Displacement");
    an->DefineGraphMesh(dim, scalnames, vecnames, plotfile);
    an->PostProcess(postProcessResolution);
}
