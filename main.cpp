
#include <cmath>

#include <iostream>
#include <fstream>
#include <string>
#include "pzgmesh.h"
#include "TPZVTKGeoMesh.h"
#include "pzanalysis.h"
#include "pzbndcond.h"
#include "TPZCPPDarcyMat.h"


#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"


using namespace std;

//------------------Model------------------------

const int dim       = 2; // Dimension of the problem
int      uNDiv      = 3;
int      vNDiv      = 4;
int      nel        = uNDiv*vNDiv; // Number of element
int     pOrder      = 4; //Polynomial order

int     matid       = 1;   // Define id for a material (weak formulation)
int      bc0        = -1;  // define id for a material (cont contour left)
int      bc1        = -2;  // define id for a material (cont contour right)
int      bc2        = -3;  // define id for a material (cont contour bottom)
int      bc3        = -4;  // define id for a material (cont contour upper)

const int dirichlet = 0;
const int neumann   = 1;

int postProcessResolution = 0; // Define post processing resolution


// @brief Function to create the geometric mesh
TPZGeoMesh *CreateGMesh(int64_t nel, int uNDiv, int vNDiv);

// @brief Function to create  the computational mesh
TPZCompMesh *CMesh(TPZGeoMesh *gmesh, int pOrder);




// Main function of the program:
int main(int argc, char *argv[])
{

    
    TPZGeoMesh *gmesh = CreateGMesh(nel, uNDiv, vNDiv); // Function to create geometry
    
    TPZCompMesh *cmesh = CMesh(gmesh, pOrder); // Function to create polynomial mesh

    // Solving the System
    int numthreads = 0;
    bool optimizeBandwidth = false;
    TPZAnalysis analysis(cmesh, optimizeBandwidth); // Creates object of analysis that manage the analysis of  problem
    
    
    TPZSkylineStructMatrix struct_mat(cmesh);
    struct_mat.SetNumThreads(numthreads);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    analysis.SetSolver(step);
    analysis.SetStructuralMatrix(struct_mat);

    
    std::cout << "Assemble matrix with NDoF = " << cmesh->NEquations() << std::endl;
    analysis.Assemble(); // Assembla the global matrix
    
#ifdef PZDEBUG
    std::ofstream filestiff("stiffness.txt");
    analysis.Solver().Matrix()->Print("K = ",filestiff,EMathematicaInput);
    
    std::ofstream filerhs("rhs.txt");
    analysis.Rhs().Print("R = ",filerhs,EMathematicaInput);
#endif

    std::cout << "Solving Matrix " << std::endl;
    analysis.Solve();
    
#ifdef PZDEBUG
        std::ofstream file("file.txt");
        analysis.Solution().Print("sol=",file,EMathematicaInput);
#endif
    
//    analysis.Run(); // Assembles the global stiffness matrix (and the load vector) and inverts the system of equations
//    TPZFMatrix<STATE> solution = cmesh->Solution(); // Taking the solution vector
//    solution.Print("Solution",cout,EMathematicaInput); // Print the solution in Mathematica format
    
    // Post processing for paraview   
    std::cout << " Post Processing " << std::endl;
    std::string plotfile("DarcyModel.vtk");
    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("p");
    scalnames.Push("k");
    vecnames.Push("v");
    
//    int dim = gmesh->Dimension();
    analysis.DefineGraphMesh(dim, scalnames, vecnames, plotfile); // Define graphic mesh
    analysis.PostProcess(postProcessResolution); // Perform post processing
    
    std::cout << "FINISHED!" << std::endl;
    
    return 0;
}


// Create geometry
TPZGeoMesh *CreateGMesh(int64_t nel, int uNDiv, int vNDiv)
{
    TPZGeoMesh * gmesh = new TPZGeoMesh; // Initialize object of class TPZGeoMesh
    
    int64_t nnodes = (uNDiv+1)*(vNDiv+1); // number of nodes
    gmesh->NodeVec().Resize(nnodes); // Resize the size of the vector of the geometric mesh

    
    // Putting us in the loop
    for (int64_t i = 0 ; i < nnodes; i++) 
    {
        TPZManVector <REAL,3> coord(3,0.);
        coord[0]=8*((float)(i/(uNDiv+1))/vNDiv)+((float)(i%(uNDiv+1))/uNDiv)*16-8;
        coord[0]/=2;
        coord[1]=coord[0]+8-((float)(i%(uNDiv+1))/uNDiv)*16;
        gmesh->NodeVec()[i].SetCoord(coord); // Set coordinate in the vector of node
        gmesh->NodeVec()[i].SetNodeId(i);    // Assign ID to one node
    }
    
    // Creating Elements
    TPZManVector<int64_t,4> topolQuad(4,0.); // Vector that will be initialized with the index of the nodes of a quadrilateral element
    TPZManVector <int64_t,4> topolLine(2,0.); // Vector that will be initialized with the index of the nodes of a one-dimensional element
    TPZVec <int64_t> TopolPoint(1); // vector that will be initialized with the index of the no of a zero-dimensional element
    int64_t id; // id of the element that will be filled by the CreateGeoElement method
    int64_t column, row;
    
    for (int64_t iel = 0; iel < nel; iel++) 
    {
        column=iel%(uNDiv);
        row=iel/(uNDiv);
        topolQuad[0]=(row)*(uNDiv+1)+column;    // Lower left vertex
        topolQuad[1]=(row)*(uNDiv+1)+column+1;  // Lower right vertex
        topolQuad[2]=(row+1)*(uNDiv+1)+column+1;// Upper right vertex
        topolQuad[3]=(row+1)*(uNDiv+1)+column;  // Upper left vertex
        
        //column==0 <=> left side
        // Cond Left contour
        if(column==0)
        {
            topolLine[0] = topolQuad[0];
            topolLine[1] = topolQuad[3];
            gmesh->CreateGeoElement(EOned, topolLine, bc0, id);
        }
        //column==uNDiv-1 <==> right side
        // Cond Right contour
        if(column==uNDiv-1)
        {
            topolLine[0] = topolQuad[1];
            topolLine[1] = topolQuad[2];
            gmesh->CreateGeoElement(EOned, topolLine, bc1, id);
        }
        //row==0<==>lower plate(V=0)
        // Cond Bottom contour
        if(row==0)
        {
            topolLine[0] = topolQuad[0];
            topolLine[1] = topolQuad[1];
            gmesh->CreateGeoElement(EOned, topolLine, bc2, id);
        }
        //row==vNDiv-1 <==> upper plate (V=1)
        // Cond Upper Contour
        if(row==vNDiv-1)
        {
            topolLine[0] = topolQuad[2];
            topolLine[1] = topolQuad[3];
            gmesh->CreateGeoElement(EOned, topolLine, bc3, id);
        }
        gmesh->CreateGeoElement(EQuadrilateral, topolQuad, matid, id);// Creates quadrilateral element
        gmesh->ElementVec()[id];
    }
    gmesh->BuildConnectivity(); // Builds mesh neighbor connectivity
    
#ifdef PZDEBUG
    std::ofstream out("geomesh.vtk"), outtxt("gmesh.txt");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out,true);// Printing the geometric mesh in vtk
    gmesh->Print(outtxt);
#endif
    
    return gmesh;
}

// Create Computational mesh
TPZCompMesh *CMesh(TPZGeoMesh *gmesh, int pOrder)
{
    
    // Create computational mesh
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder); // Set the order of approximation polynomial
    cmesh->SetDimModel(dim); // Set the dimension of the model
    cmesh->SetAllCreateFunctionsContinuous(); // Create approximation H1 space

    
    // Creating material that implements the weak formulation of the model problem
    TPZCPPDarcyMat * material = new TPZCPPDarcyMat(matid);
    
    // Inserting material into the mesh
    cmesh->InsertMaterialObject(material);
        
    // Insert left contour condition
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZMaterial * BCond0 = material->CreateBC(material, bc0, neumann, val1, val2); // Creates material that implements the left contour condition
    cmesh->InsertMaterialObject(BCond0); // Insert material into the mesh
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZMaterial * BCond1 = material->CreateBC(material, bc1, neumann, val1, val2); // Creates material that implements the right contour condition
    cmesh->InsertMaterialObject(BCond1); // Insert material into the mesh
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZMaterial * BCond2 = material->CreateBC(material, bc2, dirichlet, val1, val2); //Creates material that implements the bottom contour condition
    cmesh->InsertMaterialObject(BCond2); // Insert material into the mesh
    
    val2(0,0) = 0.0;
    val2(1,0) = 1.0;
    TPZMaterial * BCond3 = material->CreateBC(material, bc3, dirichlet, val1, val2);//creates material that implements the top contour condition
    cmesh->InsertMaterialObject(BCond3); // Insert material into the mesh
    
    //  Creates computational elements that will manage the approach space of the mesh
    cmesh->AutoBuild();
    
    return cmesh;
    
}
