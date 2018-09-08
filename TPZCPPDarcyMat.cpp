//
//  TPZCPPDarcyMat.cpp
//  PZ
//
//  Created by Manouchehr on Agust 24, 2018.
//
//


#include "TPZCPPDarcyMat.h"
#include <iostream>
#include <string>
#include "pzbndcond.h"
#include "pzaxestools.h"
#include <algorithm>
#include "pzlog.h"
#include "pzfmatrix.h"
#include <stdio.h>
#include "TPZMaterial.h"



#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.TPZCPPDarcyMat"));
#endif



/** @brief default costructor */
TPZCPPDarcyMat::TPZCPPDarcyMat(): TPZMaterial()
{
    m_Dim = 0;
    m_k_0 = 0;
    m_eta = 0;
}

/** @brief costructor based on a material id */
TPZCPPDarcyMat::TPZCPPDarcyMat(int id): TPZMaterial(id)
{
    m_Dim = 2;
    m_k_0 = 1.0;
    m_eta = 1.0;
}

/** @brief copy constructor $ */
TPZCPPDarcyMat::TPZCPPDarcyMat(const TPZCPPDarcyMat& other)
{
    m_Dim  = other.m_Dim;
    m_k_0  = other.m_k_0;
    m_eta  = other.m_eta;
}


/** @brief default destructor */
TPZCPPDarcyMat::~TPZCPPDarcyMat()
{
}


/** @brief Copy assignemnt operator $ */
TPZCPPDarcyMat& TPZCPPDarcyMat::operator = (const TPZCPPDarcyMat& other)
{
    if (this != & other) // prevent self-assignment
    {
        this->m_Dim               = other.m_Dim;
    }
    return *this;
}


/** @brief of compute permeability (Kappa) */
void TPZCPPDarcyMat::Compute_Kappa(TPZMaterialData &data, REAL &kappa)
{
    kappa = m_k_0;
}


/** @brief of contribute in 2 dimensional */
void TPZCPPDarcyMat::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE>  &ek, TPZFMatrix<STATE> &ef)
{
    m_k_0 = 1.0;
    m_eta = 1.0;
    
    // Getting the space functions
    TPZFMatrix<REAL>        &phip         =   data.phi;
    TPZFMatrix<REAL>        &grad_phi_p   =   data.dphix;
    TPZFNMatrix <9,REAL>    &axes_p	      =	  data.axes;
    
    int nphi_p = phip.Rows();

    
    // Getting the solutions and derivatives
    TPZManVector<REAL,1> p = data.sol[0];
    TPZFNMatrix <6,REAL> dp = data.dsol[0];
    
    TPZFNMatrix<6,REAL> Grad_p(2,1,0.0),Grad_phi_i(2,1,0.0),Grad_phi_j(2,1,0.0);
    Grad_p(0,0) = dp(0,0)*axes_p(0,0)+dp(1,0)*axes_p(1,0);
    Grad_p(1,0) = dp(0,0)*axes_p(0,1)+dp(1,0)*axes_p(1,1);
    
    
    // Compute permeability
    REAL k = 0.0;
    Compute_Kappa(data, k);
    
    REAL c = (k/m_eta);

    // Darcy mono-phascis flow
    for (int ip = 0; ip < nphi_p; ip++)
   {
        
        Grad_phi_i(0,0) = grad_phi_p(0,ip)*axes_p(0,0)+grad_phi_p(1,ip)*axes_p(1,0);
        Grad_phi_i(1,0) = grad_phi_p(0,ip)*axes_p(0,1)+grad_phi_p(1,ip)*axes_p(1,1);
        
        REAL dot = 0.0;
        for (int i = 0;  i < m_Dim; i++)
        {
            dot += Grad_p(i,0) * Grad_phi_i(i,0);
        }
        
            ef(ip, 0)		+=  weight *  c * dot;
        
        for (int jp = 0; jp < nphi_p; jp++)
        {
            
            Grad_phi_j(0,0) = grad_phi_p(0,jp)*axes_p(0,0)+grad_phi_p(1,jp)*axes_p(1,0);
            Grad_phi_j(1,0) = grad_phi_p(0,jp)*axes_p(0,1)+grad_phi_p(1,jp)*axes_p(1,1);
            
            REAL dot = 0.0;
            for (int i = 0;  i < m_Dim; i++)
            {
                dot += Grad_phi_j(i,0) * Grad_phi_i(i,0);
            }
            
            ek(ip, jp)		+= weight * c * dot;
        }
    }
}


/** @brief of contribute of BC_2D */
void TPZCPPDarcyMat::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
{

    TPZFMatrix<REAL>  &phip = data.phi;
    TPZManVector<REAL,1> p  = data.sol[0];

    int phrp = phip.Rows();
    short in,jn;

    const REAL BIGNUMBER = TPZMaterial::gBigNumber;

    
    // Boundaries

    // Dirichlet in Pressure
    switch (bc.Type())
    {
        case 0 : // Dp
        {
            REAL v[1];
            v[0] = bc.Val2()(0,0);    //    Pressure

            //    Darcy Equation
            for(in = 0 ; in < phrp; in++)
            {
                //    Contribution for load Vector
                ef(in,0)        += BIGNUMBER*(p[0]-v[0])*phip(in,0)*weight;    // P Pressure

                for (jn = 0 ; jn < phrp; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(in,jn)        += BIGNUMBER*phip(in,0)*phip(jn,0)*weight;    // P Pressure
                }
            }
            break;
        }
            
        // Neumman in Flux
        case 1 : // Nq
        {
            REAL v[1];
            v[0] = bc.Val2()(0,0);    //    Qn

            //    Darcy Equation
            for(in = 0 ; in < phrp; in++)
            {
                //    Normal Flux on neumman boundary
                ef(in,0)    += 1.0 *  weight * v[0] * phip(in,0);    // Qnormal
            }
            break;
        }

        default:
        {
            DebugStop();
        }
            break;
    }
}


/** Returns the Fill Data Requirement */
void TPZCPPDarcyMat::FillDataRequirements(TPZMaterialData &data)

{
    data.SetAllRequirements(false);
}


/** Returns the Fill Boundary Condition Data Requirement */
void TPZCPPDarcyMat::FillBoundaryConditionDataRequirement(int type,TPZMaterialData &data)
{
    data.SetAllRequirements(false);
}


/** Returns the print */
void TPZCPPDarcyMat::Print(std::ostream &out)
{
    out << "Material Name : "               << Name()  << "\n";
    out << "Properties for TPZCPPDarcyMat: \n";
    TPZMaterial::Print(out);
    out << "\n";
}


/** Returns the variable index associated with the name */
int TPZCPPDarcyMat::VariableIndex(const std::string &name)
{
    //	Diffusion Variables
    if(!strcmp("p",name.c_str()))				return	0;
    if(!strcmp("k",name.c_str()))				return	1;
//    if(!strcmp("vx",name.c_str()))				return	2;
//    if(!strcmp("vy",name.c_str()))				return	3;
    
    return TPZMaterial::VariableIndex(name);
}


/** Returns the number of solution variables */
int TPZCPPDarcyMat::NSolutionVariables(int var)
{
    if(var == 0)	return 1;
    if(var == 1)	return 1;
//    if(var == 2)	return 1;
//    if(var == 3)	return 1;

    return TPZMaterial::NSolutionVariables(var);
}


//	Calculate Secondary variables based on ux, uy, Pore pressure and their derivatives
void TPZCPPDarcyMat::Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout)
{
    Solout.Resize( this->NSolutionVariables(var));
    
    m_k_0 = 1.0;
    m_eta = 1.0;
    
    REAL to_Mpa     = 1; // 1.0e-6;
    REAL to_Darcy   = 1; // 1.01327e+12;
    
    // Getting the solutions and derivatives
    TPZManVector<REAL,1> p  = data.sol[0];
    TPZFNMatrix <9,REAL> dp = data.dsol[0];
    p.Print(std::cout);
    dp.Print(std::cout);

    
    // Getting the space functions
    TPZFNMatrix <9,REAL>	&axes_p	=	data.axes;
    
    // Computing Gradient of the Solution
    TPZFNMatrix<3,REAL> Grad_p(3,1,0.0);
    
    // ************************************** The value of parameters ************************
    
    // ************************	Darcy Variables ************************
    //	Pore Pressure
    if(var == 0)
    {
        Solout[0] = p[0]*to_Mpa;
        return;
    }
    
    //	Permeability
    if(var == 1)
    {
        REAL k = 0.0;
        Compute_Kappa(data, k);
        Solout[0] = k*to_Darcy;
        return;
    }
    
    
//    //	Darcy's velocity in x direction
//    if(var == 2)
//    {
//        
//        REAL k = 0.0;
//        Compute_Kappa(data, k);
//        
//        Solout[0] = -(k/m_eta) * (dp(0,0)*axes_p(0,0)+dp(1,0)*axes_p(1,0));
//        return;
//    }
//    
//    //	Darcy's velocity in y direction
//    if(var == 3)
//    {
//        
//        REAL k = 0.0;
//        Compute_Kappa(data, k);
//        
//        Solout[0] = -(k/m_eta) * (dp(0,0)*axes_p(0,1)+dp(1,0)*axes_p(1,1));
//        return;
//    }
    
}

/** @brief Unique identifier for serialization purposes */
int TPZCPPDarcyMat::ClassId() const
{
    return Hash("TPZCPPDarcyMat") ^ TPZMaterial::ClassId() << 1;
}


////////////////////////////////////////////////////////////////////
void TPZCPPDarcyMat::Write(TPZStream &buf, int withclassid) const
{
    TPZMaterial::Write(buf, withclassid);
}

////////////////////////////////////////////////////////////////////
void TPZCPPDarcyMat::Read(TPZStream &buf, void *context)
{
    TPZMaterial::Read(buf, context);
}


