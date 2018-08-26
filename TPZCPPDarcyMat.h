//
//  TPZCPPDarcyMat.cpp
//  PZ
//
//  Created by Manouchehr on Agust 24, 2018.
//
//

#ifndef TPZCPPDarcyMat_h
#define TPZCPPDarcyMat_h

#include <stdio.h>
#include "TPZMaterial.h"
#include "TPZMatWithMem.h"
#include "pzbndcond.h"
#include "pzvec.h"
#include <iostream>


class TPZCPPDarcyMat: public TPZMaterial
{
    
protected:
    
    /** @brief Problem dimension */
    int m_Dim;
    
    /** @brief Initial Permeability of the rock */
    REAL m_k_0;
    
    /** @brief Fluid viscosity */
    REAL m_eta;
    

    
public:
    
    // @brief a Defult Constructor
    TPZCPPDarcyMat();
    
    // @brief a Constructor
    TPZCPPDarcyMat(int matid);
    
    // @brief a Destructor
    ~TPZCPPDarcyMat();
    
    /** @brief Copy constructor $ */
    TPZCPPDarcyMat(const TPZCPPDarcyMat& other);
    
    /** @brief Copy assignemnt operator $ */
    TPZCPPDarcyMat & operator = (const TPZCPPDarcyMat& other);
    
    
    void Print(std::ostream & out);
    
    std::string Name() { return "TPZCPPDarcyMat"; }
    
    /** @brief compute permeability (Kappa) */
    virtual void Compute_Kappa(TPZMaterialData &data, REAL &kappa);
    
    void FillDataRequirements(TPZVec<TPZMaterialData > &datavec);
    
    void FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec);
    
    /** @brief It computes a contribution to the stiffness matrix and load vector at one integration point to simulation. */
    void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    
    int VariableIndex(const std::string &name);
    
    int NSolutionVariables(int var);
    
    void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout);

    // Save the element data to a stream
    virtual void Write(TPZStream &buf, int withclassid) const;
    
    // Read the element data from a stream
    void Read(TPZStream &buf, void *context);
    
};


#endif /* TPZCPPDarcyMat_hpp */
