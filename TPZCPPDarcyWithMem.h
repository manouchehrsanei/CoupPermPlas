//
//  TPZCPPDarcyWithMem.cpp
//  PZ
//
//  Created by Manouchehr on Agust 24, 2018.
//
//

#ifndef TPZCPPDarcyWithMem_h
#define TPZCPPDarcyWithMem_h

#include <stdio.h>
#include "TPZMaterial.h"
#include "TPZMatWithMem.h"
#include "pzbndcond.h"
#include "pzvec.h"
#include <iostream>

#include "TPZCPPDarcyMem.h"



class TPZCPPDarcyWithMem: public TPZMatWithMem<TPZCPPDarcyMem, TPZDiscontinuousGalerkin >
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
    TPZCPPDarcyWithMem();
    
    // @brief a Constructor
    TPZCPPDarcyWithMem(int matid, int dim);
    
    // @brief a Destructor
    ~TPZCPPDarcyWithMem();
    
    /** @brief Copy constructor $ */
    TPZCPPDarcyWithMem(const TPZCPPDarcyWithMem& other);
    
    /** @brief Copy assignemnt operator $ */
    TPZCPPDarcyWithMem & operator = (const TPZCPPDarcyWithMem& other);
    
    
    void Print(std::ostream & out);
    
    std::string Name() { return "TPZCPPDarcyWithMem"; }
    
    /** @brief Returns the number of state variables associated with the material */
    int NStateVariables() {return 1;}
    
    /** @brief compute permeability (Kappa) */
    virtual void Compute_Kappa(TPZMaterialData &data, TPZFNMatrix<9,STATE> &k);
    
    void FillDataRequirements(TPZVec<TPZMaterialData > &datavec);
    
    void FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec);
    
    /** @brief It computes a contribution to the stiffness matrix and load vector at one integration point to simulation. */
    void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    
    int VariableIndex(const std::string &name);
    
    int NSolutionVariables(int var);
    
    virtual void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout);
    
    
    void Solution(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleftvec, TPZVec<TPZMaterialData> &datarightvec, int var, TPZVec<STATE> &Solout, TPZCompEl * Left, TPZCompEl * Right)
    {
        DebugStop();
    }
    
    void ContributeInterface(TPZVec<TPZMaterialData> &datavec, TPZVec<TPZMaterialData> &dataleftvec, TPZVec<TPZMaterialData> &datarightvec,
                             REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
    {
        DebugStop();
    }
    void ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleftvec, TPZVec<TPZMaterialData> &datarightvec,
                             REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
    {
        DebugStop();
    }
    
    void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ef)
    {
        DebugStop();
    }
    
    void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft,
                               REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
    {
        DebugStop();
    }
    
    void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
    {
    }
    
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
    {
    }

    /** @{
     * @name Save and Load methods
     */
    
    /** @brief Unique identifier for serialization purposes */
public:
    virtual int ClassId() const;
    
    // Save the element data to a stream
    virtual void Write(TPZStream &buf, int withclassid) const;
    
    // Read the element data from a stream
    virtual void Read(TPZStream &buf, void *context);
    
public:
    
    /** @brief dimension of the model: */
    void SetDimension(int dimension)
    {
        m_Dim = dimension;
    }
    
    int Dimension() const {return m_Dim;}

    
};


#endif /* TPZCPPDarcyWithMem_hpp */
