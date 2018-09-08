//
//  TPZCPPDarcyMem.cpp
//  PZ
//
//  Created by Manouchehr on Agust 23, 2018.
//
//

#ifndef TPZCPPDarcyMem_h
#define TPZCPPDarcyMem_h

#include <stdio.h>
#include "pzreal.h"
#include "pzfmatrix.h"

////////////////////////////////////////////////////////////
//  TPZCPPDarcyMem
////////////////////////////////////////////////////////////

class TPZCPPDarcyMem 
{
    
    /// @brief of Pore Pressure
    STATE m_PorePressure;
    
    /// @brief of Gradient of Pore Pressure
    TPZFMatrix<REAL> m_GradPorePressure;
    
    /// @brief of Absolute Permeability
    TPZFNMatrix<9,REAL>  m_kappa;
    
    
public:
    
    /// @brief Default constructor
    TPZCPPDarcyMem();
    
    /// @brief Copy constructor
    TPZCPPDarcyMem(const TPZCPPDarcyMem & other);
    
    /// @brief Assignement constructor
    const TPZCPPDarcyMem & operator=(const TPZCPPDarcyMem & other);
    
    /// @brief Desconstructor
    virtual ~TPZCPPDarcyMem();
    
    /// @brief Class name
    const std::string Name()const;
    
    /// @brief Write class attributes
    virtual void Write(TPZStream &buf, int withclassid) const;
    
    /// @brief Read class attributes
    void Read(TPZStream &buf, void *context);
    
    /// @brief Print class attributes
    virtual void Print(std::ostream &out = std::cout)const;
    
    /// @brief Print class attributes
    friend std::ostream& operator<<( std::ostream& out, const TPZCPPDarcyMem & s )
    {
        s.Print(out);
        return out;
    }

    
    virtual int ClassId() const;
    
    
public:

    
    /** @brief Set Pore Pressure */
    void Set_porePressure(STATE & pore_pressure)
    {
        m_PorePressure = pore_pressure;
    }
    
    /** @brief Get Pore Pressure */
    STATE & pore_pressure()
    {
        return m_PorePressure;
    }
    
    
    /** @brief Set Gradient of Pore Pressure */
    void Set_gradporepressure(TPZFMatrix<REAL> & grad_pore_pressure)
    {
        m_GradPorePressure = grad_pore_pressure;
    }
    
    /** @brief Get Gradient of Pore Pressure */
    TPZFMatrix<REAL> grad_pore_pressure()
    {
        return m_GradPorePressure;
    }
    
    
    /** @brief Set Absolute Permeability */
    void Set_kappa_n(TPZFMatrix<REAL> & kappa)
    {
        m_kappa = kappa;
    }
    
    /** @brief Get Absolute Permeability */
    TPZFMatrix<REAL> kappa()
    {
        return m_kappa;
    }
    
    
};

#endif /* TPZCPPDarcyMem_h */
