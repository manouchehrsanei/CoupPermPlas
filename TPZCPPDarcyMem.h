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
    STATE m_PorePressure_n;
    
    /// @brief of Gradient of Pore Pressure
    TPZFMatrix<REAL> m_GradPorePressure_n;
    
    /// @brief of Absolute Permeability
    TPZFNMatrix<9,REAL>  m_kappa_n;
    
    
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

    
    /** @brief Set Pore Pressure at n (last) state */
    void Set_porePressure_n(STATE & pore_pressure_n)
    {
        m_PorePressure_n = pore_pressure_n;
    }
    
    /** @brief Get Pore Pressure at n (last) state */
    STATE & pore_pressure_n()
    {
        return m_PorePressure_n;
    }
    
    
    /** @brief Set Gradient of Pore Pressure at n (last) state */
    void Set_gradporepressure_n(TPZFMatrix<REAL> & grad_pore_pressure_n)
    {
        m_GradPorePressure_n = grad_pore_pressure_n;
    }
    
    /** @brief Get Gradient of Pore Pressure at n (last) state */
    TPZFMatrix<REAL> grad_pore_pressure_n()
    {
        return m_GradPorePressure_n;
    }
    
    
    /** @brief Set Absolute Permeability at n (last) state */
    void Set_kappa_n(TPZFMatrix<REAL> & kappa_n)
    {
        m_kappa_n = kappa_n;
    }
    
    /** @brief Get Absolute Permeability at n (last) state */
    TPZFMatrix<REAL> kappa_n()
    {
        return m_kappa_n;
    }
    
    
    
    
};

#endif /* TPZCPPDarcyMem_h */
