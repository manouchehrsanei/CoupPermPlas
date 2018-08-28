//
//  TPZCPPDarcyMem.cpp
//  PZ
//
//  Created by Manouchehr on Agust 23, 2018.
//
//

#include "TPZCPPDarcyMem.h"
#include <stdio.h>
#include "pzreal.h"
#include "pzfmatrix.h"

TPZCPPDarcyMem::TPZCPPDarcyMem()
{
     /** @brief Pore Pressure at n (last) state */
    m_PorePressure_n    = 0;
    
    /** @brief Gradient of Pore Pressure at n (last) state */
    m_GradPorePressure_n.Resize(3, 3);
    m_GradPorePressure_n.Zero();
    
     /** @brief Permeability at n (last) state */
    m_kappa_n.Resize(3, 3);
    m_kappa_n.Zero();
}

/// Copy constructor
TPZCPPDarcyMem::TPZCPPDarcyMem(const TPZCPPDarcyMem & other)
{
    if(&other != this)
    {
        m_PorePressure_n       = other.m_PorePressure_n;
        m_GradPorePressure_n   = other.m_GradPorePressure_n;
        m_kappa_n              = other.m_kappa_n;
    }
}

/// Assignement constructor
const TPZCPPDarcyMem & TPZCPPDarcyMem::operator=(const TPZCPPDarcyMem & other)
{
    // check for self-assignment
    if(&other == this)
    {
        return *this;
    }
    m_PorePressure_n       = other.m_PorePressure_n;
    m_GradPorePressure_n   = other.m_GradPorePressure_n;
    m_kappa_n              = other.m_kappa_n;
    return *this;
}

/// Desconstructor
TPZCPPDarcyMem::~TPZCPPDarcyMem()
{
    
}

/// Class name
const std::string TPZCPPDarcyMem::Name()const
{
    return "TPZCPPDarcyMem";
}

/// Write class attributes
void TPZCPPDarcyMem::Write(TPZStream &buf, int withclassid) const
{
    buf.Write(&m_PorePressure_n);

}

/// Read class attributes
void TPZCPPDarcyMem::Read(TPZStream &buf, void *context)
{
    buf.Read(&m_PorePressure_n);

}

/// Print class attributes
void TPZCPPDarcyMem::Print(std::ostream &out) const
{
    out << Name();
    out << "\n Pore Pressure             = " << m_PorePressure_n;
    out << "\n Gradient of Pore Pressure = " << m_GradPorePressure_n;
    out << "\n Absolute Permeability     = " << m_kappa_n;
}


int TPZCPPDarcyMem::ClassId() const{
    return Hash("TPZCPPDarcyMem");
}
