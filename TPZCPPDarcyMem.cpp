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
     /** @brief Pore Pressure */
    m_PorePressure    = 0;
    
    /** @brief Gradient of Pore Pressure */
    m_GradPorePressure.Resize(3, 3);
    m_GradPorePressure.Zero();
    
     /** @brief Permeability */
    m_kappa.Resize(3, 3);
    m_kappa.Zero();
}

/// Copy constructor
TPZCPPDarcyMem::TPZCPPDarcyMem(const TPZCPPDarcyMem & other)
{
    if(&other != this)
    {
        m_PorePressure       = other.m_PorePressure;
        m_GradPorePressure   = other.m_GradPorePressure;
        m_kappa              = other.m_kappa;
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
    m_PorePressure       = other.m_PorePressure;
    m_GradPorePressure   = other.m_GradPorePressure;
    m_kappa              = other.m_kappa;
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
    buf.Write(&m_PorePressure);

}

/// Read class attributes
void TPZCPPDarcyMem::Read(TPZStream &buf, void *context)
{
    buf.Read(&m_PorePressure);

}

/// Print class attributes
void TPZCPPDarcyMem::Print(std::ostream &out) const
{
    out << Name();
    out << "\n Pore Pressure             = " << m_PorePressure;
    out << "\n Gradient of Pore Pressure = " << m_GradPorePressure;
    out << "\n Absolute Permeability     = " << m_kappa;
}


int TPZCPPDarcyMem::ClassId() const{
    return Hash("TPZCPPDarcyMem");
}
