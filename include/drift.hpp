#ifndef drift_HPP
#define drift_HPP

#include "treeStructure.hpp"
#include "material_paramaters.hpp"
#include "random.hpp"
#include <math.h>


/*
void CheckBoundary(ContainerClass* particles, geometry geometry, FSize idxPart)
{   
}
*/

void Drift(geometryClass DevGeometry, ContainerClass* particles, mat_paramClass *mat_par, FSize idxPart, double tau)
{
    double x,y,z, fx, fy, fz, kx, ky, kz, charge, gamma, energy, invmass, m_c, ionization_energy;


	
    //Get positions from FMM //

    x = particles->getPositions()[0][idxPart];
    y = particles->getPositions()[1][idxPart];
    z = particles->getPositions()[2][idxPart];

    //Get momenta 

    kx = particles->getPhysicalValues(1,1)[idxPart];    
    ky = particles->getPhysicalValues(2,1)[idxPart];
    kz = particles->getPhysicalValues(3,1)[idxPart];

    // Get Charge

    charge =  particles->getPhysicalValues()[idxPart];
      


    // Get forces from FMM // F = q * E (Charge in FMM is set to 1)

    fx = - charge / (4*M_PI* mat_par->get_eps_0()  * mat_par->get_eps_infty()) *  mat_par->get_q() * mat_par->get_q() * particles->getForcesX()[idxPart] * DevGeometry.get_boxdim()*DevGeometry.get_boxdim();  
    fy = - charge / (4*M_PI* mat_par->get_eps_0()  * mat_par->get_eps_infty()) *  mat_par->get_q() * mat_par->get_q() * particles->getForcesY()[idxPart] * DevGeometry.get_boxdim()*DevGeometry.get_boxdim(); 
    fz = - charge / (4*M_PI* mat_par->get_eps_0()  * mat_par->get_eps_infty()) *  mat_par->get_q() * mat_par->get_q() * particles->getForcesZ()[idxPart] * DevGeometry.get_boxdim()*DevGeometry.get_boxdim();

    

    
    invmass = 1 / (mat_par->get_effmassX());


    gamma   = invmass / ( 2 * mat_par->get_q()) * mat_par->get_hbar()* mat_par->get_hbar()*(kx*kx + ky*ky + kz*kz) ; // parabolic energy, alfa -> 0

    if (mat_par->get_alpha() > 0)
        energy  =  1 / (mat_par->get_q()) * 0.5*(1 / mat_par->get_alpha())*(pow(1 + 4*gamma*mat_par->get_alpha(), 0.5) - 1);
    else
        energy = gamma ;
    
    m_c     = mat_par->get_effmassX()*pow(1 + 4*mat_par->get_alpha()*energy, 0.5); // constant used to compute displacement



    // position times dimensional factor to go back to dimensions of the box. 

    // Doping
    if (particles->getPhysicalValues()[idxPart] > 0)
        {
        // invmass /= 1e9;
        x += /*( 1 / m_c)*/ invmass * (kx * mat_par->get_hbar() * tau + charge * 0.5*fx * pow(tau,2))*DevGeometry.get_boxdim(); // kx*mat_par.hbar / m_c * tau + 0.5*fx / m_c*pow(tau,2);
        y += /*( 1 / m_c)*/ invmass * (ky * mat_par->get_hbar() * tau + charge * 0.5*fy * pow(tau,2))*DevGeometry.get_boxdim(); // ky*mat_par.hbar / m_c * tau + 0.5*fy / m_c*pow(tau,2);
        z += /*( 1 / m_c)*/ invmass * (kz * mat_par->get_hbar() * tau + charge * 0.5*fz * pow(tau,2))*DevGeometry.get_boxdim(); // kz*mat_par.hbar / m_c * tau + 0.5*fz / m_c*pow(tau,2);
         
        }

    else
    {
    x += /*( 1 / m_c)*/ invmass * (kx * mat_par->get_hbar() * tau + charge * 0.5*fx * pow(tau,2))*DevGeometry.get_boxdim(); // kx*mat_par.hbar / m_c * tau + 0.5*fx / m_c*pow(tau,2);
    y += /*( 1 / m_c)*/ invmass * (ky * mat_par->get_hbar() * tau + charge * 0.5*fy * pow(tau,2))*DevGeometry.get_boxdim(); // ky*mat_par.hbar / m_c * tau + 0.5*fy / m_c*pow(tau,2);
    z += /*( 1 / m_c)*/ invmass * (kz * mat_par->get_hbar() * tau + charge * 0.5*fz * pow(tau,2))*DevGeometry.get_boxdim(); // kz*mat_par.hbar / m_c * tau + 0.5*fz / m_c*pow(tau,2);
    }

    kx +=(charge * fx / mat_par->get_hbar()) * tau; // p = hbar * k. 1 ev -> k = 3.9e18, p = 4.14e-16
    ky +=(charge * fy / mat_par->get_hbar()) * tau;
    kz +=(charge * fz / mat_par->get_hbar()) * tau;

    // Periodic. 

    if ( x > DevGeometry.Xmax){
        x    = x - DevGeometry.Xmax;} // 2*DevGeometry.Xmax - x;}
       // kx = - kx;}
    else if (x < DevGeometry.Xmin){
        x    =  x + DevGeometry.Xmax;}
    //    kx = - kx;}

    if (y > DevGeometry.Ymax){
        y    = y - DevGeometry.Ymax;} // 2*DevGeometry.Ymax - y;}
    //    ky = -ky;}

    else if (y < DevGeometry.Ymin){
        y    =  y + DevGeometry.Ymax;}
    //    ky = - ky;}

    if (z > DevGeometry.Zmax){
        z    = z - DevGeometry.Zmax;}//2*DevGeometry.Zmax - z;}
    //    kz = - kz;}
    else if (z < DevGeometry.Zmin){
        z    =  z + DevGeometry.Zmax;}
    //    kz = - kz;}

    // set updated positions back 


    particles->getPositions()[0][idxPart] = x;
    particles->getPositions()[1][idxPart] = y;
    particles->getPositions()[2][idxPart] = z;

    // set updated momenta back 

    particles->getPhysicalValues(1,1)[idxPart] = kx;
    particles->getPhysicalValues(2,1)[idxPart] = ky;
    particles->getPhysicalValues(3,1)[idxPart] = kz;
    
  
    // Energy threshold, if Energy is above certain value remove particle

    gamma   = invmass / ( 2 * mat_par->get_q()) * mat_par->get_hbar()* mat_par->get_hbar()*(kx*kx + ky*ky + kz*kz) ; // parabolic energy, alfa -> 0

    if (mat_par->get_alpha() > 0)
        energy  =  1 / (mat_par->get_q()) * 0.5*(1 / mat_par->get_alpha())*(pow(1 + 4*gamma*mat_par->get_alpha(), 0.5) - 1); 
    else
        energy = gamma ;


    
//    ionization_energy = 4 ; 
//
//    if (energy > ionization_energy) 
//    {   
//        std::cout << "removing particle \n\n" << std::endl;
//        particles->removeParticles(&idxPart, 1);
//    }    
}   





#endif
