#include <TF1.h>
#include <TRandom3.h>

#include "HepMC/GenEvent.h"
#include "HepMC/GenVertex.h"
#include "HepMC/GenParticle.h"
#include "HepMC/IO_GenEvent.h"
#include "HepMC/IO_AsciiParticles.h"

#include <iostream>

using namespace std;

double InverseBeta(double Enu,double CosThetaLab);
double neuFluxIsotope(int type, double Eneu);
double neuFluxAll(double *x,double *par);

const double XcMeVtoCmsqrd = 0.389379292e-21;
const double GFERMI = 1.16639e-11;
const double MPROTON = 938.27200;
const double MNEUTRON = 939.57;
const double MELECTRON = 0.510998902;   //likewise
const double PI = 3.1415928;    //Likewise
const double DELTA = 1.2933318;  //likewise

// average energy released per fission.
// U235  = 0 // Pu239 = 1 // U238  = 2 // Pu241 = 3
// P.Huber PRD 70, 053011 (2004)
const double isoFrac[4] = {0.6977, 0.1912, 0.0765, 0.0346};
const double EperFission[4] = {201.7,210.0,205.0,212.4};
double EneuMin = 1.80, EneuMax = 12.0;
TF1 neuSpect("f0", neuFluxAll, EneuMin, EneuMax, 4);

int main(int argc, char** argv)
{
  if ( argc < 3 ) {
    cout << "Usage: " << argv[0] << " OUTPUT_HEPMC_FILENAME.dat NUMBER_OF_EVENTS [EVENT_NUMBER_BEGIN]\n";
    return 1;
  }

  const std::string outFileName = argv[1];
  const int nEvent = atoi(argv[2]);
  const int iEvent = argc > 3 ? atoi(argv[3]) : 1;
  const int seed = 12345 + iEvent;

  //const double PI = 3.1415928;    //Likewise
  //const double eVtoJ = 1.60217649*pow(10.,-19.); // conversion factor eV to J.
  const double xMax = InverseBeta(EneuMax,-1);

  gRandom = new TRandom3(seed);
  for (int i=0; i<4; ++i ) neuSpect.SetParameter(i, isoFrac[i]);

  // Start main loop
  HepMC::IO_GenEvent out(outFileName.c_str());
  //HepMC::IO_AsciiParticles out(outFileName.c_str());
  for ( int i=0; i<nEvent; ++i ) {
    cout << "@@ Processing event " << i+iEvent << endl;

    // Get the neutrino energy from the reactor
    double energy = 0, cosTheta = 0;
    while ( true ) {
      energy = neuSpect.GetRandom();
      const double xran = gRandom->Uniform(0,1);
      cosTheta = gRandom->Uniform(-1,1);
      if ( xran < InverseBeta(energy, cosTheta)/xMax ) break;
    }

    // Neutron momenta
    const HepMC::FourVector neutronP4 = [&](){
      const double E0 = energy - DELTA;
      const double Ysquared = (DELTA*DELTA-MELECTRON*MELECTRON)/2;
      const double p0 = sqrt(E0*E0-MELECTRON*MELECTRON);
      const double v0 = p0/E0;
      const double E1 = E0*(1-energy/MPROTON*(1-v0*cosTheta))-Ysquared/MPROTON;
      const double p1 = sqrt(E1*E1-MELECTRON*MELECTRON);
      const double phiR = gRandom->Uniform(0,2*PI);

      return HepMC::FourVector(-p1*sqrt(1-cosTheta*cosTheta)*cos(phiR),
                               -p1*sqrt(1-cosTheta*cosTheta)*sin(phiR),
                               energy-p1*cosTheta, hypot(p1, MNEUTRON));
    }();

    // Electron momenta
    const HepMC::FourVector electronP4 = [&](){
      const double px = -neutronP4.px();
      const double py = -neutronP4.py();
      const double pz = 0+energy-neutronP4.pz();
      const double e = sqrt(px*px+py+py+pz*pz+MELECTRON*MELECTRON);
      return HepMC::FourVector(px, py, pz, e);
    }();

    // Fill into HepMC event
    HepMC::GenEvent* event = new HepMC::GenEvent();
    event->set_event_number(i+iEvent);

    HepMC::GenParticle* pp = new HepMC::GenParticle(HepMC::FourVector(0,0,0,MPROTON), 2212, 4);
    HepMC::GenParticle* pv = new HepMC::GenParticle(HepMC::FourVector(0,0,energy,energy), -12, 4);
    pp->set_generated_mass(MPROTON);
    pv->set_generated_mass(0);

    HepMC::GenParticle* pn = new HepMC::GenParticle(neutronP4, 2112, 1);
    HepMC::GenParticle* pe = new HepMC::GenParticle(electronP4, 11, 1);
    pn->set_generated_mass(MNEUTRON);
    pe->set_generated_mass(MELECTRON);

    HepMC::GenVertex* v0 = new HepMC::GenVertex(HepMC::FourVector(0,0,0,0)); // No vertex smearing yet, to be decided
    v0->add_particle_in(pp);
    v0->add_particle_in(pv);
    v0->add_particle_out(pn);
    v0->add_particle_out(pe);
    event->add_vertex(v0);

    event->set_beam_particles(pp, pv);
    event->set_signal_process_vertex(v0);

    out.write_event(event);
  }

  return 0;
}

double InverseBeta(double Enu,double CosThetaLab){ //
/********************************************************************************
lifted from IBDgenerator.cc 
********************************************************************************/

  // Cross section constants.  Some for overall scale are just
  // to allow absolute comparison to published articel.
  //
  double XCunits = XcMeVtoCmsqrd;
  double CosThetaC = (0.9741+0.9756)/2;
  //
  // Radiative correction constant
  //
  double RadCor = 0.024;
  //
  // check for threshold
  //
  double EminBeta =
    ((MPROTON+DELTA+MELECTRON)*(MPROTON+MELECTRON+DELTA)
     -MPROTON*MPROTON)/2/MPROTON;
  if(Enu<EminBeta)return 0;
  //
  // overall scale
  //
  double Sigma0 = GFERMI*GFERMI*CosThetaC*CosThetaC/PI*(1+RadCor);
  //
  // couplings
  //
  double f = 1.00;
  double f2 = 3.706;
  double g = 1.26;
  //
  //
  double E0 = Enu - DELTA;
  if(E0<MELECTRON)E0=MELECTRON;
  double p0 = sqrt(E0*E0-MELECTRON*MELECTRON);
  double v0 = p0/E0;
  //
  //  order 1 terms
  //
  double Ysquared = (DELTA*DELTA-MELECTRON*MELECTRON)/2;
  double E1 = E0*(1-Enu/MPROTON*(1-v0*CosThetaLab))-Ysquared/MPROTON;
  if(E1<MELECTRON)E1=MELECTRON;
  double p1 = sqrt(E1*E1-MELECTRON*MELECTRON);
  double v1 = p1/E1;

  double Gamma =
    2*(f+f2)*g*((2*E0+DELTA)*(1-v0*CosThetaLab)-
                MELECTRON*MELECTRON/E0)
    +(f*f+g*g)*(DELTA*(1+v0*CosThetaLab)+MELECTRON*MELECTRON/E0)
    +(f*f+3*g*g)*((E0+DELTA)*(1-CosThetaLab/v0)-DELTA)
    +(f*f-g*g)*((E0+DELTA)*(1-CosThetaLab/v0)-DELTA)*v0*CosThetaLab;

  double XC =
    ((f*f+3*g*g)+(f*f-g*g)*v1*CosThetaLab)*E1*p1
    -Gamma/MPROTON*E0*p0;
  XC *= Sigma0/2;
  XC *= XCunits;
  return XC;
}

double neuFluxIsotope(int type, double Eneu) {
// This class returns the neutrino flux spectrum without oscillations
// Parametrisation by P. Vogel and J. Engel PRD39, 3378 (1989)
// updated with
// P. Huber and T. Schwetz Phys. Rev. D70, 053011 (2004) which has six parameter
// parametrisation used by D-C. 
/*
  U235  = 0
  Pu239 = 1
  U238  = 2
  Pu241 = 3
*/

// modified by ydkim  2012. 4. 5
// P. Huber arxiv:1106.0687 (2011) v7 (2012) except U238
  const double a[4][6] = {
    {4.367, -4.577, 2.100, -5.294E-1, 6.186E-2, -2.777E-3},
    {4.757, -5.392, 2.563, -6.596E-1, 7.820E-2, -3.536E-3},
    {4.833E-1, 1.927E-1, -1.283E-1, -6.762E-3, 2.233E-3, -1.536E-4}, // 238U from Mueller's paper, prc83, 054615 (2011)
    {2.990, -2.882, 1.278, -3.343E-1, 3.905E-2, -1.754E-3}
  };

  double f=0;
  for (int i=0;i<6;i++) f += a[type][i]*pow(Eneu,(double)i);
  return exp(f);
}

double neuFluxAll(double *x,double *par)
{
  double f = 0;
  for (int i=0;i<4;i++) {
    double z = neuFluxIsotope(i,x[0])*par[i];
    f += z;
  }
  return f;
}


