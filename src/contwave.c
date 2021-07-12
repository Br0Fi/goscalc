//additions for compatibility:
#include"besselgen.c"

#include <stdexcept>


using namespace std;

/**
 *@file contwave.c
 *@brief Calculates wave function of the free atomic electron \cite Frigge2011
 */


//****************************************************************
// This function calculates the non relativistic contiuum radial
// wavefunction (Rel) on a logarhytmic mesh (rr) by solving the
// Schroedinger equation for a given potential (atompot).
// And in a second step the normalization of Rel is done
// according to Cowan, "The theory of atomic structure and
// spectra" p.519. 'Pel(r)=Pi**(-0.5)*(2e)**(-0,25)*r*Xel(r)' were
// 'r*Xel(r)' is normalized to unit amplitude at infinity
// 'alphal/cos(deltal)/q'. "e" is here the energy in Hartree.
// To calculate these amplitude the wavefunction is matsh to the
// analytic expression fitReL=alpha jL(qr)+beta nL(qr), with
// jL, (nL)=sperical Bessel, (Neuman) functions and according to
// W. Nolting "Grundkurs theoretische Physik 5/2" for Lim r--> oo
// ReL=alphaL/cos(deltaL)*sin(qr-L*PI/2+deltaL)/q/r.
//
// The file besselgen.cpp contains the calculation of the Bessel-
// and Neumannfunctions, which are needed for norming the
// wavefunctions.
//
// This code is part of a programm to calculate the energy
// differential cross-section. the main function is defined in the
// file dsigma.cpp.
//
//  Martin Frigge                                     26.06.2011
//  **************************************************************


// adams extrapolation and interpolation formulas for outward integration
// See: Abramowitz,  M.  and  I.A.  Stegun:  "Handbook  of  mathematical  functions
// with  Formulas, Graphs and mathematical tables", Dover (New York, 1965), p. 896
double aeo(vector<double> &y,int &j){
    double output;
    output=(4.16666666667e-2)*(55.*y[j]-59.*y[j-1]+37.*y[j-2]-9.*y[j-3]);
    return output;
}
double aio(vector<double> &y,int &j){
    double output;
    output=(4.16666666667e-2)*(9.*y[j+1]+19.*y[j]-5.*y[j-1]+y[j-2]);
    return output;
}

//******************  decalaration of fuctions  ********************************
double SPHFUN(int N, double X,double &neumann);

//calculation of continuum radial wavefunctions
void contwave(double amesh,int emax,int llmax, double einc,double emin,int mmax,
              vector<double> &rr, vector<double> &atompot,
              vector<vector<vector< double> > > &Rel){


    vector<double> bessel(3),neumann(3);
    double dbessel,dneumann; //derivation of bessel- and neumannfunctions
    double dR, R,r, qr,q, beta, alpha,delta;

    //double M0=9.1094e-31;              // kg -> Ruhemasse e-
    //double H_2=6.58211899E-16*1.05457266E-34; // eVs*Js=Kg*eV*m*m  Plancksches Wirkungsquantum
    double PI=3.1415926535;
    //double A0=5.2918e-11;              //  Bohrscher Radius in m
    //double Ry=13.6058;                  // Rydberg-Konstante in eV


    int fitpoint=0; //XXXchanged: fitpoint -> fitpoint = 0
    double gamma,sls,al,energy,als;//,fss
    vector<double> u(mmax),up(mmax),upp(mmax), cf(mmax), dv(mmax), fr(mmax),frp(mmax);
    al=log(amesh);
    als=al*al;


    //******************* start calcualtion ****************************
    // find the machpoint j, where V(j)=0
    for(int j=0; j<mmax; ++j){
        if(approx_equal(atompot[j],0.0, 1.0e-8)){
            fitpoint=j;
            j=mmax;
        }
    }
    if(fitpoint==0){
      //if there is no point close enough to zero found, reduce the critirium and
        //try outmost point. If it's still not close to zero, set error condition
      if(approx_equal(atompot[mmax-1],0.0,1.0e-5)){fitpoint=mmax-1;}
      else{throw std::runtime_error("atomic potential does not approach zero fast enough.");}
    }
    else{
      fitpoint=fitpoint+10;
      if(fitpoint>=mmax-1)
          fitpoint=mmax-1;
    }

    //convert energy from eV to Hartree;	1eV /(2* 13.605eV)=0.036749325 (Hartree)
    einc=einc*0.036749325;
    emin=emin*0.036749325;

    for(int l=0; l<=llmax; l++){
        for(int e=0; e<=emax; e++){
            energy=e*einc+emin;//in Hartree
            //fss=0.;
            gamma=l+1.;
            sls=l*(l+1.);
            for(int j=0; j<mmax; j++){
               cf[j]=als*sls + 2.*als*(atompot[j]-energy)*rr[j]*rr[j];
//                cf[j]=als*sls + 2.*als*(0.-energy)*rr[j]*rr[j];
            }




            //calculate dv/dr for darwin correction
            dv[0]=(-50.*atompot[0]+96.*atompot[1]-72.*atompot[2]+32.*atompot[3]-6.*atompot[4])/(24.*al*rr[0]);
            dv[1]=(-6.*atompot[0]-20.*atompot[1]+36.*atompot[2]-12.*atompot[3]+2.*atompot[4])/(24.*al*rr[1]);

            for(int j=2; j<mmax-2; j++){
                dv[j]=(2.*atompot[j-2]-16.*atompot[j-1]+16.*atompot[j+1]-2.*atompot[j+2])/(24.*al*rr[j]);
            }
            dv[mmax-1]=0.;
            dv[mmax-2]=0;

            //relativistic coefficient arrays for u (fr) and up (frp).
            for(int j=0; j<mmax; j++){
//                fr[j]=als*(rr[j]*rr[j])*(pow(-fss*(atompot[0][j]-energy),2.) + 0.5*fss*dv[j]/(rr[j]*(1.+0.5*fss*(energy-atompot[0][j]))));
//                frp[j]=-al*rr[j]*0.5*fss*dv[j]/(1.+0.5*fss*(energy-atompot[0][j]));
                fr[j]=0.;
                frp[j]=0.;
            }

            for(int j=0;j<4;j++){
                u[j]=pow(rr[j],gamma);
                up[j]=al*gamma*pow(rr[j],gamma);
                upp[j]=(al+frp[j])*up[j]+(cf[j]+fr[j])*u[j];
            }
            //outward integration using predictor once, corrector twice
            for(int j=3; j<mmax-1; j++){
                u[j+1]=u[j]+aeo(up,j);
                up[j+1]=up[j]+aeo(upp,j);
                for(int it=0; it<2; it++){
                    upp[j+1]=(al+frp[j+1])*up[j+1]+(cf[j+1]+fr[j+1])*u[j+1];
                    up[j+1]=up[j]+aio(upp,j);
                    u[j+1]=u[j]+aio(up,j);
                }
            }


            //calculate R and dR/dr at fitpoint(V(r=fitpoint)=0)
            r=rr[fitpoint];  //in au
            R=u[fitpoint]/r;
            q=sqrt(2.*energy);
            qr=q*r;
            dR=(u[fitpoint+1]/rr[fitpoint+1]-u[fitpoint-1]/rr[fitpoint-1])/
               (rr[fitpoint+1]-rr[fitpoint-1]);

            if(l==0){

                // calculate fitparameter alpha, beta and delta;
                // fitReL=alpha jL(qr)+beta nL(qr); with jL, (nL)=sperical Bessel, (Neuman) function
                // Lim r--> oo ReL=alphaL/cos(deltaL)*sin(qr-L*PI/2+deltaL)/q/r according to Nolting
                // 5/2 "Grundkurs theoretische Physik"
                beta=(dR+R/r-R*q/tan(qr))*r*sin(qr);
                alpha=beta/tan(qr)+R*qr/sin(qr);
                delta=atan(-beta/alpha);

                for(int j=0; j<mmax; j++){
                    Rel[j][e][l]=u[j]/rr[j]*
                                   (q/alpha*cos(delta)/pow(PI,0.5)*pow(2./energy,0.25));
                }
            }
            else{//for l>0

                // calculate spherical bessel- and neumannfunctions
                bessel[0]=SPHFUN(l,q*rr[fitpoint+1],neumann[0]);
                bessel[1]=SPHFUN(l, qr,neumann[1]);
                bessel[2]=SPHFUN(l, q*rr[fitpoint-1],neumann[2]);
                dbessel=(bessel[0]-bessel[2])/
                        (rr[fitpoint+1]-rr[fitpoint-1]);
                dneumann=(neumann[0]-neumann[2])/
                        (rr[fitpoint+1]-rr[fitpoint-1]);

                // calcualte fitparameter
                beta=(dR-R*dbessel/bessel[1])/(dneumann-neumann[1]*dbessel/bessel[1]);
                alpha=(R-beta*neumann[1])/bessel[1];
                delta=atan(-beta/alpha);

                for(int j=0;j<mmax;j++){
                    Rel[j][e][l]=u[j]/rr[j]*
                                   (q/alpha*cos(delta)/pow(PI,0.5)*pow(2./energy,0.25));
                }
            }//END of l>0
        }//END of e-loop
    }//END of l-loop
}//END of contwave
