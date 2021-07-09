#include"hankel_trafo.hpp"

#include<cmath>
#include<cassert>
#include<vector>
#include <fstream> //TODO remove

//This is a program to run an example of the spherical Bessel transformation also known as the hankel transformation.
//The original version is NumSBT (aanz_v2), which was written in Fortran90 and published by Peter Koval and J.D. Talman in Computer Physics Communications in 2010. (CPC license)
	// v2: https://data.mendeley.com/datasets/y294ttxyw4/1
	// v3: https://data.mendeley.com/datasets/m3fc83rytv/1
	// (corresponding article: https://www.sciencedirect.com/science/article/pii/S0010465508003329 ; here the link to the program is broken)

//The multithreading possibilities implemented in the original program are not incorporated here.
//there might be reason to implement them later
//on a very basic machine using multithreding lead to the time required for the transform to be decreased by about half.
//depending on the parameters used this might or might not make a significant difference in the overall run-time of goscalc.

//samples should preferably be a power of 2
//currently an extrapolation to a lattice double the original size is used.
	//This might not be necessary as the wave functions in the waveup-file approach 0 fast enough.
//TODO change unsigned ints into size_t (better style) (here and in hankel_trafo.hpp)

//ff is the input function f(r)
//rr is the input mesh of f(r)
//gg is the output function g(k)
//kk is the ouput mesh; kk[0] must be given as start value; Delta kk = Delta rr
//li is the l-value of the Bessel-function
//np is a correction factior, standard use is np=0
//samples is the r mesh size, samples2 = 2*samples (this was renamed from nr and nr2)


namespace constants{
  constexpr int fftw_estimate {64};
  constexpr int fftw_backward {1};
  constexpr bool with_sqrt_pi_2 {false}; //if true, result g is scaled by 1/sqrt(pi/2)
  constexpr int direction {1}; //forward hankel transform. not tested for backwards
  constexpr int np {0}; //correctional factor for use with the extrapolation.
                        //generally np=0, but can be chosen differently for
                        //different potentials.
}

namespace{

  void sbt_mesh(const unsigned int samples,
      std::vector<double> &kk, const double rhomin, const double rhomax,
      const double kmax){
      //Initialize the momentum space mesh

      double dr, cf, kpmin;
      dr = (rhomax-rhomin)/(samples-1);
      kpmin = std::log(kmax)-rhomax+rhomin;
        //set up the r and k space meshes
      kk[0] = std::exp(kpmin);
      cf = std::exp(dr);
      for(unsigned int i=1; i<samples; ++i){
        kk[i] = cf*kk[i-1];
      }
  }

  //Spherical bessel functions
  //computes a table of j_l(x) for fixed xx, Eq. (39)
  std::vector<double> XJL(const double xx, const unsigned int lc){
    //rewritten to have xj as a return value
    std::vector<double> xj(lc+1); //the return value
    double aam, aa, bbm, bb, sa, sb, qqm, aap, bbp, qq, cc;

    for(unsigned int i=0; i<=lc; ++i){
      xj[i] = 0.0;
    }
    if(fabs(xx) < 1.0E-10){
      xj[0]=1.0;
      return xj;
    }

    if(xx < 0.75*lc){
      aam = 1.0;
      aa = (2.*lc+1.)/xx;
      bbm = 0.0;
      bb = 1.0;
      sa = -1.0;
      qqm = 1.0e10;

      for(unsigned int k=0; k<51; k++){
        sb = (2*(lc+k)+1)/xx;
        aap = sb*aa+sa*aam;
        bbp = sb*bb+sa*bbm;
        aam=aa;
        bbm=bb;
        aa=aap;
        bb=bbp;
        qq=aa/bb;
        if(fabs(qq-qqm)<1.0E-15) break;
        qqm=qq;
      }
      xj[lc] = 1.0;
      if(lc>0) xj[lc-1] = qq;
      if(lc>1){
        for(unsigned int l=lc-1; l>0; --l){
          xj[l-1] = (2*l+1)*xj[l]/xx - xj[l+1];
        }
      }
      cc = (sin(xx)/xx)/xj[0];
      for(unsigned int l=0; l<=lc; ++l) xj[l] = cc*xj[l];
    }
    else{
      xj[0] = sin(xx)/xx;
      if(lc>0) xj[1] = xj[0]/xx - cos(xx)/xx;
      if(lc>1){
        for(unsigned int l=1; l<lc; ++l) xj[l+1] = (2*l+1)*xj[l]/xx - xj[l-1];
      }
    }

    return xj;
  }
}

//constructor:
Hankel_trafo::Hankel_trafo(const std::vector<double>& rr, const unsigned int l_max_con, const double k_max_con) :
        l_max(l_max_con>2 ? l_max_con : 2), //l_max needs to be at least 2 to prevent crashes
        rr_fixed(rr),
        samples(static_cast<unsigned int>(rr.size())), //TODO casts size_t to unsigned int. => change samples, samples2 into size_t
        samples2(2*samples),
        rhomin_fixed(std::log(rr_fixed[0])),
        rhomax_fixed(std::log(rr_fixed[samples-1])) {

      kk = std::vector<double>(samples);
      sbt_mesh(samples, kk, rhomin_fixed, rhomax_fixed, k_max_con);


      //sbt_plan
        //internal:
        std::vector<std::vector<double>> j_ltable;
        std::vector<double> xj;
        double factor, phi, phi1, phi2, phi3, rad, tt, dt, kappamin, dr,
          rmin, kmin, rhomin;
        const std::complex<double> ci (0.0, 1.0);

        if(samples<2){
          throw hankel_trafo_exception("samples in sbt_plan are <=1.");
        }


        //http://fftw.org/fftw3_doc/One_002dDimensional-DFTs-of-Real-Data.html#One_002dDimensional-DFTs-of-Real-Data
        //demands a length of 2*(n+1) for the real arrays with the first n
        //being padding.
        // here: n = samples*2 = samples2
        // r2c_in is n=samples2 long. according to the link above,
          //r2c_out has to be n/2+1 long. This is in complex numbers making the
          //necessary allocated space equal to samples2+2
        //sizeof(double)*2 gives allocated space for complex arrays
        temp1 = reinterpret_cast<std::complex<double>*>(
                  fftw_malloc(samples2 * sizeof(double)*2));
        temp2 = reinterpret_cast<std::complex<double>*>(
                  fftw_malloc(samples2 * sizeof(double)*2));
        r2c_in = reinterpret_cast<double*>(                 //real
                  fftw_malloc((samples2)*sizeof(double)));
        r2c_out = reinterpret_cast<std::complex<double>*>(  //complex
                    fftw_malloc((samples2 + 2)*sizeof(double)));
        c2r_in = reinterpret_cast<std::complex<double>*>(   //complex
                    fftw_malloc((samples2 + 2)*sizeof(double)));
        c2r_out = reinterpret_cast<double*>(                //real
                    fftw_malloc(samples2*sizeof(double)));

        //check for NULL pointers:
        if(!temp1 || !temp2 || !r2c_in || !r2c_out || !c2r_in || !c2r_out){
          throw hankel_trafo_exception("err: sbt_plan: fftw_malloc caused"
            "a NULL pointer.");
        }

        plan12 = fftw_plan_dft_1d(samples2,
            reinterpret_cast<fftw_complex*>(temp1),
            reinterpret_cast<fftw_complex*>(temp2),
            constants::fftw_backward, constants::fftw_estimate);
        plan_r2c = fftw_plan_dft_r2c_1d(samples2,
            r2c_in, reinterpret_cast<fftw_complex*>(r2c_out),
            constants::fftw_estimate);
        plan_c2r = fftw_plan_dft_c2r_1d(samples2,
            reinterpret_cast<fftw_complex*>(c2r_in), c2r_out,
            constants::fftw_estimate);


        sbt_rr3 = std::vector<double>(samples);
        sbt_kk3 = std::vector<double>(samples);
        premult = std::vector<double>(2*samples);
        smallr = std::vector<double>(samples);
        postdiv = std::vector<double>(samples);
        mult_table1 = std::vector<std::vector<std::complex<double>>>
            (samples, std::vector<std::complex<double>>(l_max+1));
        mult_table2 = std::vector<std::vector<std::complex<double>>>
            (samples+1, std::vector<std::complex<double>>(l_max+1));

        for(unsigned int i=0; i<samples; ++i){
          sbt_rr3[i] = std::pow(rr_fixed[i], 3);
          sbt_kk3[i] = std::pow(kk[i], 3);
        }

        rmin = rr_fixed[0];
        kmin = kk[0];
        rhomin = std::log(rmin);
        kappamin = std::log(kmin);

        //originally: dr = std::log(rr_fixed[1]/rr_fixed[0]);
          //this lead to numerical problems with samples=2^13
        dr = std::log(rr_fixed[samples-1]/rr_fixed[0]) / (samples-1);
        dt = 2.0*M_PI/(samples2*dr);

        sbt_rmin = rmin;
        sbt_kmin = kmin;
        sbt_rhomin = rhomin;
        sbt_kappamin = kappamin;

        //check FFTW
        //fourier transform of delta peak at x=0 should be vertical line at y=1.
        //therefore the sum over y at samples*2 data points should be equal to
        //the number of data points.
        //uncomment if you want to perform the test
        /*
        for(unsigned int i=0; i<samples2; ++i){
          temp1[i] = 0.0;
        }
        temp1[0] = 1.0;
        fftw_execute(plan12);

        double xx = 0;
        for(unsigned int i=0; i<samples2; ++i){
          xx += std::real(temp2[i]);
        }
        if(fabs(static_cast<double>(samples2)-xx) > 1.0E-10){
          throw hankel_trafo_exception("err: sbt_plan: problem with fftw:"
            "sum(real(temp2)) doesn't pass test for transform of delta peak.");
        }
        */


        // Obtain the r values for the extended mesh, and the values r_i^1.5
        //  in the arrays smallr and premult
        factor = std::exp(dr);
        for(size_t i=0; i<samples; ++i){
          smallr[i] = rr_fixed[0] * pow(factor , static_cast<double>(static_cast<int>(i) - static_cast<int>(samples)));
          //smallr[i]=rr_f0/fac^(samples-i)
        }
        /* same as [A]
        // old iterative version, equivalent:
        smallr[samples-1] = rr_fixed[0]/factor;
        for(unsigned int i=samples-2; i+1 > 0; --i){
        //just a circuituos way of counting to 0 with unsigned ints without generating warnings
          smallr[i]=smallr[i+1]/factor;
        }
        */
        // r_i^1.5 :
        factor = std::exp(1.5*dr);
        for(size_t i=0; i<samples2; ++i){
          premult[i] = pow(factor , static_cast<double>(static_cast<int>(i) - static_cast<int>(samples)));
        }
        //fills premult from samples to 2*samples-1 and then from samples-1 to 0.
        // [A] This was done iteratively in fortran, was changed to more intuitive version,
            //without noticable change in computation time (not extensively tested)

        /* old iterative version, equivalent in results
        premult[samples] = 1.0;
        for(unsigned int i=1; i<samples; ++i){
          premult[samples+i] = factor*premult[samples+i-1];
        }
        premult[samples-1] = 1.0/factor;
        for(unsigned int i=1; i<samples; i++){
          premult[samples-i-1] = premult[samples-i]/factor;
        }
        */

        //Obtain the values 1/k_i^1.5 in the array postdivide
        postdiv[0] = (constants::with_sqrt_pi_2 ? 1.0/sqrt(M_PI/2) : 1.0);

        for(unsigned int i=1; i<samples; ++i){
          postdiv[i]= postdiv[i-1]/factor;
        }


        //construct the array of M_l(t) times the phase
        for(unsigned int it=1; it<=samples; ++it){
          tt = (it-1)*dt; //as specified in talman2009 chapter 3: 0<=n<=N-1
          phi3 = (kappamin+rhomin)*tt; //see talman2009 Eq. (33)
          rad = sqrt(std::pow(10.5, 2) + std::pow(tt,2));
          phi = atan((2.0*tt)/21.0);
          phi1 = -10.0*phi-std::log(rad)*tt+tt+sin(phi)/(12.0*rad)
            -sin(3.0*phi)/(360.0*std::pow(rad, 3)) + sin(5.0*phi)/(1260.0*std::pow(rad, 5))
            -sin(7.0*phi)/(1680.0*std::pow(rad, 7));
          for(unsigned int ix=1; ix<=10; ++ix){ //ix corresponds to p (?)
            //phi = 2*tt/(2.0*ix-1); //probably an error in the fortran-code as this
              //would be completely inconsequential
            phi1 = phi1+atan((2.0*tt)/(2.0*ix-1)); //see talman2009 Eqs. (27, 28)
          }
          phi2 = -atan(tanh(M_PI*tt/2)); //see talman2009 Eq. (20)
          //phi2 = -atan(sinh(M_PI*tt/2)/cosh(M_PI*tt/2)); //see talman2009 Eq. (20)
          //replaced with -atan(tanh(...)), because for
          //it=1750 or tt=452.3 inf is divided by inf, which gives NaN.
            //this is because sinh(900) and cosh(900) overflow the double.
          phi = phi1+phi2+phi3;
          mult_table1[it-1][0] = sqrt(M_PI/2) *
                                  std::exp(ci*phi)/std::complex<double>(samples, 0);
            //see talman2009 Eq. (18) [rather (19)]
          if(it==1){
            mult_table1[0][0] = 0.5*mult_table1[0][0];
          }
          phi = -phi2-atan(2.0*tt);
          mult_table1[it-1][1] = std::exp(2.0*ci*phi)*mult_table1[it-1][0];
              //Apply talman2009 Eq. (24):
          for(unsigned int lk=1; lk<l_max; ++lk){
            phi = -atan(2*tt/(2*lk+1));
            mult_table1[it-1][lk+1] = std::exp(2.0*ci*phi)*mult_table1[it-1][lk-1];
          }
        }

        //make the initialization for the calculation at small k values
        //  for 2N mesh values

        j_ltable = std::vector<std::vector<double>>
            (samples2, std::vector<double>(l_max+1));
        xj = std::vector<double>(l_max+1);

        //construct a table of j_l values
        double xx;
        for(unsigned int i=0; i<samples2; ++i){
          xx = std::exp(rhomin+kappamin+i*dr);
          xj = XJL(xx, l_max);
          for(unsigned int ll=0; ll<=l_max; ++ll){
            j_ltable[i][ll] = xj[ll];
          }
        }

        for(unsigned int ll=0; ll<=l_max; ++ll){
          for(unsigned int i=0; i<samples2; ++i){
            temp1[i] = j_ltable[i][ll];
          }
          fftw_execute(plan12);
          for(unsigned int i=0; i<=samples; ++i){
            mult_table2[i][ll] = conj(temp2[i])/std::complex<double>(samples2, 0);
          }
        }
        if(constants::with_sqrt_pi_2){
          for(unsigned int i=0; i<=samples; ++i){
            for(unsigned int j=0; j<=l_max; ++j){
              mult_table2[i][j] = mult_table2[i][j]/sqrt(M_PI/2);
            }
          }
        }
}


//destructor
Hankel_trafo::~Hankel_trafo(){

    //free the arrays allocated with fftw_malloc
    fftw_free(temp1);
    fftw_free(temp2);
    fftw_free(r2c_in);
    fftw_free(r2c_out);
    fftw_free(c2r_in);
    fftw_free(c2r_out);

    fftw_destroy_plan(plan12);
    fftw_destroy_plan(plan_r2c);
    fftw_destroy_plan(plan_c2r);
}



//uses the subroutines below to calculate the hankel transformation
Value_pairs Hankel_trafo::perform_hankel_trafo(const std::vector<double>& ff, const unsigned int li){
  //li is the l-value of the bessel function
  //ff is the integrand to be transformed. See talman2009 for definition (f(r)*r^2)

  std::vector<double> gg(samples);
  Value_pairs result(samples);


  //sbt_execute

    unsigned int kdiv;
    double factor, bigC, dr; //bigC renamed from C in fortran

    if(constants::direction==1){
      //ff in coordinate space (on mesh rr)
      //gg computed in momentum space (on kk)
      dr = std::log(rr_fixed[1]/rr_fixed[0]); //logarithmic step size
      bigC = ff[0]/std::pow(rr_fixed[0], constants::np+li);
    }
    else if(constants::direction==-1){
      //ff in momentum space (on mesh kk)
      //gg computed in coordinate space (on rr)
      dr = std::log(kk[1]/kk[0]);
      bigC = ff[0]/std::pow(kk[0], constants::np+li);
    }
    else{
      throw hankel_trafo_exception("hankel trafo direction must equal 1 or -1.");
    }

    // Make the calculation for LARGE k values extend the input
    // to the doubled mesh, extrapolating the input as C r^(np+li)
    //std::ofstream outFileR("test_output_r.txt"); //TODO remove; just replaces the file a bunch of times, but that's fine, I can just test for the last occurance
    for(unsigned int i=0; i<samples; ++i){
      r2c_in[i] = bigC * premult[i] * std::pow(smallr[i], constants::np + li);
    }
    for(unsigned int i=samples; i<samples2; ++i){
      r2c_in[i] = premult[i]*ff[i-samples];

    }
    for(unsigned int i=0; i<samples; ++i){
      //outFileR << ff[i] << "\n"; //TODO remove
    }
    //outFileR << smallr[i] << "\n"; //TODO remove
    
    fftw_execute(plan_r2c);

    // obtain the large k results in the vector gg
    for(unsigned int i=0; i<samples; ++i){
      temp1[i] = conj(r2c_out[i]) * mult_table1[i][li];
    }
    for(unsigned int i=samples; i<samples2; ++i){
      temp1[i] = 0.0;
    }
    fftw_execute(plan12);
    factor = std::pow(sbt_rmin/sbt_kmin, 1.5);
    for(unsigned int i=0; i<samples; ++i){
      gg[i] = factor * std::real(temp2[i+samples]) * postdiv[i];
    }

    //obtain the SMALL k results in the array c2r_out
    if(constants::direction==1){
      for(unsigned int i=0; i<samples; ++i){
        r2c_in[i] = sbt_rr3[i] * ff[i];
      }
    }
    else{	
      for(unsigned int i=0; i<samples; ++i){
        r2c_in[i] = sbt_kk3[i] * ff[i];
      }
    }

    for(unsigned int i=samples; i<samples2; ++i){
      r2c_in[i] = 0.0;
    }
    fftw_execute(plan_r2c);
    for(unsigned int i=0; i<samples+1; ++i){
      c2r_in[i] = conj(r2c_out[i]) * mult_table2[i][li];
      //outFileR << real(c2r_in[i]) << "\n"; //TODO remove
      //outFileI << imag(c2r_in[i]) << "\n"; //TODO remove
    }
    // testing start; TODO: remove:
/*
    #include <fstream>
    std::ofstream outFile("test_output.txt");
    for (const auto &test_v_elem.real : c2r_in) outFile << test_v_elem << "\n";
*/
    // testing end
    fftw_execute(plan_c2r);
    for(unsigned int i=0; i<samples; ++i){ //was 2 seperate loops in fortran
      c2r_out[i] = c2r_out[i] * dr;
      r2c_in[i] = fabs(gg[i] - c2r_out[i]);
    }
    kdiv = 0; //Position of minimal value in r2c_in 	TODO change to size_t
    for(unsigned int i=0; i<samples; ++i){
      if(r2c_in[i]<r2c_in[kdiv]){
          kdiv = i;
      }
    }
    for(unsigned int i=0; i<=kdiv; ++i){
      gg[i] = c2r_out[i];
    }

  for(unsigned int i=0; i<samples; ++i){
    result(i, 0) = kk[i];
    result(i, 1) = gg[i];
  }

  return result;
}

