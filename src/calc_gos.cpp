#include"calc_gos.hpp"

#include<assert.h>

#include<armadillo>
#include<wignerSymbols.h>

#include"goscalc_constants.hpp"




std::vector<double> calc_energy_losses(const Goscalc_shared_data& gsd){
  const unsigned no_energy_losses = gsd.cd.energy_free_steps+1;
  std::vector<double> energy_losses(no_energy_losses);
  for (unsigned int i_e = 0; i_e < no_energy_losses; i_e++){
    const double energy_free = gsd.cd.energy_free_start+i_e*gsd.cd.energy_free_increase;
    energy_losses[i_e] = abs(gsd.cd.bound_energy) + energy_free;
  }

  return energy_losses;
}



Value_pairs calculate_gos(unsigned int index_eloss, Hankel_trafo& trafo_plan, const Goscalc_shared_data& gsd){
  const unsigned int l = gsd.cd.l_bound;


  Value_pairs result;

  for (unsigned int lfree=0; lfree<gsd.free_waves[index_eloss].size(); lfree++){
    Value_pairs integrand = gsd.prepd_bound_wave*gsd.free_waves[index_eloss][lfree];

    //see leapman1980 eq. 4 and eq. 7
    int lower_valid_lambda = std::abs(static_cast<int>(lfree-l));
    int upper_valid_lambda = static_cast<int>(lfree+l);





    for (int lambda = lower_valid_lambda; lambda<=upper_valid_lambda; lambda++){
      double wigner_symb = WignerSymbols::wigner3j(lfree,lambda,l,0,0,0);
      if (wigner_symb == 0){continue;} //TODO comparing doubles with == is bad practice

      const Value_pairs rad_int = trafo_plan.perform_hankel_trafo(integrand.column_as_vector(1), static_cast<unsigned int>(lambda));

      if ((lambda == lower_valid_lambda) && (lfree == 0)){
        result = (2*lfree+1)*(2*lambda+1)*pow(wigner_symb,2)*rad_int*rad_int;
      } else {
        result += (2*lfree+1)*(2*lambda+1)*pow(wigner_symb,2)*rad_int*rad_int;
      }

    }
  }




  const arma::colvec inv_q2 = arma::pow(result.column_as_colvec(0),-2); //getestet. tut
  const double energy_loss = gsd.results_energy_losses[index_eloss];
  const arma::colvec prefactor = cst::m_e2dhbarh2*energy_loss*inv_q2; //unitless(!)




  return prefactor*result;


}

/*For testing purposes include 

#include<fstream> 
#include<vector>

and use the following to put out a table for a specific energy loss and free l:
to create a table:
    if(index_eloss==0&&lfree==3){
      board = std::vector<std::vector<double>>(integrand.size()+1, std::vector<double>(upper_valid_lambda - lower_valid_lambda + 1));
    }
to fill the table:
      if(index_eloss==0 && lfree==3){
        //output = rad_int.column_as_vector(1);
        //board[0][lambda-lower_valid_lambda] = lambda; //uninitialisiert (0) für wigner=0
        //for(unsigned int i=1; i<output.size()+1; ++i){
          //board[i][lambda-lower_valid_lambda] = output[i-1];
        //}
        std::ofstream ofsin ("input.dat", std::ofstream::out);
        gsd.free_waves[0][0].column_as_colvec(1).raw_print(ofsin);
        ofsin.close();
      }
to export the table to a file:
  if(index_eloss==0){
    //std::ofstream ofs ("test_sbt.txt", std::ofstream::out);
    //for(unsigned int i=0; i<board.size(); ++i){
      //for(unsigned int j=0; j<board[0].size(); ++j){
        //ofs << board[i][j];
        //if(j<board[0].size()-1){ofs << "\t";}
      //}
      //ofs << "\n";
    //}
    //ofs.close();
    //std::ofstream ofs2 ("result.txt", std::ofstream::out);
    //result.column_as_colvec(1).raw_print(ofs2);
    //ofs2.close();
    std::ofstream ofsr ("r.dat", std::ofstream::out);
    gsd.free_waves[0][0].column_as_colvec(0).raw_print(ofsr);
    ofsr.close();
  }
*/

