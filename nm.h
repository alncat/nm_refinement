#ifndef MMTBX_NM_H
#define MMTBX_NM_H

#include <cctbx/sgtbx/space_group.h>
#include <scitbx/array_family/versa.h>
#include <vector>
#include <scitbx/array_family/shared_algebra.h>
#include <scitbx/array_family/versa_algebra.h>
#include <cctbx/import_scitbx_af.h>
#include <cmath>
#include <cctbx/adptbx.h>
#include <cctbx/xray/scattering_type_registry.h>
#include <assert.h>
#include <math.h>
#include <iostream>
#include <mmtbx/error.h>
#include <cctbx/xray/targets.h>
#include <scitbx/matrix/outer_product.h>
#include <iotbx/pdb/hierarchy.h>
#include <fstream>

using namespace std;
namespace mmtbx { namespace nm {
namespace af=scitbx::af;
using scitbx::vec3;
using scitbx::mat3;
using scitbx::sym_mat3;
typedef af::shared<vec3<double> > normal_mode;

class nm_init {
private:
    af::shared<normal_mode> modes;
    af::shared<normal_mode> zero_modes;
    af::shared< std::size_t > nm_index;
    std::size_t count_zero_modes = 0;

public:
    nm_init(   
                const char* filename,
                std::size_t n_modes,
                af::shared<iotbx::pdb::hierarchy::atom> const& atoms,
                bool zero_mode_input_flag,
                bool zero_mode_flag)
    {
//      modes(atoms.size(), af::init_functor_null<vec3<normal_modes> >());
        vec3<double> null_mode_i(0.0, 0.0, 0.0);
        normal_mode null_mode;
        nm_index.resize(atoms.size(), 0);
//        for(std::size_t i = 0; i < atoms.size(); i++){
//            null_mode.push_back(null_mode_i);
//        }
//        vec3<normal_modes> null_mode;
//        for(std::size_t i = 0; i < 3; i++){
//           null_mode[i] = null_mode_i;
//       }
        std::size_t nmode_start = 0;
        if(zero_mode_flag and (not zero_mode_input_flag)){
            nmode_start = 6;
        }
        else {
            nmode_start = 0;
        }
        ifstream file(filename, ios::in | ios::binary);
        std::size_t total_evc;
        std::size_t align_mark = 0;
        if( file.is_open() ){
            for(std::size_t i = 0; i <= atoms.size(); i++){
                char * serial = new char[8];
                file.read(serial, 8);
                file.seekg(5, ios::cur);
                char * resname = new char[3];
                file.read(resname, 3);
                file.seekg(4, ios::cur);
                char * atmname = new char[4];
                file.read(atmname, 4);
                for(std::size_t j = align_mark; j < atoms.size(); j++){
                    if( atoms[j].parent()->data->resname == resname and atoms[j].data->name == atmname){
                        nm_index[i] = j;
                        align_mark = j+1;
                        break;
                    }
                }
                if(strncmp(serial, "END", 3) == 0){
                    total_evc = i;
                    break;
                    delete[] resname;
                    delete[] atmname;
                    delete[] serial;
                }
                delete[] resname;
                delete[] atmname;
                delete[] serial;
            }
            for(std::size_t i = 0; i < nmode_start; i++){
                normal_mode null_mode(atoms.size(), af::init_functor_null<vec3<double> >());
                modes.push_back(null_mode);
            }
            for(std::size_t i = nmode_start; i < n_modes; i++){
                normal_mode tmp_mode(atoms.size(), af::init_functor_null<vec3<double> >());
                for(std::size_t j = 0; j < total_evc; j++){
                    double vx, vy, vz;
                    file.read( (char *) &vx, 8);
                    file.read( (char *) &vy, 8);
                    file.read( (char *) &vz, 8);
                    std::size_t k = nm_index[j];
                    tmp_mode[k][0] = vx;
                    tmp_mode[k][1] = vy;
                    tmp_mode[k][2] = vz;
                }
//                cout << modes[i-1][0][0] << modes[i-1][0][1] << modes[i-1][0][2] <<endl;
                double vx, vy, vz;
                file.read((char *) &vx, 8);
                file.read((char *) &vy, 8);
                file.read((char *) &vz, 8);
                if( vx!=111111111. and vy!=-222222222. and vz!=333333333. ){
                    cout << "error in reading eigenvector file" <<endl;
                    cout << "please check consistency!!!!!!" << endl;
                    cout << vx << vy << vz <<endl;
                    exit(EXIT_FAILURE);
                }
                else{
                    cout << "finish reading mode " << i <<endl;
                    modes.push_back(tmp_mode);
                }

//                for( std::size_t w = 0; w <= i; w++ ){
//                    cout<< modes[w][0][0] << modes[w][0][1] << modes[w][0][2] << endl;
//                }
            }//generate non-zero mode
//            if( nmode_start != 0){
//                gen_zero_modes(modes);
//            } move to python
        }
        else {
            cout << "cannot open eigenvector file" << filename << endl;
            exit(EXIT_FAILURE);
        }
    }
    void schmidt(af::shared<normal_mode> & pall)
    {
        cout << pall.size() << endl;
        MMTBX_ASSERT(pall.size() == 6*(count_zero_modes+1));
        double anorm = 0.0;
        std::size_t padd = count_zero_modes*6;
        cout << padd << endl;
        cout << pall[0+padd].size() << endl;
        for(std::size_t i = 0; i < pall[0+padd].size() ; i++){
            anorm += pall[padd][i][0]*pall[padd][i][0];
        }
        MMTBX_ASSERT(anorm != 0);
        anorm = 1/sqrt(anorm);
        
        for(std::size_t i = 0; i < pall[padd].size() ; i++){
            pall[padd][i][0] *= anorm;
        }
        double rec[6][6];
        for(std::size_t i = 2; i < pall.size()/(count_zero_modes+1) ; i++){
            for(std::size_t j = 0; j < i-1; j++){
                rec[j][i] = 0.0;
                for(std::size_t k = 0; k < pall[i+padd].size(); k++){
                    rec[j][i] += pall[j+padd][k][0]*pall[i+padd][k][0] + pall[j+padd][k][1]*pall[i+padd][k][1] + pall[j+padd][k][2]*pall[i+padd][k][2];
                }
            }
            for(std::size_t k = 0; k < pall[i+padd].size(); k++){
                double aaa = 0.0;
                for(std::size_t j = 0; j < i - 1; j++){
                    aaa += (pall[j+padd][k][0] + pall[j+padd][k][1] + pall[j+padd][k][2])*rec[j][i];
                }
                pall[i+padd][k][0] -= aaa;
                pall[i+padd][k][1] -= aaa;
                pall[i+padd][k][2] -= aaa;
            }

            anorm = 0.0;
            for(std::size_t k = 0; k < pall[i+padd].size(); k++){
                anorm += pall[i+padd][k].length_sq();
            }
            MMTBX_ASSERT(anorm != 0);
            anorm = 1/sqrt(anorm);

            for(std::size_t k = 0; k < pall[i+padd].size(); k++){
                pall[i+padd][k][0] /= anorm;
                pall[i+padd][k][1] /= anorm;
                pall[i+padd][k][2] /= anorm;
            }
        }
    }
    void gen_zero_modes(af::shared<vec3<double> > const& sites_cart,
                        af::shared<double> const& weights)
    {
        MMTBX_ASSERT(sites_cart.size() == weights.size() );
        double xcm = 0.0;
        double ycm = 0.0;
        double zcm = 0.0;
        double tmass = 0.0;
        af::shared<vec3<double> > sites_cart_new(sites_cart.size(), af::init_functor_null<vec3<double> >());
        for( std::size_t i = 0; i < sites_cart.size() ; i++ ){
            double weight = weights[i];
            vec3<double> site = sites_cart[i];
            sites_cart_new[i] = sites_cart[i];
            xcm += site[0] * weight;
            ycm += site[1] * weight;
            zcm += site[2] * weight;
            tmass += weight;
        }
        MMTBX_ASSERT(tmass != 0);
        xcm /= tmass;
        ycm /= tmass;
        zcm /= tmass;
        for( std::size_t i = 0; i < sites_cart_new.size() ; i++){
            sites_cart_new[i][0] -= xcm;
            sites_cart_new[i][1] -= ycm;
            sites_cart_new[i][2] -= zcm;
        }
        for( std::size_t i = 0; i < 6 ; i++ ){
            af::shared<vec3<double> > null_mode( sites_cart.size(), af::init_functor_null<vec3<double> >());
            zero_modes.push_back(null_mode);
        }
        std::size_t padd = count_zero_modes*6;
        for( std::size_t i = 0; i < 6 ; i++ ){
            for( std::size_t j = 0 ; j < sites_cart.size(); j++ ){
                zero_modes[i+padd][j][0] = 0.0;
                zero_modes[i+padd][j][1] = 0.0;
                zero_modes[i+padd][j][2] = 0.0;
            }
        }
        for( std::size_t i = 0; i < sites_cart.size() ; i++ ){
            double sqrt_weight = sqrt(weights[i]);
            zero_modes[0+padd][i][0] = sqrt_weight;
            zero_modes[1+padd][i][1] = sqrt_weight;
            zero_modes[2+padd][i][2] = sqrt_weight;
            zero_modes[3+padd][i][1] = sqrt_weight*sites_cart_new[i][2];
            zero_modes[3+padd][i][2] = -sqrt_weight*sites_cart_new[i][1];
            zero_modes[4+padd][i][0] = -sqrt_weight*sites_cart_new[i][2];
            zero_modes[4+padd][i][2] = sqrt_weight*sites_cart_new[i][0];
            zero_modes[5+padd][i][0] = sqrt_weight*sites_cart_new[i][1];
            zero_modes[5+padd][i][1] = -sqrt_weight*sites_cart_new[i][0];
        }
        schmidt(zero_modes);
        count_zero_modes++;
        cout << count_zero_modes << endl;
    }
    af::shared<vec3<double> > return_modes( std::size_t i ) { return modes[i]; }    
    af::shared<vec3<double> > return_zero_modes( std::size_t i ) { return zero_modes[i]; }
    void print_eigenvector( std::size_t i )
    {
        for( std::size_t j = 0; j < modes[i].size(); j++ ){
            cout<< modes[i][j][0] << modes[i][j][1] << modes[i][j][2] << endl;
        }
    }
};

//af::versa<af::shared<sym_mat3<double> >, af::c_grid<2> > init_nm_adp(af::shared<normal_mode> const& modes,
//                                                                    af::shared<double> const& weights,
//                                                                    std::size_t n_modes,
//                                                                    bool zero_mode_flag)
//{
//    af::versa<af::shared<sym_mat3<double> >, af::c_grid<2> > adp_nma;
//    adp_nma.resize(af::c_grid<2>(n_modes, n_modes));
////    af::shared<sym_mat3<double> > adp_nma_i(modes[0].size(), af::init_functor_null<sym_mat3<double> >());
//    for(std::size_t n = 0; n < n_modes ; n++){
//        for(std::size_t m = n; m < n_modes ; m++){
//            MMTBX_ASSERT(modes[n].size() == weights.size());
//            MMTBX_ASSERT(modes[n].size() == modes[m].size());
//            af::shared<sym_mat3<double> > adp_nma_i(modes[n].size(), af::init_functor_null<sym_mat3<double> >());
//            if(zero_mode_flag){
//                if( n < 6 and m > 7 )
//                    continue;
//            }
//            for( std::size_t i = 0; i < modes[n].size() ; i++){
//                sym_mat3<double> adp;
//                adp[0] = modes[n][i][0]*modes[m][i][0];
//                adp[1] = modes[n][i][1]*modes[m][i][1];
//                adp[2] = modes[n][i][2]*modes[m][i][2];
//                adp[3] = (modes[n][i][0]*modes[m][i][1] + modes[n][i][1]*modes[m][i][0])/2.0;
//                adp[4] = (modes[n][i][0]*modes[m][i][3] + modes[n][i][3]*modes[m][i][0])/2.0;
//                adp[5] = (modes[n][i][2]*modes[m][i][3] + modes[n][i][3]*modes[m][i][2])/2.0;
//                adp = adp/weights[i];
//                adp_nma_i[i] = adp;
//            }
//            adp_nma(n, m) = adp_nma_i;
//        }
//    }
//    return adp_nma;
//}

//class uaniso_from_s {
//private:
//    af::versa<double, af::c_grid<2> > s;
//    af::versa<double, af::c_grid<2> > sigma;

}}//namespace mmtbx::nm
#endif //MMTBX_NM_H
