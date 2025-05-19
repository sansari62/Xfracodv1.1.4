#include <stdafx.h>
#include <CommonPara.h>

using namespace CommonPara_h::comvar;



void save(ofstream& file10)
{       
    const char* File_ID = "FRACOD SAVED FILE";
    
    file10 << File_ID << std::endl;
    file10 << numbe << " " << no << " " << numbe_old << std::endl;
    file10 << title << std::endl;

    for (int mm = 0; mm < 10; mm++)
    {
        rock1[mm].save_to_file(file10);        
    }
    s2us.save_to_file(file10);
    symm.save_to_file(file10);
    s4.save_to_file(file10);
    s5u.save_to_file(file10);

    for (int m = 0; m < numbe; ++m)
    {
        elm_list[m].save_to_file(file10,m);
    }   

    for (int m = 0; m < no; ++m) {
        tips[m].save_to_file(file10);
    }

    for (int m = 0; m < 20; ++m) {
        s8[m].save_to_file(file10);
    }
    file10 << numbe << no << delta << w0 << w1 << ni << nc << numbe_old << std::endl;
    dispwin.save_to_file(file10);
    if(ntunnel>0)
     tunnl.save_to_file(file10);
    if(nellipse>0)
        file10 << nellipse << std::endl;
    for (int m = 0; m < nellipse; ++m) {
        ellip_list[m].save_to_file(file10);
    }
   
    insituS.save_to_file(file10);

    file10 << w0 << " " << w1 << " " << ni << std::endl;
    file10 << mcyc0 << " " << mcyc<< std::endl;
    file10 << lastinput << " " << ktipgrow << " " << StopReturn << " " << line << " " << ID_dtip << std::endl;

    s15.save_to_file(file10);

    watercm.save_to_file(file10);

    file10 << mat_lining << std::endl;

    file10 << n_it << std::endl;
    file10 << k_num << " " << d_max << std::endl;
    file10 << ihist << std::endl;

    for (int m = 0; m < ihist; ++m) {       //ihist and lhist are not included here , should  think about Sara!
        mpoint_list[m].save_to_file(file10);
    }
    file10 << lhist << std::endl;
    for (int m = 0; m < lhist; ++m) {
        mline_list[m].save_to_file(file10);
    }
    creep.save_to_file(file10);
    file10 << mf << std::endl;
    for (int m = 0; m < mf; ++m)
    {
        init_point[m].save_to_file(file10);
    }  
    for (int m = 0; m < numbe; ++m)
    {
        file10<< joint[m].aperture0 <<" "<< joint[m].aperture_r<<endl;
    }
    file10 << perm.viscosity << " " << perm.density << " " << perm.perm0 << std::endl;

    file10 << factors.factor_f << " " << factors.factor_e << " " << factors.tolerance << std::endl;

    file10 << exca.ID_Exca << " " << exca.d_wall << " " << exca.rand_e << std::endl;
   
    file10.close();
    return;
}
