#include <stdafx.h>
#include <CommonPara.h>

using namespace CommonPara_h::comvar;



void save(ofstream& file10)
{       
    const char* File_ID = "FRACOD SAVED FILE";
    
    file10 << File_ID << std::endl;
    file10 << numbe << " " << no << " " << numbe_old << " " 
        << title<< std::endl;

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




void restore()
{
    std::ifstream inputFile("filename.dat");
    if (!inputFile.is_open()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    inputFile >> numbe >> no >> numbe_old >> title;

    for (int mm = 0; mm < 10; ++mm)
    {
        rock1[mm].read_from_file(inputFile);
    }

    s2us.read_from_file(inputFile);
    symm.read_from_file(inputFile);
    
    s4.read_from_file(inputFile);

    s5u.read_from_file(inputFile);

    for (int m = 0; m < numbe; ++m)
    {
        elm_list[m].read_from_file(inputFile, m);
    }

    for (int m = 0; m < no; ++m) {
        tips[m].read_from_file(inputFile);
    }

    for (int m = 1; m <= 20; ++m) {
        s8[m].read_from_file(inputFile);
    }

    inputFile >> numbe >> no >> delta >> w0 >> w1 >> ni >> nc >> numbe_old;
    dispwin.read_from_file(inputFile);

    tunnl.read_from_file(inputFile);

    /*int ntunnel;
    inputFile >> ntunnel;
    tunnl_list.resize(ntunnel);
    for (int m = 0; m < ntunnel; ++m)
        {
             tunnl_list[m].read_from_file1(inputFile);
        }*/


        //int nellipse;
    inputFile >> nellipse;
    ellip_list.resize(nellipse); //size or resize? Sara!
    for (int m = 0; m < nellipse; ++m) {
        ellip_list[m].read_from_file(inputFile);
    }

    /*for (int m = 1; m <= 20; ++m) {     //should think more  in tunnel
        inputFile >> xpp[m] >> ypp[m] >> dd[m];
    }*/

    insituS.read_from_file(inputFile);

    inputFile >> w0 >> w1 >> ni;
    inputFile >> mcyc0 >> mcyc;

    inputFile >> lastinput >> ktipgrow >> StopReturn >> line >> ID_dtip;

    s15.read_from_file(inputFile);

    watercm.read_from_file1(inputFile);

    inputFile >> mat_lining;

    inputFile >> n_it;

    inputFile >> k_num >> d_max;

    for (int m = 0; m < 10; ++m) {       //ihist and lhist are not included here , should  think about
              
        inputFile >> mpoint_list[m].xmon >> mpoint_list[m].ymon;      

    }

    for (int m = 0; m < 10; ++m) {        
        inputFile >> mline_list[m].x1l >> mline_list[m].y1l >> mline_list[m].x2l >> 
            mline_list[m].y2l >> mline_list[m].npl;
    }

    creep.read_from_file(inputFile);

    inputFile >> mf;
    for (int m = 0; m < 500; ++m)
        init_point[m].read_from_file(inputFile);

    for (int m = 0; m < numbe; ++m) {
        inputFile >> joint[m].aperture0 >> joint[m].aperture_r;  
    }

    inputFile >> perm.viscosity >> perm.density >> perm.perm0;

    inputFile >> factors.factor_f >> factors.factor_e >> factors.tolerance;

    inputFile >> exca.ID_Exca >> exca.d_wall >> exca.rand_e;

    //call SendWindowText('Save file reading error, file might be saved by earlier version'//CHAR(0))
    //    StopReturn = .true.
    //    return
    //    !stop
    //    GOTO 800

    //    600 call SendWindowText('File saved by earlier version, no crack initiation assumed'//CHAR(0))
    //        StopReturn = .true.
    //        return
    //        !stop

    lastinput = "endf";

    return;
}



