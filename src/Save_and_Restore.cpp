#include <stdafx.h>
#include <CommonPara.h>

using namespace CommonPara_h::comvar;



void save(ofstream& file10)
{       
    const char* File_ID = "FRACOD_SAVED_FILE";
    
    file10 << File_ID << std::endl;
    file10 << numbe << " " << no << " " << numbe_old << " " 
        << title<< std::endl;
    file10 << pi << " " << irock<< std::endl;
    for (int mm = 0; mm < 10; mm++)
    {
        rock1[mm].save_to_file(file10);        
    }
    s2us.save_to_file(file10);
    symm.save_to_file(file10);
    s4.save_to_file(file10);
    file10 << "ends4ishere" << endl;
    s5u.save_to_file(file10);

    for (int m = 0; m < numbe; ++m)
    {
        elm_list[m].save_to_file(file10,m);
    }   

    for (int m = 0; m < no; ++m) {
        tips[m].save_to_file(file10);
    }

    for (int m = 1; m <= 20; ++m) {
        s8[m].save_to_file(file10);
    }
    file10 << numbe << " " << no << " " << delta << " " << w0 <<
        " " << w1 << " " << ni << " " << nc << " " <<
        numbe_old << std::endl;
    dispwin.save_to_file(file10);
    file10 << ntunnel<< std::endl;
    if(ntunnel>0)
     tunnl.save_to_file(file10);
    file10 << nellipse << std::endl;
    if(nellipse>0)
    {
        for (int m = 0; m < nellipse; ++m) {
            ellip_list[m].save_to_file(file10);
        }
    }   
   
    insituS.save_to_file(file10);

    file10 << w0 << " " << w1 << " " << ni << std::endl;
    file10 << mcyc0 << " " << mcyc<< std::endl;
    file10 << lastinput << " " << ktipgrow << " " << StopReturn << " " << line << 
        " " << ID_dtip << std::endl;

    s15.save_to_file(file10);
    file10 << water_mod<<endl;
    if (water_mod)
        watercm.save_to_file(file10);

    file10 << mat_lining << std::endl;
    file10 << n_it << std::endl;
    file10 << k_num << " " << d_max << std::endl;
    file10 << ihist << std::endl;

    for (int m = 0; m < ihist; ++m) { 
        mpoint_list[m].save_to_file(file10);
    }
    file10 << lhist << std::endl;
    for (int m = 0; m < lhist; ++m) {
        mline_list[m].save_to_file(file10);
    }
    file10 << creep.ID_creep<<endl;
    if(creep.ID_creep!=0)
        creep.save_to_file(file10);
    file10 << mf << std::endl;
    if(mf>0)
    {
        for (int m = 0; m < mf; ++m)
        {
            init_point[m].save_to_file(file10);
        }
    }
    
    file10 << perm.viscosity << " " << perm.density << " " << perm.perm0 << std::endl;

    file10 << factors.factor_f << " " << factors.factor_e << " " << factors.tolerance << std::endl;

    file10 << exca.ID_Exca << " " << exca.d_wall << " " << exca.rand_e << std::endl;
   
    file10.close();
    return;
}




void restore(string filename)
{
    std::ifstream inputFile(filename);
    if (!inputFile.is_open()) {
        std::cerr << "Error in opening save_file!" << std::endl;
        return;
    }
    string lineData;
    getline(inputFile, lineData);
    if (lineData != "FRACOD_SAVED_FILE")
    {
        MessageBox(nullptr, L"File is not a FRACOD saved file!press OK to quit.",
            L"Message!", MB_OK);
        exit(0);
    } 
    getline(inputFile, lineData);
    std::stringstream ss(lineData);
    ss >> numbe >> no >> numbe_old >> title;
    inputFile >> pi>> irock;
    for (int mm = 0; mm < 10; ++mm)
    {
        rock1[mm].read_from_file(inputFile);
    }

    s2us.read_from_file(inputFile);
    symm.read_from_file(inputFile);    
    s4.read_from_file(inputFile);
    string key;
    inputFile >> key;
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

    inputFile >> numbe >> no >> delta >> w0 >> w1 >> ni >> nc >> 
        numbe_old;
    dispwin.read_from_file(inputFile);
    inputFile >> ntunnel;
    if (ntunnel>0)
    {
        tunnl.read_from_file(inputFile);
    }    
    inputFile >> nellipse;
    if (nellipse > 0)
    {
        for (int m = 0; m < nellipse; ++m) {
            ellip_list[m].read_from_file(inputFile);
        }
    }
    insituS.read_from_file(inputFile);
    inputFile >> w0 >> w1 >> ni;
    inputFile >> mcyc0 >> mcyc;
    inputFile >> lastinput >> ktipgrow >> StopReturn >> line >> ID_dtip;

    s15.read_from_file(inputFile);
    inputFile >> water_mod;
    if(water_mod)
        watercm.read_from_file1(inputFile);

    inputFile >> mat_lining >> n_it;//instead of seprate reaading  
    inputFile >> k_num >> d_max;
    inputFile >> ihist;
    for (int m = 0; m < ihist; ++m) { 
              
        inputFile >> mpoint_list[m].xmon >> mpoint_list[m].ymon; 
    }
    inputFile >> lhist;
    for (int m = 0; m < lhist; ++m) {        
        inputFile >> mline_list[m].x1l >> mline_list[m].y1l >> mline_list[m].x2l >> 
            mline_list[m].y2l >> mline_list[m].npl;
    }
    inputFile >> creep.ID_creep;
    if (creep.ID_creep != 0)
        creep.read_from_file(inputFile);
    inputFile >> mf;
    if (mf > 0)
    {
        for (int m = 0; m < mf; ++m)
            init_point[m].read_from_file(inputFile);
    }    

    inputFile >> perm.viscosity >> perm.density >> perm.perm0;

    inputFile >> factors.factor_f >> factors.factor_e >> factors.tolerance;

    inputFile >> exca.ID_Exca >> exca.d_wall >> exca.rand_e;
    
    //lastinput = "endf";
    return;
}


