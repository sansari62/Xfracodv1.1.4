#include <stdafx.h>
#include <CommonPara.h>
#include <WinInterface.h>
#include<Work.h>
#include<ExcavationCracks.h>
#include<InputCheck.h>
#include<iostream>
#include<functional>
#include<chrono>
#include <unordered_map>
#include<Save_and_Restore.h>
#include <charconv>
#include<Tip.h>



using namespace CommonPara_h::comvar;







void processEdge()
{
       int material = 1;
       float gradsy = 0, gradny = 0;        
       int num, kode;
       float xbeg, ybeg, xend, yend, bvs, bvn;
       string lineData;
      
       try
       {
           getline(inFile, lineData);
           std::stringstream ss(lineData);

           
           ss >> num >> xbeg >> ybeg >> xend >> yend >> kode >> bvs >> bvn >> material >> gradsy >> gradny;
         
           file2 << "\n Straight boundary -- from: x,y = "
               << std::fixed << std::setprecision(4) << std::setw(10) << xbeg << ", "
               << std::fixed << std::setprecision(4) << std::setw(10) << ybeg << std::endl
               << "                    to        = "
               << std::fixed << std::setprecision(4) << std::setw(10) << xend << ", "
               << std::fixed << std::setprecision(4) << std::setw(10) << yend << std::endl
               << "      shear stress or disp    = "
               << std::scientific << std::setprecision(4) << std::setw(10) << bvs << std::endl
               << "       normal stress or disp  = "
               << std::scientific << std::setprecision(4) << std::setw(10) << bvn << std::endl
               << "      shear value gradient in y direction = "
               << std::scientific << std::setprecision(4) << std::setw(10) << gradsy << std::endl
               << "      normal value gradient in y direction = "
               << std::scientific << std::setprecision(4) << std::setw(10) << gradny << std::endl
               << "         number of elements   = " << std::setw(4) << num << std::endl
               << "         material              = " << std::setw(4) << material << std::endl;


       }
       catch (std::ifstream::failure e)
       {
           std::cerr << "Exception opening/reading/closing file:in processEdge\n";

       }

        switch (kode)
        {
        case 1:
        case 11:
            bund_list[nb] = Edge(material, num, kode, xbeg, ybeg, xend, yend, bvs, bvn, 0, 0, gradsy, gradny);           
            break;
        case 2:
        case 12:
            bund_list[nb] = Edge(material, num, kode, xbeg, ybeg, xend, yend, 0, 0, bvs, bvn, gradsy, gradny);
            
            break;
        case 3:
        case 13:
            bund_list[nb] = Edge(material, num, kode, xbeg, ybeg, xend, yend, 0, bvn, bvs, 0, gradsy, gradny);            
            break;
        case 4:
        case 14:
            bund_list[nb] = Edge(material, num, kode, xbeg, ybeg, xend, yend, bvs, 0, 0, bvn, gradsy, gradny);          
           
        }       
        nb++;
        return;
}




void processGost()
{
    int material = 1;
    float gradsy = 0, gradny = 0;
    int num, kode;
    float xbeg, ybeg, xend, yend, bvs, bvn;

    try
    {

        inFile >> num >> xbeg >> ybeg >> xend >> yend >> kode >> bvs >> bvn >> material >> gradsy >> gradny;  
        file2 << "Gost element     -- from: x,y = "
            << std::fixed << std::setprecision(4) << std::setw(10) << xbeg << ", "
            << std::fixed << std::setprecision(4) << std::setw(10) << ybeg << std::endl
            << "                    to        = "
            << std::fixed << std::setprecision(4) << std::setw(10) << xend << ", "
            << std::fixed << std::setprecision(4) << std::setw(10) << yend << std::endl
            << "      shear stress or disp    = "
            << std::scientific << std::setprecision(4) << std::setw(10) << bvs << std::endl
            << "       normal stress or disp  = "
            << std::scientific << std::setprecision(4) << std::setw(10) << bvn << std::endl
            << " shear value gradient in y direction = "
            << std::scientific << std::setprecision(4) << std::setw(10) << gradsy << std::endl
            << " normal value gradient in y direction = "
            << std::scientific << std::setprecision(4) << std::setw(10) << gradny << std::endl
            << "         number of elements   = " << std::setw(4) << num << std::endl
            << "         material              = " << std::setw(4) << material << std::endl;

    }
    catch (std::ifstream::failure e)
    {
        MessageBox(nullptr, L"Error in input File,Gost definition!", L"Error", MB_OK);

    }

    Edge new_sBound(material, num, 7, xbeg, ybeg, xend, yend, bvs, bvn, 0, 0, gradsy, gradny);
    bund_list[nb] = new_sBound;       
    nb++;
    return;
}





void ToLowerCase(string&  str) {
  
    for (char& c: str) 
    {
        c = std::tolower(c);
    }
}




void processFracture()
{
    int material = 1;      
    int num, itype, jmat;
    float xbeg, ybeg, xend, yend;
    string lineData;
    try
    {
       getline(inFile, lineData);
       std::stringstream ss(lineData);
       ss >> num >> xbeg >> ybeg >> xend >> yend >> itype >> jmat >> material;
       file2 << "Fracture ---------- from: x,y = "
            << std::fixed << std::setprecision(4) << std::setw(10) << xbeg << ", "
            << std::fixed << std::setprecision(4) << std::setw(10) << ybeg << std::endl
            << "                    to        = "
            << std::fixed << std::setprecision(4) << std::setw(10) << xend << ", "
            << std::fixed << std::setprecision(4) << std::setw(10) << yend << std::endl
            << "         number of elements   = " << std::setw(4) << num << std::endl
            << "         material              = " << std::setw(4) << material << std::endl;
            
        Fracture frac(xbeg, ybeg, xend, yend, material, num, jmat);         //dfine new fracture
        frac_list[nf] = frac;
        frac_list[nf].bound_type = 5;
        nf++;
    }
    catch (std::ifstream::failure e)
    {
        MessageBox(nullptr, L"Error in input File,Fracture definition!", L"Error", MB_OK);
    }   
     return;
}




 void processLinearInterface()
    {
     int num, mat1, mat2;
     float xbeg, ybeg, xend, yend, bvs = 0, bvn = 0;
     string lineData;
        try
        {    
            getline(inFile, lineData);
            std::stringstream ss(lineData);
            ss >> num >> xbeg >> ybeg >> xend >> yend >> mat1 >> mat2;            
            file2 << "\n Interface -------- from: x, y =  " << xbeg << "," << ybeg << "\n" <<
                "                    to      =  " << xend << "," << yend << "\n" <<
                "         number of elements   = " << num << "\n";
            
            Edge_interface new_lin_intfce(xbeg, ybeg, xend, yend, mat1, mat2,num); 
            lin_intrfce_list[npli] = new_lin_intfce;
            npli++;
            multi_region = true;

        }
        catch (std::ifstream::failure e)
        {
            MessageBox(nullptr, L"Error in input File,linear interface definition!", L"Error", MB_OK);

        }
    }



    void processArcInterface()
    {
        string lineData;
        int num = 0, mat1 = 0, mat2 = 0;
        float xcen = 0, ycen = 0, diam = 0, ang1 = 0, ang2 = 0, bvs = 0, bvn = 0;
        try
        {            
            getline(inFile, lineData);
            std::stringstream ss(lineData);
            ss >> num >> xcen >> ycen >> diam >> ang1 >> ang2 >> mat1 >> mat2; 

            file2 << "\n ARC Int- centre position: x,y  = "
                << std::fixed << std::setprecision(4) << std::setw(10) << xcen << ", "
                << std::fixed << std::setprecision(4) << std::setw(10) << ycen << std::endl
                << "         diameter              = "
                << std::fixed << std::setprecision(4) << std::setw(10) << diam << std::endl
                << "         start angle           = "
                << std::fixed << std::setprecision(4) << std::setw(10) << ang1 << std::endl
                << "         end angle             = "
                << std::fixed << std::setprecision(4) << std::setw(10) << ang2 << std::endl
                << "         shear stress or disp  = "
                << std::scientific << std::setprecision(4) << std::setw(10) << bvs << std::endl
                << "         normal stress or disp = "
                << std::scientific << std::setprecision(4) << std::setw(10) << bvn << std::endl
                << "         material one side     = " << std::setw(4) << mat1 << std::endl
                << "         material one side     = " << std::setw(4) << mat2 << std::endl
                << "         material other side   = " << std::setw(4) << num << std::endl;

            
            Arch_interface newArcIntrfce(xcen, ycen, diam / 2., ang1 * pi / 180.,
                ang2 * pi / 180., num, mat1, mat2);
            arc_intrface_list[nq] = newArcIntrfce;
            int kode = 6;   
            bvs = 0;
            bvn = 0;
            nq = nq + 1;
            multi_region = true;
        }
        catch (std::ifstream::failure e)
        {
            MessageBox(nullptr, L"Error in input File,Arc definition!", L"Error", MB_OK);

        }
        return;                   
    
    }




    void processWaterPressure()
    {
        string id;
        string lineData;
        water_mod = true;
        try
        {
            getline(inFile, lineData);
            std::stringstream ss(lineData);

            ss >> id ;           
            ToLowerCase(id);
           // id = id.substr(0, 4);
            if (id == "hole")
            {
                getline(inFile, lineData);
                std::stringstream ss(lineData);
                
                int iwhole = watercm.iwhole;               
                ss >> watercm.w_xc[iwhole] >> watercm.w_yc[iwhole] >> watercm.w_d[iwhole]
                    >> watercm.wph[iwhole];

                watercm.ID_range = 1;
                watercm.iwhole++;
                file2 << "Water pressure in hole" << endl
                    << "     center x = " << watercm.w_xc[iwhole] << endl
                    << "     center y = " << watercm.w_yc[iwhole] << endl
                    << "     diameter = " << watercm.w_d[iwhole] << endl
                    << "     pressure = " << watercm.wph[iwhole] << endl << endl;
            }
            else if (id == "rect")
            {
                getline(inFile, lineData);
                std::stringstream ss(lineData);
                int iwrect = watercm.iwrect;         

                ss >> watercm.w_x1[iwrect] >> watercm.w_x2[iwrect] >> watercm.w_y1[iwrect] >>
                    watercm.w_y2[iwrect] >> watercm.wpr[iwrect];

                watercm.ID_range = 2;
                watercm.iwrect++;
                file2 << "Water pressure in rectangular opening" << endl
                    << "     x1 = " << watercm.w_x1[iwrect] << endl
                    << "     x2 = " << watercm.w_x2[iwrect] << endl
                    << "     y1 = " << watercm.w_y1[iwrect] << endl
                    << "     y2 = " << watercm.w_y2[iwrect] << endl
                    << "     pressure = " << watercm.wpr[iwrect] << endl << endl;
            }
        }
        catch (std::ifstream::failure e)
        {
            MessageBox(nullptr, L"Error in input File,water definition!", L"Error", MB_OK);

        }
        return;
    }


    void processRect()
    {
        string lineData;
        try
        {
            getline(inFile, lineData);
            std::stringstream ss(lineData);           
            int iwrect = watercm.iwrect;

             ss >> watercm.w_x1[iwrect] >> watercm.w_x2[iwrect] >> watercm.w_y1[iwrect] >>
            watercm.w_y2[iwrect] >> watercm.wpr[iwrect];

            watercm.ID_range = 2;
            watercm.iwrect++;
            file2 << "Water pressure in rectangular opening" << endl
                << "     x1 = " << watercm.w_x1[iwrect] << endl
                << "     x2 = " << watercm.w_x2[iwrect] << endl
                << "     y1 = " << watercm.w_y1[iwrect] << endl
                << "     y2 = " << watercm.w_y2[iwrect] << endl
                << "     pressure = " << watercm.wpr[iwrect] << endl << endl;
        }
        catch (std::ifstream::failure e)
        {
            MessageBox(nullptr, L"Error in input File,rect definition!", L"Error", MB_OK);

        }
        return;
    }




    void processHole()
    {
        string lineData;
        try
        {
            getline(inFile, lineData);
            std::stringstream ss(lineData);

            int iwhole = watercm.iwhole;
            ss >> watercm.w_xc[iwhole] >> watercm.w_yc[iwhole] >> watercm.w_d[iwhole]
                >> watercm.wph[iwhole];

            watercm.ID_range = 1;
            watercm.iwhole++;
            file2 << "Water pressure in hole" << endl
                << "     center x = " << watercm.w_xc[iwhole] << endl
                << "     center y = " << watercm.w_yc[iwhole] << endl
                << "     diameter = " << watercm.w_d[iwhole] << endl
                << "     pressure = " << watercm.wph[iwhole] << endl << endl;
        }
        catch (std::ifstream::failure e)
        {
            MessageBox(nullptr, L"Error in input File,hole definition!", L"Error", MB_OK);

        }
        return;

    }




    void processPermeability()
    {
        string lineData;
        try
        {
            getline(inFile, lineData);
            std::stringstream ss(lineData);        
            ss >> perm.viscosity >> perm.density >> perm.perm0;
            file2 << "permeability parameters" << endl
                << "     fluid viscosity = " << perm.viscosity << " kg/(m s)" << endl
                << "     fluid density = " << perm.density << " kg/m3" << endl
                << "     intact rock hydraulic conductivity = " << perm.perm0 << " m/s" << endl << endl;
        }
        catch (std::ifstream::failure e)
        {
            MessageBox(nullptr, L"Error in input File,permeability definition!", L"Error", MB_OK);
        }
        return;
    }




    bool isValidNumber(const std::string& str) {
        int value;
        auto [ptr, ec] = std::from_chars(str.data(), str.data() + str.size(), value);
        return ec == std::errc() && ptr == str.data() + str.size();
    }


    void start_calculation()
    {   
        string lineData;        
        nc++;
        string  tem = "";
        try
        {
            getline(inFile, lineData);
            std::stringstream ss(lineData);
            ss >> tem;
        }
        catch (std::ifstream::failure e)
        {
            MessageBox(nullptr, L"Error in input File,touc definition!", L"Error", MB_OK | MB_ICONERROR);
        }
        
        if (tem != " ") 
        {
            if (!isValidNumber(tem)) {
                MessageBox(nullptr, L"Expected a number after 'cycle'.Please add the number of cycles on the next line in the input file!", L"Error", MB_OK | MB_ICONERROR);
                exit(EXIT_FAILURE);
            }
            
            //file2 << tem << endl;              
            mcyc0 += stoi(tem);               
        }
        //----------------------------
        if (mcyc != 0) return;            
        s5u.dtol = 0.0001;            
           
        if (nc == 1)
        {
            inputcheck();
            input_tip_check();         
            file2 << endl;
            file2 << " " << endl
                << " boundary element data." << endl
                << endl
                << "element   kode   X(center)  Y(center)  length  angle  us or sigma-s   un or sigma-n   material No." << endl;

            float amin = 10000.0;

            // Calculate min
            for (int m = 0; m < numbe; ++m)
            {
                amin = min(amin, elm_list[m].a);
            }
            s5u.dtol = amin / 1000.0;
            //if (water_mod)
               // water();   
            // Write  elementdata to file
            for (int m = 0; m < numbe; ++m)
            {
                BoundaryElement& elm = elm_list[m];    
                float size = 2.0 * elm.a;
                float angle = 180.0 * static_cast<float>(atan2(elm.sinbet, elm.cosbet)) / pi;
                file2 << std::setw(7) << m + 1
                    << std::setw(6) << elm.kod
                    << std::setw(9) << std::fixed << std::setprecision(4) << elm.xm
                    << std::setw(9) << std::fixed << std::setprecision(4) << elm.ym
                    << std::setw(9) << std::fixed << std::setprecision(4) << size
                    << std::setw(8) << std::fixed << std::setprecision(2) << angle
                    << std::setw(13) << std::scientific << s4.b0[2 * m]
                    << std::setw(13) << std::scientific << s4.b0[2 * m + 1]
                    << std::setw(5) << elm.mat_no << std::endl;
            }

            float aks_bb = 0, akn_bb = 0, phi_bb = 0; float coh_bb = 0.0;
            float phid_bb = 0, ap_bb = 0, apr_bb = 0;

            // Call stiffness_bb for Tensile fractures
            stiffness_bb(aks_bb, akn_bb, phi_bb, coh_bb, phid_bb, ap_bb, apr_bb, 1, 1);
            file2 << "New fresh fracture properties --- Tensile fractures" << std::endl
                << "               ks = " << std::scientific << std::setprecision(4) << aks_bb << std::endl
                << "               kn = " << std::scientific << std::setprecision(4) << akn_bb << std::endl
                << "              phi = " << std::fixed << std::setprecision(1) << (phi_bb * 180 / pi) << std::endl
                << "        cohesion = " << std::scientific << std::setprecision(4) << coh_bb << std::endl
                << "   dilation angle = " << std::fixed << std::setprecision(1) << (phid_bb * 180 / pi) << std::endl
                << " initial aperture = " << std::scientific << std::setprecision(4) << ap_bb << std::endl
                << "residual aperture = " << std::scientific << std::setprecision(4) << apr_bb << std::endl;

            // Call stiffness_bb for Shear fractures
            stiffness_bb(aks_bb, akn_bb, phi_bb, coh_bb, phid_bb, ap_bb, apr_bb, 2, 1);
            file2 << "\n New fresh fracture properties --- Shear fractures" << std::endl
                << "               ks = " << std::scientific << std::setprecision(4) << aks_bb << std::endl
                << "               kn = " << std::scientific << std::setprecision(4) << akn_bb << std::endl
                << "              phi = " << std::fixed << std::setprecision(1) << (phi_bb * 180 / pi) << std::endl
                << "        cohesion = " << std::scientific << std::setprecision(4) << coh_bb << std::endl
                << "   dilation angle = " << std::fixed << std::setprecision(1) << (phid_bb * 180 / pi) << std::endl
                << " initial aperture = " << std::scientific << std::setprecision(4) << ap_bb << std::endl
                << "residual aperture = " << std::scientific << std::setprecision(4) << apr_bb << std::endl;
        }
                
        return;
    }




   void saveData()
   {   
       string lineData;
       
       //string id ;
       string tem = "  ";
       try
       {
           getline(inFile, lineData);           
           std::stringstream ss(lineData);
           ss >> tem;         
           if (tem == "  ")
           {
               MessageBox(nullptr, L"Incomplete input - Give a file name for saving", L"Error", MB_OK);
           }
           else
           {
               ofstream file10(tem, ios::binary);
               save(file10);
               file10.close();
           }
       }
       catch (std::ifstream::failure e)
       {
           MessageBox(nullptr, L"Error in input File,save definition!", L"Error", MB_OK);

       }      
        return;
    }




   void endFile()
   {       
       lastinput = "endf";  

       return;
   }




   void processExcavation_Cracks()
   {
       exca.ID_Exca = 1;
       string lineData;
       try
       {
           getline(inFile, lineData);
           std::stringstream ss(lineData);      
           ss >> exca.d_wall >> exca.rand_e;
           file2 << "Excavation induced cracks: " << std::endl;
           file2 << "     Distance into wall = " << exca.d_wall << std::endl;
           file2 << "     Percentage of points with cracks = " << exca.rand_e * 100 << "%" << std::endl;
       }
       catch (std::ifstream::failure e)
       {
           MessageBox(nullptr, L"Error in input File,escavation_crack definition!", L"Error", MB_OK);

       }           

       return;
   }




   void error_Handling()
   {
       std::cerr << "Error reading data from file" << std::endl;
       inFile.close();
       return ;
   }
    



   void symmetryProcessing()
   {     

        string lineData;
        try
        {
            getline(inFile, lineData);
            std::stringstream ss(lineData);
            ss >> symm.ksym >> symm.xsym >> symm.ysym;
        }
        catch (std::ifstream::failure e)
        {
            MessageBox(nullptr, L"Error in input File,symm definition!", L"Error", MB_OK);

        }                

        switch (symm.ksym)
        {
        case 0:
            file2 << "No symmetry" << std::endl;
            break;
        case 1:
            file2 << "Symmetry against vertical line x = " << symm.xsym << std::endl;
            break;

        case 2:
            file2 << "Symmetry against horizontal line y = " << symm.ysym << std::endl;
            break;

        case 3:
            file2 << "Symmetry against point x = " << symm.xsym << std::endl <<
                " y = " << symm.ysym << std::endl;
            break;

        case 4:
            file2 << "Symmetry against vertical line x = " << symm.xsym << std::endl <<
                " and horizontal line y = " << symm.ysym << std::endl;
            break;
            }
              
        return;
   }



   void frac_energy_processing()
   {
       string lineData;
       float gic, giic;
       int material = 1;      
       try
       {
           getline(inFile, lineData);           
           std::stringstream ss(lineData);
           ss >> gic >> giic >> material;
           int mm = material;
           file2 << "ROCK MATERIAL = " << material << std::endl;
           file2 << "    Mode I fracture energy = " << std::scientific << std::setprecision(3) << gic << std::endl;
           file2 << "    Mode II fracture energy = " << std::scientific << std::setprecision(3) << giic << std::endl;
           rock1[mm].akic = sqrt(gic * rock1[mm].e / (1 - rock1[mm].pr * rock1[mm].pr));
           rock1[mm].akiic = sqrt(giic * rock1[mm].e / (1 - rock1[mm].pr * rock1[mm].pr));
       }
       catch (std::ifstream::failure e)
       {
           MessageBox(nullptr, L"Error in input File,touc definition!", L"Error", MB_OK);

       }       
       return;
   }



   void frac_toughness_process()
   {
       string lineData;      
       int mm = 0;
       int material = 1;
       float akic0 = 0, akiic0 = 0;      
       try
       {
           getline(inFile, lineData);
           std::stringstream ss(lineData);
           ss >> akic0 >> akiic0 >> material ;// >> material;    I removed temp vars here and assign directly the rock attrbs. Sara! check again for errors
           file2 << "ROCK MATERIAL = " << material << std::endl;
           file2 << "    Mode I toughness = " << akic0 << std::endl;
           file2 << "    Mode II toughness = " << akiic0 << std::endl;
           mm = material;
           rock1[mm].akic = akic0;
           rock1[mm].akiic = akiic0;
       }
       catch (std::ifstream::failure e)
       {
           MessageBox(nullptr, L"Error in input File,toughness definition!", L"Error", MB_OK);

       }
       return;
   }




   void elastic_modul_processing()
   {
       string lineData;       
       int material = 1;
       int mm = 0;
       float pr0, e0;       

       try
       {
           getline(inFile, lineData);
           std::stringstream ss(lineData);
           ss>>pr0 >> e0 >> material;
           mm = material;
           rock1[mm].e = e0;
           rock1[mm].pr = pr0;
           file2 << "ROCK MATERIAL = " << mm << std::endl
            << " Poisson's ratio = " << std::fixed << std::setprecision(2) << rock1[mm].pr << std::endl
            << " Young's modulus = " << std::scientific << std::setprecision(4) << rock1[mm].e << std::endl;
       }
       catch (std::ifstream::failure e)
       {
           MessageBox(nullptr, L"Error in input File,modulus definition!", L"Error", MB_OK);

       }
       return;
   }





   void dbouProcess()
   {
       file9.open(dir + L"/Cbound.dat", std::ios::in | std::ios::out | std::ios::trunc);
       if (!file9) {
           std::cerr << "Error opening boundary file!" << std::endl;
           return ;
       }
       //---------------straight boundary value increment ------------------      
           float x1, x2, y1, y2, dss, dnn;
           string lineData;
           try
           {
               getline(inFile, lineData);
               std::stringstream ss(lineData);           
               ss >> x1 >> x2 >> y1 >> y2 >> dss >> dnn;               
               file2 << "boundary value increment, x1,x2,y1,y2,dss,dnn =" << endl;
               file2 << std::fixed << std::setprecision(2) << std::setw(10) << x1 << std::setw(10) << x2 << std::setw(10) << y1 << std::setw(10)
                   << y2 << std::setw(10) << dss << std::setw(10) << dnn << std::endl;
               file9 << "dbou" << std::endl;
               file9 << x1 << " " << x2 << " " << y1 << " " << y2 << " " << dss << " " << dnn << std::endl;
              
           }
           catch (std::ifstream::failure e) 
           {
               MessageBox(nullptr, L"Error in input File,dboundary definition!", L"Error", MB_OK);

           }
           insituS.incres = 3;    
           file9.close();
           return;
   }






   void darcProcess()
   {
       //---------------arch boundary value increment ------------------
      
           float xcen, ycen, diam1, diam2, ang1, ang2, dss, dnn;
           string lineData;
           if (!file9) {
               std::cerr << "Error opening file boundary!" << std::endl;
               return;
           }
           try
           {
               getline(inFile, lineData);
               std::stringstream ss(lineData);
                ss >> xcen >> ycen >> diam1 >> diam2 >> ang1 >> ang2 >> dss >> dnn;
               file2 << "Arch boundary value increment, xc, xy, d1, d2, ang1, ang2, dss, dnn =\n"
                   << std::fixed << std::setprecision(2) << std::setw(10)<<xcen << std::setw(10) << ycen << " " << diam1 << " " << diam2 << " "
                   << ang1 << " " << ang2 << " " << dss << " " << dnn << std::endl;
              
               file9 << "darc" << std::endl;
               file9 << xcen << " " <<  ycen << " " <<  diam1 << " " <<  diam2 << " " << ang1 << " " <<  ang2 <<
                   " " << dss << " " << dnn << std::endl;
           }
           catch (std::ifstream::failure e)
           {
               MessageBox(nullptr, L"Error in input File,darc definition!", L"Error", MB_OK);

           }         
           insituS.incres = 3;
     
       return;
   }





   void  processMonitoring_points()
   {
       string lineData;
       MonitoringPoint mpoint;
       
       try
       {
           getline(inFile, lineData);
           std::stringstream ss(lineData);
           ss >> mpoint.xmon >> mpoint.ymon;
           mpoint_list[ihist] = mpoint;
           ihist++;
           if (ihist > 19)
           {
               file2 << "ERROR: Too many monitoring points (max=19) No monitoring point set for "
                   << mpoint.xmon << " " << mpoint.ymon << std::endl;
               ihist--;
               return;     
           }
           file2 << "history file = hist" << ihist << ".dat" << std::endl
               << "monitoring point = " << mpoint.xmon << " " << mpoint.ymon << std::endl;
           wstring filename =  monit_dir + L"/hist" + std::to_wstring(ihist) + L".dat";
            mon_files[ihist].open(filename);
           if (!mon_files[ihist]) {
               MessageBox(nullptr, L"Error in creating ihist File!", L"Error", MB_OK);
               }
           

           mon_files[ihist] << " -- Monitoring point close to a boundary will be calculated at the boundary and marked as AT BOUNDARY ELEMENT\n"
               << " -- in this case, the boundary normal and shear stress will be provided by sigxx, sigyy, and sigxy\n"
               << " -- if the boundary is vertical or horizontal, the normal stress = sigxx or sigyy, shear stress = sigxy\n"
               << " -- if the boundary is inclined, the normal stress = sigxx = sigyy, shear stress = sigxy\n\n"
               << "  Cycle(iteration)     Time         xp       yp        sigxx         sigyy       sigxy         dx          dy     max creep velocity\n";
       }
       catch (std::ifstream::failure e)
       {
           MessageBox(nullptr, L"Error in input File: monitoring point definition!", L"Error", MB_OK);
       }
       return;       

   }



   void  processMonitoring_lines()
   {
       MonitoringLine monlin;
       string lineData;
       try
       {
           getline(inFile, lineData);
           std::stringstream ss(lineData);
           ss >> monlin.x1l >> monlin.y1l >> monlin.x2l >> monlin.y2l >> monlin.npl;
           mline_list[lhist] = monlin;
           lhist++;

           if (lhist > 9)
           {
               file2 << "ERROR: Too many monitoring lines (max=9). No monitoring line set for "
                   << monlin.x1l << " " << monlin.y1l << " " << monlin.x2l << " " << monlin.y2l << std::endl;
               lhist--;
               return;
           }
           file2 << "history file = hist_line" << lhist << ".dat" << std::endl
               << "monitoring line = " << monlin.npl << " " << monlin.x1l << " " << monlin.y1l << " \n"
               << monlin.x2l << " " << monlin.y2l << std::endl;

           wstring filename = monit_dir + L"/hist_line" + std::to_wstring(lhist) + L".dat";
           ml_files[lhist].open(filename);
            if (!ml_files[lhist]) {
            MessageBox(nullptr, L"Error in creating hist_line File!", L"Error", MB_OK);
                }
            ml_files[lhist] << " -- Monitoring points in a line will ignore the existence of boundary elements or fractures\n"
               << " -- User should check if any points are too close to the existing elements\n"
               << " -- incorrect results could be resulted at these points\n\n"
               << "  Cycle(iteration)     Time     Line     xp       yp        sigxx         sigyy       sigxy         dx          dy              max creep velocity\n";
       
       }
       catch (std::ifstream::failure e)
       {
           MessageBox(nullptr, L"Error in input File: monitoring line definition!", L"Error", MB_OK);

       }
       return;

   }




   void processCreepParameters()
   {
       
       string lineData;
       try
       {
           getline(inFile, lineData);
           std::stringstream ss(lineData); 
           ss >> creep.v1 >> creep.nn1 >> creep.v2 >> creep.nn2 >>
           creep.totalT >> creep.deltaT_min >> creep.deltaT_max;
           int ID_creep = 1;   //not sure Sara!
           if (creep.totalT < 20.0)
           {
               creep.totalT = 20;
           }
           file2 << "Mode I creep velocity parameter Vmax = " << creep.v1 << "  n = " << creep.nn1 << "\n"
               << "Mode II creep velocity parameter Vmax = " << creep.v2 << "  n = " << creep.nn2 << "\n"
               << "Total creep time (seconds) = " << creep.totalT << "\n"
               << "Minimum time step (seconds) = " << creep.deltaT_min << "\n"
               << "Maximum time step (seconds) = " << creep.deltaT_max << std::endl;
       }
       catch (std::ifstream::failure e)
       {
           MessageBox(nullptr, L"Error in input File,creep definition!", L"Error", MB_OK);
       }
       return;
   }





   void processInsituStress()
   {
       float pxx0, pyy0, pxy0;  
       string lineData;
       try
       {
           getline(inFile, lineData);
           std::stringstream ss(lineData);
           ss >> pxx0 >> pyy0 >> pxy0;
           insituS.dsxx = pxx0; 
           insituS.dsxy = pxy0; 
           insituS.dsyy = pyy0;
           insituS.incres = 1;

           file2 << "In situ stresses:\n"
               << " xx-component of field stress = " << scientific << pxx0 << "\n"
               << " yy-component of field stress = " << scientific << pyy0 << "\n"
               << " xy-component of field stress = " << scientific << pxy0 << std::endl<<endl;
       }
       catch (std::ifstream::failure e)
       {
           MessageBox(nullptr, L"Error in input File,stress definition!", L"Error", MB_OK);

       }
       return;
   }






   void processGravity()
   {
       float dens_rock, gy, sh_sv_ratio; // y_surf - surface height
       string lineData;
       try
       {
           getline(inFile, lineData);
           std::stringstream ss(lineData);
           ss >> dens_rock >> gy >> sh_sv_ratio >> g.y_surf;
           g.sky = dens_rock * gy;
           g.skx = g.sky * sh_sv_ratio;
           file2 << "Gravitational & tectonic stress:\n"
               << " Rock density = " << dens_rock << "\n"
               << " Gravity acceleration = " << gy << "\n"
               << " Horizontal/vertical stress ratio = " << sh_sv_ratio << "\n"
               << " Ground surface height (y) = " << g.y_surf << std::endl;
       }
       catch (std::ifstream::failure e)
       {
           MessageBox(nullptr, L"Error in input File,gravity definition!", L"Error", MB_OK);

       }
       return;

   }




   void process_ContactSurface_Properties()
   {
       string lineData;
       getline(inFile, lineData);
       int jmat = 1;
       
       std::stringstream ss(lineData);

       if (ss)
       {
           ss >> jmat >> s8[jmat].aks0 >> s8[jmat].akn0 >> s8[jmat].phi0 >> s8[jmat].coh0 >>
               s8[jmat].phid0 >> s8[jmat].apert0 >> s8[jmat].apert_r;
       }
      else
      {
                   s8[jmat].phid0 = 0;
                   s8[jmat].apert0 = 1e-6;
                   s8[jmat].apert_r = 1e-6;
      }
      //Sara still not sure about this
      if (s8[jmat].apert0 == 0 )      
          s8[jmat].apert0 = 1e-6;

      if (s8[jmat].apert_r == 0)
          s8[jmat].apert_r = 1e-6;
      

           file2 << "   joint material = " << jmat << std::endl
            << "               ks = " << std::scientific << std::setprecision(4) << s8[jmat].aks0 << std::endl
            << "               kn = " << std::scientific << std::setprecision(4) << s8[jmat].akn0 << std::endl
            << "              phi = " << std::fixed << std::setprecision(1) << s8[jmat].phi0 << std::endl
            << "        cohesion = " << std::scientific << std::setprecision(4) << s8[jmat].coh0 << std::endl
            << "   dilation angle = " << std::fixed << std::setprecision(1) << s8[jmat].phid0 << std::endl
            << " initial aperture = " << std::scientific << std::setprecision(4) << s8[jmat].apert0 << std::endl
            << "residual aperture = " << std::scientific << std::setprecision(4) << s8[jmat].apert_r << std::endl;
                         
          
           s8[jmat].phi0 = s8[jmat].phi0 / 180.0 * pi;
           s8[jmat].phid0 = s8[jmat].phid0 / 180.0 * pi;
           j_material = jmat;
            
       return;       
   }




   void processRockStrength_properties()
   {
       irock = 1;
       int material = 1;
       float rphi_tem =0, rcoh_tem = 0, rst_tem = 0;
       string lineData;
       try
       {
           getline(inFile, lineData);
           std::stringstream ss(lineData);
           ss >> rphi_tem >> rcoh_tem >> rst_tem >> material;
           int mm = material;
           rock1[mm].rphi = rphi_tem;
           rock1[mm].rcoh = rcoh_tem;
           rock1[mm].rst = rst_tem;

           file2 << "Intact rock strength: phi = " << rphi_tem << "\n"
               << "                      coh = " << std::scientific << std::fixed <<rcoh_tem << "\n"
               << "                     sigt = " << std::scientific << std::fixed << rst_tem << std::endl;
           rock1[mm].rphi = rock1[mm].rphi / 180.0 * pi;
       }
       catch (std::ifstream::failure e)
       {
           MessageBox(nullptr, L"Error in input File,rock definition!", L"Error", MB_OK);

       }
       return;
   }



   
   void processRectWindow()
   {
       dispwin.ID_win = 1;
       string lineData;
       try
       {
           getline(inFile, lineData);
           std::stringstream ss(lineData);
           ss >> dispwin.xll >> dispwin.xur >> dispwin.yll >> dispwin.yur >>
               dispwin.numx >> dispwin.numy;               
          
           file2 << "Plot window (rect):  x = " << dispwin.xll << " " << dispwin.xur << "\n"
               << "                     y = " << dispwin.yll << " " << dispwin.yur << "\n"
               << "  Grid number: x and y = " << dispwin.numx << " " << dispwin.numy << std::endl<<"\n"<< "\n";
       }
       catch (std::ifstream::failure e)
       {
           MessageBox(nullptr, L"Error in input File,window definition!", L"Error", MB_OK);;

       }
       return;
       
   }



   void  processCircWindow()
   {
       dispwin.ID_win = 2;
       string lineData;
       try
       {
           getline(inFile, lineData);
           std::stringstream ss(lineData);
            ss >> dispwin.xc0 >> dispwin.yc0 >> dispwin.radium >> dispwin.numr >> dispwin.numa;
           file2 << "Plot window (circle):  centre = " << dispwin.xc0 << " " << dispwin.yc0 << "\n"
               << "                       radium = " << dispwin.radium << "\n"
               << "       rid number: r and seta = " << dispwin.numr << " " << dispwin.numa << std::endl;

           dispwin.xll = dispwin.xc0 - dispwin.radium;
           dispwin.xur = dispwin.xc0 + dispwin.radium;
           dispwin.yll = dispwin.yc0 - dispwin.radium;
           dispwin.yur = dispwin.yc0 + dispwin.radium;
       }
       catch (std::ifstream::failure e)
       {
           MessageBox(nullptr, L"Error in input File,window definition!", L"Error", MB_OK);

       }
       return;
   }




   void processFracinitWindow()
   {
       int iwin = 1;
       string lineData;
       try
       {
           getline(inFile, lineData);
           std::stringstream ss(lineData);
            ss >> s5u.xmin >> s5u.xmax >> s5u.ymin >> s5u.ymax;

           file2 << "Fracture initiation window:  X = " << s5u.xmin << " " << s5u.xmax << "\n"
               << "                             Y = " << s5u.ymin << " " << s5u.ymax << std::endl;

           // Check if window parameters are not specified
           if (iwin == 0 && s5u.xmin == 0 && s5u.xmax == 0 && s5u.ymin == 0 && s5u.ymax == 0)
           {
               s5u.xmin = -100000.0;
               s5u.xmax = 100000.0;
               s5u.ymin = -100000.0;
               s5u.ymax = 100000.0;
           }
       }
       catch (std::ifstream::failure e)
       {
           MessageBox(nullptr, L"Error in input File,window definition!", L"Error", MB_OK);

       }
       return;
   }




   void processTunnel()
   {
       float xcen, ycen, diam;
       try
       {
           if (!(inFile >> xcen >> ycen >> diam))
           {
               MessageBox(nullptr, L"Error in input File,tunnel definition!", L"Error", MB_OK);
           }
           
           tunnl.diameter[ntunnel] = diam;  //need to optimize Sara!
           tunnl.x_cent[ntunnel] = xcen;
           tunnl.y_cent[ntunnel] = ycen;
           ntunnel++;
       }
       catch (std::ifstream::failure e)
       {
           std::cerr << "Exception opening/reading/closing file: in processTunnel\n";

       }
       return;
      
   }



   void processArch()
   {
       int material = 1;
       float xcen = 0.0, ycen = 0.0, diam = 0.0, ang1 = 0.0, ang2 = 0.0,
           bvs = 0.0, bvn = 0.0,
           gradsy = 0.0,
           gradny = 0.0;
       int num, kode;
       string lineData;
       
      
       try
       {
           getline(inFile, lineData);
           //removeCommas(lineData);
           std::stringstream ss(lineData);
           ss >> num >> xcen >> ycen >> diam >> ang1 >> ang2 >> kode >> bvs >> bvn >>
               material >> gradsy >> gradny;
           float angle = 0;
           if (symm.ksym == 0)
               angle = 360.0;
           else if (symm.ksym == 1 || symm.ksym == 2 || symm.ksym == 3)
               angle = 180.0;
           else if (symm.ksym == 4)
               angle = 90.0;

           if (abs(ang2 - ang1) == angle)
           {
               
               tunnl.diameter[ntunnel] = diam;  //need to optimize Sara!
               tunnl.x_cent[ntunnel] = xcen;
               tunnl.y_cent[ntunnel] = ycen;
               ntunnel++;
           }

           file2 << "ARC ---- centre position: x,y  = " << xcen << " " << ycen << "\n"
               << "         diameter              = " << diam << "\n"
               << "         start angle           = " << ang1 << "\n"
               << "         end angle             = " << ang2 << "\n"
               << "         shear stress or disp  = " << bvs << "\n"
               << "         normal stress or disp = " << bvn << "\n"
               << "         shear value gradient in y direction = " << gradsy << "\n"
               << "         normal value gradient in y direction = " << gradny << "\n"
               << "         number of elements    = " << num << "\n"
               << "         material              = " << material << std::endl;

          

           //define a new arc and initialize it with user data

           float ang11 = ang1 * pi / 180.0;
           float ang12 = ang2 * pi / 180.0;

           switch (kode)
           {
           case 1:
           case 11:
               arc_list[na] = Arch(xcen, ycen, diam / 2.0, ang11, ang12, bvs, bvn, 0, 0, gradsy,
                   gradny, num, material, kode);   
               break;

           case 2:
           case 12:
               arc_list[na] = Arch(xcen, ycen, diam / 2.0, ang11, ang12, 0, 0, bvs, bvn, gradsy, 
                   gradny, num, material, kode);   
               break;

           case 3:
           case 13:
               arc_list[na] = Arch(xcen, ycen, diam / 2.0, ang11, ang12, 0, bvn, bvs, 0, gradsy,
                   gradny, num, material, kode);    
               break;

           case 4:
           case 14:
               arc_list[na] = Arch(xcen, ycen, diam / 2.0, ang11, ang12, bvs, 0, 0, bvn, gradsy,
                   gradny, num, material, kode);                 
           }
           na++;
       }
       catch (std::ifstream::failure e)
       {
           std::cerr << "Exception opening/reading/closing file: in processArch\n";

       }
       return;
   }




   void processEllipticalOpens()
   {
       int num = 0;
       int material = 1;
       float xcen, ycen, diam1, diam2, kode, bvs, bvn, gradsy, gradny;
       string lineData;
       try
       {
           getline(inFile, lineData);
           std::stringstream ss(lineData);
           ss >>num>> xcen >> ycen >> diam1 >> diam2 >> kode >> bvs >> bvn >> material >> gradsy >> gradny;                    
           Elliptical_opening ep(xcen, ycen, diam1, diam2);
           
           ellip_list[nellipse] = ep;    
          
           file2 << "ARC ---- centre position: x,y  = " << xcen << " " << ycen << "\n"
               << "         diameter a              = " << diam1 << "\n"
               << "         diameter b              = " << diam2 << "\n"
               << "         shear stress or disp  = " << bvs << "\n"
               << "         normal stress or disp = " << bvn << "\n"
               << "         shear value gradient in y direction = " << gradsy << "\n"
               << "         normal value gradient in y direction = " << gradny << "\n"
               << "         number of elements    = " << num << "\n"
               << "         material              = " << material << std::endl;

           float seta1 = 0.0, seta2 = 0.0;

           switch (symm.ksym)
           {
           case 0:
               seta1 = 0;
               seta2 = 2.0 * pi; 
               break;
           case 1:
               seta1 = 0;
               seta2 = pi;
               break;
           case 2:               
           case 3:
               seta1 = 0;
               seta2 = pi;
               break;
           case 4:
               seta1 = 0;
               seta2 = pi/2.0;
           }
           float dseta = (seta2 - seta1) / num;

           for (int m = 0; m < num; ++m)
           {
               float seta_b = seta1 + m * dseta;
               float seta_e = seta1 + (m + 1) * dseta;
               float xbeg = xcen + 0.5 * diam1 * cosf(seta_b);
               float ybeg = ycen + 0.5 * diam2 * sinf(seta_b);
               float xend = xcen + 0.5 * diam1 * cosf(seta_e);
               float yend = ycen + 0.5 * diam2 * sinf(seta_e);
               
               int jmat = 1; // jmat here is meaningless
               int itype = 0; // meaningless
               ellip_list[nellipse].def_boundary_elements_for_Geoform(1,xbeg, ybeg,
                   xend, yend, bvs, bvn, gradsy, gradny,itype, jmat);
           }
       }
       catch (std::ifstream::failure e)
       {
           std::cerr << "Exception reading file: in processEllipticalOpens\n";

       }
       nellipse++;
       return;

   }






  inline  void settProcess()
   {
      string lineData;
      try
      {
          getline(inFile, lineData);
          std::stringstream ss(lineData);
          ss >> factors.tolerance;
          file2 << "Fracture tip merging tolerance distance is set to be " <<
              factors.tolerance << std::endl;
      }
      catch (std::ifstream::failure e)
      {
          std::cerr << "Exception opening/reading/closing file: in processEllipticalOpens\n";

      }
      return;

   }


  inline void seteProcess() {
       inFile >> factors.factor_e;
       file2 << "Elastic fracture propagation checking level is set to be " <<
           factors.factor_e << "Kc" << std::endl;
   }

  inline void setfProcess() {
       inFile >> factors.factor_f;
       file2 << "Fracture initiation cut-off level is set to be " << factors.factor_f
           << "% of max FoS" << std::endl;
   }


   void randProcess() {
       string linedata;
       s15.i_rand = 1, s15.l_rand = 1;

       try
       {
           getline(inFile, linedata);
           std::stringstream ss(linedata);
           ss >> s15.f_ini0 >> s15.l_rand;
           file2 << "Fracture initiation starts at " << s15.f_ini0 << " strength"
               << std::endl;
           file2 << "Random level (1 - 100%; 0 - no)" << s15.l_rand << std::endl;
       }
       catch (std::ifstream::failure e)
       {
           std::cerr << "Exception opening/reading/closing file:in input rand\n";

       }
   }
   inline void titl_Process() {
       string tem;
       getline(inFile, tem);
       file2 << tem << std::endl;
   }

   inline void boundProcess()
   {
       s15.i_bound = 1;
       file2 << "Fracture initiation at boundaries is allowed" << std::endl;
   }

   inline void inteProcess()
   {
       s15.i_intern = 1;
       file2 << "Fracture initiation in intact rock is allowed" << std::endl;
   }

   inline void isizProcess()
   {
       inFile >> s15.a_ini;
       file2 << "Fracture initiation element size is set to be " <<
           std::setprecision(3) << s15.a_ini << std::endl;
   }




   inline void iterProcess() {
       inFile >> n_it;
       file2 << "Iteration cycle number = " << n_it << std::endl;
   }




   inline void numeProcess()
   {
       string lineData;
       try
       {
           getline(inFile, lineData);
           std::stringstream ss(lineData);
           ss >> k_num >> d_max;
           file2 << "Set numerical stability parameters\n"
               << " boundary element k_num = " << k_num << "\n"
               << " Max fracture disp = " << d_max << std::endl;
       }
       catch (std::ifstream::failure e)
       {
           std::cerr << "Exception opening/reading/closing file\n";

       }
   }



   inline void  mlinProcess()
   {
       try
       {
           inFile >> mat_lining;
           file2 << "MATERIAL NO = " << mat_lining <<
               " is set to be the concrete lining without insitu stresses" << std::endl;
       }
       catch (std::ifstream::failure e)
       {
           std::cerr << "Exception opening/reading/closing file\n";

       }
   }



   inline void  dstrProcess()
   {
       string lineData;
       try
       {
           getline(inFile, lineData);
           std::stringstream ss(lineData);
           ss >> insituS.dsxx >> insituS.dsyy >> insituS.dsxy;

           file2 << "Insitu stress increment, dsxx, dsyy, dsxy = " << insituS.dsxx << " " <<
               insituS.dsyy << " " << insituS.dsxy << std::endl;
           insituS.incres = 1;
       }
       catch (std::ifstream::failure e)
       {
           std::cerr << "Exception opening/reading/closing file in dstr\n";

       }
   }







   void input()
   {
      
       string id;
       string tem;
       string message;
       string lineData;
      
      
       //auto start = std::chrono::high_resolution_clock::now();
       using func = function<void()>;     

       unordered_map<string, func> mp{

        {"titl", titl_Process},
        {"symm", symmetryProcessing},
        {"toug",frac_energy_processing},
        {"touk",frac_toughness_process},
        {"modu", elastic_modul_processing},
        {"stre" ,processInsituStress},
        {"swin",processRectWindow },
        {"arch",processArch},
        {"frac",processFracture},
        {"cycl",start_calculation},
        {"rand", randProcess},
        {"endf",endFile},
        {"iter",iterProcess },
        {"prop",process_ContactSurface_Properties},
        {"rock",processRockStrength_properties},
        {"inte",inteProcess},
        {"exca",processExcavation_Cracks},
        {"dbou",dbouProcess},
        {"darc",darcProcess},
        {"save",saveData},
        {"elli",processEllipticalOpens},{ "edge",processEdge},{"gost",processGost},
        {"moni",processMonitoring_points},
        {"monl",processMonitoring_lines},
        {"nume",numeProcess},
        {"mlin",mlinProcess},
        {"cree",processCreepParameters},
        {"setf",setfProcess},
        {"sete",seteProcess},
        {"sett",settProcess},
        {"boun",boundProcess},
        {"isiz",isizProcess},
        {"grav",processGravity},
        {"dstr",dstrProcess},
        {"rwin",processCircWindow},
        {"iwin",processFracinitWindow},
        {"tunn",processTunnel},
        {"lint",processLinearInterface},{ "aint",processArcInterface},{"wate",processWaterPressure},
        {"perm",processPermeability},
        {"rect",processRect},
        {"hole",processHole}
       };  



       try
       {
           while (std::getline(inFile, lineData))   
           {
               if (lineData.empty() || lineData=="\t")
               {
                   continue;
               }
               else 
               if (lineData.front() == '*')
               {
                   std::getline(inFile, lineData); // Read and discard the next line
                   continue;
               }
               
               line++;
               if (lastinput == "endf")
                   return;
               id = lineData.substr(0, 4);
               ToLowerCase(id);
               line++;
               lastinput = id;
               auto it = mp.find(id);
               if (it != mp.end()) {
                   // Call the corresponding function
                   it->second();
                   if (id == "cycl")
                       return;
               }
               else {
                   
                   std::wstring fullMessage = L"Error in input file: line " + std::to_wstring(line-1) + L" .Do you want to continue?";

                   int response = MessageBox(nullptr, fullMessage.c_str(), L"Confirmation", MB_OKCANCEL | MB_ICONQUESTION);
                   if (response == IDOK) {
                       continue;
                   }
                   else if (response == IDCANCEL) {
                       logfile << "exits with: " << IDCANCEL << endl;
                       exit(EXIT_FAILURE);
                   }

               }
           }
           if (inFile.eof()) {
               lastinput = "endf";
               return;
           }
       }
       catch (std::ifstream::failure e)
       {
           std::cerr << "Exception opening/reading/closing file in input\n";
           exit(EXIT_FAILURE);

       }
       
       return;
          
   }