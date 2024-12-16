#include <common.h>
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


using namespace CommonPara_h::comvar;



std::ifstream inFile;

void processEdge()
{
       int material = 1;
       float gradsy = 0, gradny = 0;        
       int num, kode;
       float xbeg, ybeg, xend, yend, bvs, bvn;
                
       try
       {
           
           inFile >> num >> xbeg >> ybeg >> xend >> yend >> kode >> bvs >> bvn >> material >> gradsy >> gradny;
           file2 << "Straigt boundary-- from: x,y = " << xbeg << "," << ybeg << std::endl;
           file2 << "                    to        = " << xend << "," << yend << std::endl;
           file2 << "      shear stress or disp    = " << bvs << std::endl;
           file2 << "       normal stress or disp  = " << bvn << std::endl;
           file2 << "      shear value gradient in y direction = " << gradsy << std::endl;
           file2 << "      normal value gradient in y direction = " << gradny << std::endl;
           file2 << "         number of elements   = " << num << std::endl;
           file2 << "         material              = " << material << std::endl;
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
        file2 << "Gost element        -- from: x,y = " << xbeg << "," << ybeg << std::endl;
        file2 << "                    to        = " << xend << "," << yend << std::endl;
        file2 << "      shear stress or disp    = " << bvs << std::endl;
        file2 << "       normal stress or disp  = " << bvn << std::endl;
        file2 << "      shear value gradient in y direction = " << gradsy << std::endl;
        file2 << "      normal value gradient in y direction = " << gradny << std::endl;
        file2 << "         number of elements   = " << num << std::endl;
        file2 << "         material              = " << material << std::endl;
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
    try
    {
        inFile >> num >> xbeg >> ybeg >> xend >> yend >> itype >> jmat >> material;       
        file2 << "Fracture ------ from: x, y =  " << xbeg << "," << ybeg << "," <<
            "     to     =  " << xend << "," << yend << "\n" <<
            "      number of elements   = " << num << "\n" << "        material     =   " << material << "\n";
        
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
        try
        {    
            inFile >> num >> xbeg >> ybeg >> xend >> yend >> mat1 >> mat2;            
            file2 << "Interface -------- from: x, y =  " << xbeg << "," << ybeg << "\n" <<
                "                    to      =  " << xend << "," << yend << "\n" <<
                "         number of elements   = " << num << "\n";
            
            Edge_interface new_lin_intfce(xbeg, ybeg, xend, yend, mat1, mat2,num); 
            lin_intrfce_list[npli] = new_lin_intfce;
            npli++;
        }
        catch (std::ifstream::failure e)
        {
            MessageBox(nullptr, L"Error in input File,linear interface definition!", L"Error", MB_OK);

        }
    }



    void processArcInterface()
    {
        try
        {
            int num = 0, mat1 = 0, mat2 = 0;
            float xcen = 0, ycen = 0, diam = 0, ang1 = 0, ang2 = 0, bvs = 0, bvn = 0;
            inFile >> num >> xcen >> ycen >> diam >> ang1 >> ang2 >> mat1 >> mat2;           
            file2 << "ARC Int- centre position: x,y  = " << xcen << "," << ycen << "\n" <<
                "         diameter              = " << diam << "\n" <<
                "         start angle           = " << ang1 << "\n" <<
                "         end angle             = " << ang2 << "\n" <<
                "         shear stress or disp  = " << bvs << "\n" <<
                "         normal stress or disp = " << bvn << "\n" <<
                "         number of elements     = " << num << "\n" <<
                "         material one side     = " << mat1 << "\n" <<
                "         material other side   = " << mat2 << "\n";

            
            Arch_interface newArcIntrfce(xcen, ycen, diam / 2., ang1 * pi / 180.,
                ang2 * pi / 180., num, mat1, mat2);
            arc_intrface_list[nq] = newArcIntrfce;
            int kode = 6;   
            bvs = 0;
            bvn = 0;
            nq = nq + 1;
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
        try
        {
            getline(inFile, lineData);
            std::stringstream ss(lineData);

            ss >> id ;           
            ToLowerCase(id);
            id = id.substr(0, 4);
            if (id == "hole")
            {
                
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
            if (id == "rect")
            {
                
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




    void processPermeability()
    {
    
        try
        {
            inFile >> perm.viscosity >> perm.density >> perm.perm0;
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
            MessageBox(nullptr, L"Error in input File,touc definition!", L"Error", MB_OK);
        }
        
        if (tem != "    ") 
            {
                  file2 << tem << endl;              
                    mcyc0 += stoi(tem);               
            }
        //----------------------------
        if (mcyc != 0) return;            
        s5u.dtol = 0.0001;            
           
        if (nc == 1)
        {
            inputcheck();
            input_tip_check();
          //  file7 << dispwin.xll << " " << dispwin.xur << " " << dispwin.yll << " " <<
             //   dispwin.yur << endl;           

             //!-------------------------
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

            water();   
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

            float aks_bb = 0, akn_bb = 0, phi_bb = 0; float coh_bb = 0.0f;
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
            file2 << "New fresh fracture properties --- Shear fractures" << std::endl
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
       string id;
       string message;
       //try
       //{
       //    inFile >> id;

       //    // Error handling     not sure Sara!
       //    if (id == "ERR")
       //    {
       //        int errorCode;
       //        inFile >> errorCode;
       //        file2 << "Input data format error, line=" << line << "***" << id << "***" << endl;
       //        file2 << "Error code: " << errorCode << endl;
       //        // Read the rest of the line
       //        getline(inFile, message);
       //        //call SendWindowText(message//CHAR(0))
       //        outFile << message << endl;
       //        return;
       //    }
       //    else if (id == "INSF")
       //    {
       //        outFile << "Insufficient input data, line=" << line << "***" << id << "***" << endl;
       //        // Read the rest of the line
       //        getline(inFile, message);
       //        //call SendWindowText(message//CHAR(0))
       //        file2 << message << endl;
       //        return;
       //    }
       //    // End of input file
       //    //call SendWindowText('End of input file'//CHAR(0))
       //    lastinput = "endf";
       //    // file2 << "End of input file" << endl;}

       //}
       //catch (std::ifstream::failure e)
       //{
       //    std::cerr << "Exception opening/reading/closing file:in endFile func\n";

       //}

       return;
   }




   void processExcavation_Cracks()
   {
       exca.ID_Exca = 1;
       try
       {
           inFile >> exca.d_wall >> exca.rand_e;
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
       std::ofstream file9("bound.dat");
       //---------------straight boundary value increment ------------------
      
           float x1, x2, y1, y2, dss, dnn;
           try
           {
               inFile >> x1 >> x2 >> y1 >> y2 >> dss >> dnn;               
               file2 << "boundary value increment, x1,x2,y1,y2,dss,dnn =" << std::endl
                   << x1 << x2 << y1 << y2 << dss << dnn;
               file9 << "('dbou')";
               file9 << x1 << x2 << y1 << y2 << dss << dnn;   
              
           }
           catch (std::ifstream::failure e) 
           {
               MessageBox(nullptr, L"Error in input File,dboundary definition!", L"Error", MB_OK);

           }
           insituS.incres = 3;    
           
           return;
   }


   void darcProcess()
   {
       //---------------arch boundary value increment ------------------
      
           float xcen, ycen, diam1, diam2, ang1, ang2, dss, dnn;
           try
           {
               inFile >> xcen >> ycen >> diam1 >> diam2 >> ang1 >> ang2 >> dss >> dnn;
               file2 << "Arch boundary value increment, xc, xy, d1, d2, ang1, ang2, dss, dnn =\n"
                   << xcen << " " << ycen << " " << diam1 << " " << diam2 << " "
                   << ang1 << " " << ang2 << " " << dss << " " << dnn << std::endl;
               std::ofstream file9("darc.dat");

               file9 << xcen << " " << ycen << " " << diam1 << " " << diam2 << " "
                   << ang1 << " " << ang2 << " " << dss << " " << dnn << std::endl;
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
       getline(inFile, lineData);
       ihist++;
       MonitoringPoint mpoint;    
       std::stringstream ss(lineData);

       try
       {
           ss >> mpoint.xmon >> mpoint.ymon;
           if (ihist > 19)
           {
               file2 << "ERROR: Too many monitoring points (max=19) No monitoring point set for "
                   << mpoint.xmon << " " << mpoint.ymon << std::endl;
               ihist--;
               return;     
           }
           file2 << "history file = hist" << ihist << ".dat" << std::endl
               << "monitoring point = " << mpoint.xmon << " " << mpoint.ymon << std::endl;
           mpoint_list.push_back(mpoint);     

           std::ofstream histPFile("hist" + std::to_string(ihist) + ".dat");           
           histPFile << " -- Monitoring point close to a boundary will be calculated at the boundary and marked as AT BOUNDARY ELEMENT\n"
               << " -- in this case, the boundary normal and shear stress will be provided by sigxx, sigyy, and sigxy\n"
               << " -- if the boundary is vertical or horizontal, the normal stress = sigxx or sigyy, shear stress = sigxy\n"
               << " -- if the boundary is inclined, the normal stress = sigxx = sigyy, shear stress = sigxy\n\n"
               << "  Cycle(iteration)     Time         xp       yp        sigxx         sigyy       sigxy         dx          dy     max creep velocity\n";
       }
       catch (std::ifstream::failure e)
       {
           MessageBox(nullptr, L"Error in input File,monitoring point definition!", L"Error", MB_OK);
       }
       return;       

   }



   void  processMonitoring_lines()
   {
       lhist++;
       MonitoringLine monlin;
       try
       {
           inFile >> monlin.x1l >> monlin.y1l >> monlin.x2l >> monlin.y2l >> monlin.npl;

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

           std::ofstream histLFile("hist_line" + std::to_string(lhist) + ".dat");
           if (!histLFile) {
               // Handle file open error
           }
           histLFile << " -- Monitoring points in a line will ignore the existence of boundary elements or fractures\n"
               << " -- User should check if any points are too close to the existing elements\n"
               << " -- incorrect results could be resulted at these points\n\n"
               << "  Cycle(iteration)     Time     Line     xp       yp        sigxx         sigyy       sigxy         dx          dy              max creep velocity\n";
       }
       catch (std::ifstream::failure e)
       {
           MessageBox(nullptr, L"Error in input File,monitoring line definition!", L"Error", MB_OK);

       }
       return;

   }




   void processCreepParameters()
   {
       try
       {           
           inFile >> creep.v1 >> creep.nn1 >> creep.v2 >> creep.nn2 >>
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
       try
       {
           inFile >> pxx0 >> pyy0 >> pxy0;

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
       try
       {
           inFile >> dens_rock >> gy >> sh_sv_ratio >> g.y_surf;
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
            
       return;       
   }




   void processRockStrength_properties()
   {
       irock = 1;
       int material = 1;
       float rphi_tem, rcoh_tem, rst_tem;
       try
       {
           inFile >> rphi_tem >> rcoh_tem >> rst_tem >> material;
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
       try
       {
           inFile >> dispwin.xll >> dispwin.xur >> dispwin.yll >> dispwin.yur >>
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
       try
       {
           inFile >> dispwin.xc0 >> dispwin.yc0 >> dispwin.radium >> dispwin.numr >> dispwin.numa;
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
       try
       {
           inFile >> s5u.xmin >> s5u.xmax >> s5u.ymin >> s5u.ymax;

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
       getline(inFile,lineData);
       std::stringstream ss(lineData);
      
       try
       {
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
       try
       {
           inFile >> xcen >> ycen >> diam1 >> diam2 >> kode >> bvs >> bvn >> material >> gradsy >> gradny;
               
           nellipse++;
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
               seta1 = 0.0;
               seta2 = 2.0 * pi;
           case 2:
               seta1 = -pi / 2.0;
               seta2 = pi / 2.0;
           case 3:
               seta1 = 0.0;
               seta2 = pi;
           case 4:
               seta1 = 0.0;
               seta2 = pi;
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

               int num0 = 1;
               int jmat = 1; // jmat here is meaningless
               int itype = 0; // meaningless
               ellip_list[nellipse].def_boundary_elements_for_Geoform(num0,xbeg, ybeg,
                   xend, yend, bvs, bvn, gradsy, gradny,
                   itype, jmat);
           }
       }
       catch (std::ifstream::failure e)
       {
           std::cerr << "Exception opening/reading/closing file: in processEllipticalOpens\n";

       }
       return;

   }






  inline  void settProcess()
   {
       inFile >> factors.tolerance;
       file2 << "Fracture tip merging tolerance distance is set to be " <<
           factors.tolerance << std::endl;

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
           scientific << std::setprecision(3) << s15.a_ini << std::endl;
   }




   inline void iterProcess() {
       inFile >> n_it;
       file2 << "Iteration cycle number = " << n_it << std::endl;
   }




   inline void numeProcess()
   {
       try
       {
           if (!(inFile >> k_num >> d_max)) {
              
           }
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
       try
       {

           inFile >> insituS.dsxx >> insituS.dsyy >> insituS.dsxy;

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
      
       string strPathFile = filepath + "/Example" +
           to_string(test_id) + ".dat";
       inFile.open(strPathFile.c_str(), std::ios_base::in);
       if (!inFile.is_open())
       {
           MessageBox(nullptr, L"File path is wrong!", L"Error", MB_OK);
           exit(EXIT_FAILURE);
       }
      
       auto start = std::chrono::high_resolution_clock::now();
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
        { "dbou",dbouProcess},
        {"darc",darcProcess},
        {"save",saveData},
        {"elli",processEllipticalOpens},{ "edge",processEdge},{"gost",processEdge},
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
        {"perm",processPermeability}
       };       

       try
       {
           while (std::getline(inFile, lineData))   
           {
               if (lineData.empty())			
               {
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
               }
               else {
                   MessageBox(nullptr, L"Error in input file.", L"Error", MB_OK);
                  exit(EXIT_FAILURE);
               }
           }
       }
       catch (std::ifstream::failure e)
       {
           std::cerr << "Exception opening/reading/closing file in input\n";
           exit(EXIT_FAILURE);

       }
       
       return;
          
   }