#include <common.h>
#include <CommonPara.h>
#include <WinInterface.h>
#include<Work.h>
#include<ExcavationCracks.h>
#include<InputCheck.h>
#include<iostream>

using namespace CommonPara_h::comvar;


void processEdge(ifstream& inFile, ofstream& outFile,char id,string type)
{
       int material = 1;
       float gradsy = 0, gradny = 0;        
       int num, kode;
       float xbeg, ybeg, xend, yend, bvs, bvn;
                
       try
       {
           
           inFile >> num >> xbeg >> ybeg >> xend >> yend >> kode >> bvs >> bvn >> material >> gradsy >> gradny;

           //int mm = material;

           outFile << type << "        -- from: x,y = " << xbeg << "," << ybeg << std::endl;
           outFile << "                    to        = " << xend << "," << yend << std::endl;
           outFile << "      shear stress or disp    = " << bvs << std::endl;
           outFile << "       normal stress or disp  = " << bvn << std::endl;
           outFile << "      shear value gradient in y direction = " << gradsy << std::endl;
           outFile << "      normal value gradient in y direction = " << gradny << std::endl;
           outFile << "         number of elements   = " << num << std::endl;
           outFile << "         material              = " << material << std::endl;
       }
       catch (std::ifstream::failure e)
       {
           std::cerr << "Exception opening/reading/closing file:in processEdge\n";

       }

       
        if (id == 'S')
        {
            switch (kode)
            {
            case 1:
            case 11:
                bund_list[nb] = Edge(material, num, kode, xbeg, ybeg, xend, yend, bvs, bvn, 0, 0, gradsy, gradny);
                //bund_list[nb] = new_sBound1;   //sn[nb] = bvn;   ss[nb] = bvs;
                break;
            case 2:
            case 12:
                bund_list[nb] = Edge(material, num, kode, xbeg, ybeg, xend, yend, 0, 0, bvs, bvn, gradsy, gradny);
                //bund_list[nb] = new_sBound2;      //dn[nb] = bvn; ds[nb] = bvs;

                break;
            case 3:
            case 13:
                bund_list[nb] = Edge(material, num, kode, xbeg, ybeg, xend, yend, 0, bvn, bvs, 0, gradsy, gradny);
                //bund_list[nb] = new_sBound3;     // bbdn[nb] = bvn;     bbss[nb] = bvs;
                break;
            case 4:
            case 14:
                bund_list[nb] = Edge(material, num, kode, xbeg, ybeg, xend, yend, bvs, 0, 0, bvn, gradsy, gradny);
                //bund_list[nb] = new_sBound4;             //bbdn[nb] = bvn;       bbss[nb] = bvs;
                break;
            }
        }
        // for GostElement
        else
        {
            Edge new_sBound(material, num, 7, xbeg, ybeg, xend, yend, bvs, bvn, 0, 0, gradsy, gradny);
            bund_list[nb] = new_sBound;   //sn[nb] = bvn;   ss[nb] = bvs;
        }
        nb++;
        return;
}



void ToLowerCase(string&  str) {
   /* for (int i = 0; str[i]; i++) {
        str[i] = tolower(str[i]);
    }*/
    // Convert each character to lowercase
    for (char& c: str) 
    {
        c = std::tolower(c);
    }
}




void processFracture(ifstream& inFile, ofstream& outFile)
{
    int material = 1;
      
    int num, itype, jmat;
    float xbeg, ybeg, xend, yend;
    //Sara! there are also kod ,bvs and bvn to setting to 0 but not sure which vars 
    try
    {
        inFile >> num >> xbeg >> ybeg >> xend >> yend >> itype >> jmat >> material;       
        outFile << "Fracture ------ from: x, y =  " << xbeg << "," << ybeg << "," <<
            "     to     =  " << xend << "," << yend << "\n" <<
            "      number of elements   = " << num << "\n" << "        material     =   " << material << "\n";
        
        Fracture frac(xbeg, ybeg, xend, yend, material, num, jmat);         //dfine new fracture
        frac_list[nf] = frac;
        frac_list[nf].bound_type = 5;
        nf++;
    }
    catch (std::ifstream::failure e)
    {
        std::cerr << "Exception opening/reading/closing file:in processFracture\n";

    }   
     return;
}




 void processLinearInterface(ifstream & inFile, ofstream& outFile)
    {
        try
        {
            int num, mat1, mat2;
            float xbeg, ybeg, xend, yend, bvs = 0, bvn = 0;

            inFile >> num >> xbeg >> ybeg >> xend >> yend >> mat1 >> mat2;
            
            outFile << "Interface -------- from: x, y =  " << xbeg << "," << ybeg << "\n" <<
                "                    to      =  " << xend << "," << yend << "\n" <<
                "         number of elements   = " << num << "\n";
            
            Edge_interface new_lin_intfce(xbeg, ybeg, xend, yend, mat1, mat2,num);  //new linear intrface added
            lin_intrfce_list[npli] = new_lin_intfce;
            npli++;
        }
        catch (std::ifstream::failure e)
        {
            std::cerr << "Exception opening/reading/closing file:in processLinearInterface\n";

        }
    }



    void processArcInterface(ifstream& inFile, ofstream& outFile)
    {
        try
        {
            int num = 0, mat1 = 0, mat2 = 0;
            float xcen = 0, ycen = 0, diam = 0, ang1 = 0, ang2 = 0, bvs = 0, bvn = 0;

            inFile >> num >> xcen >> ycen >> diam >> ang1 >> ang2 >> mat1 >> mat2;
           
            outFile << "ARC Int- centre position: x,y  = " << xcen << "," << ycen << "\n" <<
                "         diameter              = " << diam << "\n" <<
                "         start angle           = " << ang1 << "\n" <<
                "         end angle             = " << ang2 << "\n" <<
                "         shear stress or disp  = " << bvs << "\n" <<
                "         normal stress or disp = " << bvn << "\n" <<
                "         number of elements     = " << num << "\n" <<
                "         material one side     = " << mat1 << "\n" <<
                "         material other side   = " << mat2 << "\n";

            
            //new arc interface was defined 
            Arch_interface newArcIntrfce(xcen, ycen, diam / 2., ang1 * pi / 180.,
                ang2 * pi / 180., num, mat1, mat2);
            arc_intrface_list[nq] = newArcIntrfce;
            int kode = 6;    //not sure why these paras considered Sara!
            bvs = 0;
            bvn = 0;
            nq = nq + 1;
        }
        catch (std::ifstream::failure e)
        {
            std::cerr << "Exception opening/reading/closing file:in processArcInterface\n";

        }
        return;                   
    
    }




    void processWaterPressure(ifstream& inFile, ofstream& outFile)
    {
        string id;
        //inFile.ignore(); // Skip the rest of the line
        string lineData;
        try
        {
            getline(inFile, lineData);
            std::stringstream ss(lineData);

            ss >> id ;
           // inFile >> id;
            // Convert id to lowercase
            ToLowerCase(id);
            id = id.substr(0, 4);

            if (id == "hole")
            {
                
                int iwhole = watercm.iwhole;     //alias for index
               // inFile.ignore();                   // Skip to next line
                ss >> watercm.w_xc[iwhole] >> watercm.w_yc[iwhole] >> watercm.w_d[iwhole]
                    >> watercm.wph[iwhole];

                watercm.ID_range = 1;
                watercm.iwhole++;
                outFile << "Water pressure in hole" << endl
                    << "     center x = " << watercm.w_xc[iwhole] << endl
                    << "     center y = " << watercm.w_yc[iwhole] << endl
                    << "     diameter = " << watercm.w_d[iwhole] << endl
                    << "     pressure = " << watercm.wph[iwhole] << endl << endl;
            }
            if (id == "rect")
            {
                
                int iwrect = watercm.iwrect;     //alias for index
                //inFile.ignore();                     // Skip to next line

                ss >> watercm.w_x1[iwrect] >> watercm.w_x2[iwrect] >> watercm.w_y1[iwrect] >>
                    watercm.w_y2[iwrect] >> watercm.wpr[iwrect];

                watercm.ID_range = 2;
                watercm.iwrect++;
                outFile << "Water pressure in rectangular opening" << endl
                    << "     x1 = " << watercm.w_x1[iwrect] << endl
                    << "     x2 = " << watercm.w_x2[iwrect] << endl
                    << "     y1 = " << watercm.w_y1[iwrect] << endl
                    << "     y2 = " << watercm.w_y2[iwrect] << endl
                    << "     pressure = " << watercm.wpr[iwrect] << endl << endl;
            }
        }
        catch (std::ifstream::failure e)
        {
            std::cerr << "Exception opening/reading/closing file:in processWaterPressure\n";

        }
        return;
    }




    void processPermeability(ifstream& inFile, ofstream& outFile)
    {
    
        try
        {
            inFile >> perm.viscosity >> perm.density >> perm.perm0;
            outFile << "permeability parameters" << endl
                << "     fluid viscosity = " << perm.viscosity << " kg/(m s)" << endl
                << "     fluid density = " << perm.density << " kg/m3" << endl
                << "     intact rock hydraulic conductivity = " << perm.perm0 << " m/s" << endl << endl;
        }
        catch (std::ifstream::failure e)
        {
            std::cerr << "Exception opening/reading/closing file:in processPermeability\n";

        }
        return;
    }




    void start_calculation(string lineData,ofstream& file7, ofstream& outFile)
    {   
        nc++;
        string  id = "";
        string  tem = "";
        std::stringstream ss(lineData);

        ss >> id >> tem;
        if (tem != "         ") 
            {
               // outFile << tem << endl;  //Sara! maybe we need to add cycle =  to the outfile explaining tem
             
                    mcyc0 += stoi(tem);
               
            }
            //----------------------------
        if (mcyc != 0) return;
            
        s5u.dtol = 0.0001;            
           
        if (nc == 1)
        {
            inputcheck();
            input_tip_check();
            file7 << dispwin.xll << " " << dispwin.xur << " " << dispwin.yll << " " <<
                dispwin.yur << endl;
            // }

             //!-------------------------
            outFile << endl;
            outFile << " " << endl
                << " boundary element data." << endl
                << endl
                << "element   kode   X(center)  Y(center)  length  angle  us or sigma-s   un or sigma-n   material No." << endl;

            float amin = 10000.0;

            // Calculate amin
            for (int m = 0; m < numbe; ++m)
            {
                amin = min(amin, elm_list[m].a);
            }
            s5u.dtol = amin / 1000.0;

            water();    //Sara! for now I commented water later think about

            // Write  elementdata to file
            for (int m = 0; m < numbe; ++m)
            {
                BoundaryElement& elm = elm_list[m];     //elm here is alias 
                float size = 2.0 * elm.a;
                float angle = 180.0 * static_cast<float>(atan2(elm.sinbet, elm.cosbet)) / pi;
                outFile << std::setw(7) << m + 1
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
            outFile << "New fresh fracture properties --- Tensile fractures" << std::endl
                << "               ks = " << std::scientific << std::setprecision(4) << aks_bb << std::endl
                << "               kn = " << std::scientific << std::setprecision(4) << akn_bb << std::endl
                << "              phi = " << std::fixed << std::setprecision(1) << (phi_bb * 180 / pi) << std::endl
                << "        cohesion = " << std::scientific << std::setprecision(4) << coh_bb << std::endl
                << "   dilation angle = " << std::fixed << std::setprecision(1) << (phid_bb * 180 / pi) << std::endl
                << " initial aperture = " << std::scientific << std::setprecision(4) << ap_bb << std::endl
                << "residual aperture = " << std::scientific << std::setprecision(4) << apr_bb << std::endl;

            // Call stiffness_bb for Shear fractures
            stiffness_bb(aks_bb, akn_bb, phi_bb, coh_bb, phid_bb, ap_bb, apr_bb, 2, 1);
            outFile << "New fresh fracture properties --- Shear fractures" << std::endl
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




   void saveData(string lineData, ofstream& outFile)
   {   
       string id ;
       string tem = "  ";
       try
       {
           //inFile.ignore(); // Skip the rest of the line
           std::stringstream ss(lineData);
           ss >> id >> tem;
          // inFile >> id >> tem;
         //  ToLowerCase(id);
           if (tem == "  ")
           {
               outFile << "Incomplete input - Give a file name for saving" << endl;
               // call SendWindowText('Incomplete input - Give a file name for saving'//CHAR(0))
           }
           else
           {
               ofstream file10(tem, ios::binary);
               //save(file10);
               file10.close();
           }
       }
       catch (std::ifstream::failure e)
       {
           std::cerr << "Exception opening/reading/closing file10:in saveData\n";

       }      
        return;
    }




   void endFile(ifstream& inFile, ofstream& outFile)
   {
       string id;
       string message;
       try
       {
           inFile >> id;

           // Error handling     not sure Sara!
           if (id == "ERR")
           {
               int errorCode;
               inFile >> errorCode;
               outFile << "Input data format error, line=" << line << "***" << id << "***" << endl;
               outFile << "Error code: " << errorCode << endl;
               // Read the rest of the line
               getline(inFile, message);
               //call SendWindowText(message//CHAR(0))
               outFile << message << endl;
               return;
           }
           else if (id == "INSF")
           {
               outFile << "Insufficient input data, line=" << line << "***" << id << "***" << endl;
               // Read the rest of the line
               getline(inFile, message);
               //call SendWindowText(message//CHAR(0))
               outFile << message << endl;
               return;
           }
           // End of input file
           //call SendWindowText('End of input file'//CHAR(0))
           lastinput = "endf";
           // outFile << "End of input file" << endl;}

       }
       catch (std::ifstream::failure e)
       {
           std::cerr << "Exception opening/reading/closing file:in endFile func\n";

       }

       return;
   }




   void processExcavation_Cracks(ifstream& inFile, ofstream& outFile)
   {
       exca.ID_Exca = 1;
       try
       {
           inFile >> exca.d_wall >> exca.rand_e;
           outFile << "Excavation induced cracks: " << std::endl;
           outFile << "     Distance into wall = " << exca.d_wall << std::endl;
           outFile << "     Percentage of points with cracks = " << exca.rand_e * 100 << "%" << std::endl;
       }
       catch (std::ifstream::failure e)
       {
           std::cerr << "Exception opening/reading/closing file:in processExcavation_Cracks\n";

       }           

       return;
   }




   void error_Handling(ifstream & inFile)
   {
       std::cerr << "Error reading data from file" << std::endl;
       inFile.close();
       return ;
   }
    



   void symmetryProcessing(string lineData, ofstream& outFile)
   {      
           std::stringstream ss(lineData);
           
           ss >> symm.ksym >> symm.xsym >> symm.ysym;          

           switch (symm.ksym)
           {
           case 0:
               outFile << "No symmetry" << std::endl;
               break;
           case 1:
               outFile << "Symmetry against vertical line x = " << symm.xsym << std::endl;
               break;

           case 2:
               outFile << "Symmetry against horizontal line y = " << symm.ysym << std::endl;
               break;

           case 3:
               outFile << "Symmetry against point x = " << symm.xsym << std::endl <<
                   " y = " << symm.ysym << std::endl;
               break;

           case 4:
               outFile << "Symmetry against vertical line x = " << symm.xsym << std::endl <<
                   " and horizontal line y = " << symm.ysym << std::endl;
               break;
           }
              
       return;
   }



   void frac_energy_processing(string lineData, ofstream& outFile)
   {
       float gic, giic;
       int material = 1;
       std::stringstream ss(lineData);
       try
       {
           ss >> gic >> giic >> material;
           int mm = material;
           outFile << "ROCK MATERIAL = " << material << std::endl;
           outFile << "    Mode I fracture energy = " << std::scientific << std::setprecision(3) << gic << std::endl;
           outFile << "    Mode II fracture energy = " << std::scientific << std::setprecision(3) << giic << std::endl;
           rock1[mm].akic = sqrt(gic * rock1[mm].e / (1 - rock1[mm].pr * rock1[mm].pr));
           rock1[mm].akiic = sqrt(giic * rock1[mm].e / (1 - rock1[mm].pr * rock1[mm].pr));
       }
       catch (std::ifstream::failure e)
       {
           std::cerr << "Exception opening/reading/closing file:frac_energy_process\n";

       }       
       return;
   }



   void frac_toughness_process(string lineData, ofstream& outFile)
   {
       int mm = 0;
       int material = 1;
       float akic0 = 0, akiic0 = 0;
       std::stringstream ss(lineData);
       try
       {
           ss >> akic0 >> akiic0 >> material ;// >> material;    I removed temp vars here and assign directly the rock attrbs. Sara! check again for errors
           outFile << "ROCK MATERIAL = " << material << std::endl;
           outFile << "    Mode I toughness = " << akic0 << std::endl;
           outFile << "    Mode II toughness = " << akiic0 << std::endl;
           mm = material;
           rock1[mm].akic = akic0;
           rock1[mm].akiic = akiic0;
       }
       catch (std::ifstream::failure e)
       {
           std::cerr << "Exception opening/reading/closing file:infrac_toughness_process\n";

       }
       return;
   }




   void elastic_modul_processing(string lineData, ofstream& outFile)
   {
       int material = 1;
       int mm = 0;
       float pr0, e0;
       std::stringstream ss(lineData);

       try
       {
           ss>>pr0 >> e0 >> material;
           mm = material;
           rock1[mm].e = e0;
           rock1[mm].pr = pr0;
           outFile << "ROCK MATERIAL = " << mm << std::endl
            << " Poisson's ratio = " << std::fixed << std::setprecision(2) << rock1[mm].pr << std::endl
            << " Young's modulus = " << std::scientific << std::setprecision(4) << rock1[mm].e << std::endl;
       }
       catch (std::ifstream::failure e)
       {
           std::cerr << "Exception opening/reading/closing file: in elastic_modul_process\n";

       }
       return;
   }




   void processBoundaries(string id, ifstream& inFile, ostream& outFile, ofstream& file9)
   {
       //file9 should be here defined or pass  Sara!

       //---------------straight boundary value increment ------------------
       if (id == "dbou")
       {
           float x1, x2, y1, y2, dss, dnn;
           try
           {
               inFile >> x1 >> x2 >> y1 >> y2 >> dss >> dnn;

               // std::setw(10) << std::fixed << std::setprecision(2)    ? Sara!
               outFile << "boundary value increment, x1,x2,y1,y2,dss,dnn =" << std::endl
                   << x1 << x2 << y1 << y2 << dss << dnn;
               file9 << "('dbou')";
               file9 << x1 << x2 << y1 << y2 << dss << dnn;   //initue.dss  , insitu.dnn
              
           }
           catch (std::ifstream::failure e) 
           {
               std::cerr << "Exception opening/reading/closing file\n";

           }
           insituS.incres = 3;     //Sara!
           
       }

       //---------------arch boundary value increment ------------------
       else if (id == "darc")
       {
           float xcen, ycen, diam1, diam2, ang1, ang2, dss, dnn;
           try
           {

               inFile >> xcen >> ycen >> diam1 >> diam2 >> ang1 >> ang2 >> dss >> dnn;
               outFile << "Arch boundary value increment, xc, xy, d1, d2, ang1, ang2, dss, dnn =\n"
                   << xcen << " " << ycen << " " << diam1 << " " << diam2 << " "
                   << ang1 << " " << ang2 << " " << dss << " " << dnn << std::endl;
               std::ofstream file9("darc.dat");

               file9 << xcen << " " << ycen << " " << diam1 << " " << diam2 << " "
                   << ang1 << " " << ang2 << " " << dss << " " << dnn << std::endl;
           }
           catch (std::ifstream::failure e)
           {
               std::cerr << "Exception opening/reading/closing file: in processBoundaries \n";

           }
          
          
           insituS.incres = 3;
       }
       return;
   }



   void  processMonitoring_points(string lineData,ofstream & outFile)
   {
       ihist++;
       MonitoringPoint mpoint;     //temp object to keep value,later add it to the mon_point list
       std::stringstream ss(lineData);

       try
       {
           ss >> mpoint.xmon >> mpoint.ymon;
           if (ihist > 19)
           {
               outFile << "ERROR: Too many monitoring points (max=19) No monitoring point set for "
                   << mpoint.xmon << " " << mpoint.ymon << std::endl;
               ihist--;
               return;     //label100 int the fortran code
           }
           outFile << "history file = hist" << ihist << ".dat" << std::endl
               << "monitoring point = " << mpoint.xmon << " " << mpoint.ymon << std::endl;
           mpoint_list.push_back(mpoint);      //add new monitoring point

           std::ofstream histPFile("hist" + std::to_string(ihist) + ".dat");           
           histPFile << " -- Monitoring point close to a boundary will be calculated at the boundary and marked as AT BOUNDARY ELEMENT\n"
               << " -- in this case, the boundary normal and shear stress will be provided by sigxx, sigyy, and sigxy\n"
               << " -- if the boundary is vertical or horizontal, the normal stress = sigxx or sigyy, shear stress = sigxy\n"
               << " -- if the boundary is inclined, the normal stress = sigxx = sigyy, shear stress = sigxy\n\n"
               << "  Cycle(iteration)     Time         xp       yp        sigxx         sigyy       sigxy         dx          dy     max creep velocity\n";
       }
       catch (std::ifstream::failure e)
       {
           std::cerr << "Exception opening/reading/closing file: in"
               <<"processMonitoring_points\n";
       }
       return;       

   }



   void  processMonitoring_lines(ifstream& inFile, ostream& outFile)
   {
       lhist++;
       MonitoringLine monlin;
       try
       {
           inFile >> monlin.x1l >> monlin.y1l >> monlin.x2l >> monlin.y2l >> monlin.npl;

           if (lhist > 9)
           {
               outFile << "ERROR: Too many monitoring lines (max=9). No monitoring line set for "
                   << monlin.x1l << " " << monlin.y1l << " " << monlin.x2l << " " << monlin.y2l << std::endl;
               lhist--;
               return;
           }
           outFile << "history file = hist_line" << lhist << ".dat" << std::endl
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
           std::cerr << "Exception opening/reading/closing file: in processMonitoring_lines\n";

       }
       return;

   }




   void processCreepParameters(ifstream& inFile, ostream& outFile)
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
           outFile << "Mode I creep velocity parameter Vmax = " << creep.v1 << "  n = " << creep.nn1 << "\n"
               << "Mode II creep velocity parameter Vmax = " << creep.v2 << "  n = " << creep.nn2 << "\n"
               << "Total creep time (seconds) = " << creep.totalT << "\n"
               << "Minimum time step (seconds) = " << creep.deltaT_min << "\n"
               << "Maximum time step (seconds) = " << creep.deltaT_max << std::endl;
       }
       catch (std::ifstream::failure e)
       {
           std::cerr << "Exception opening/reading/closing file: processCreepParameters\n";
       }
       return;
   }




   void processInsituStress(ifstream& inFile, ostream& outFile)
   {
       float pxx0, pyy0, pxy0;  //!pxx0 etc are used to avoid interference with pxx 
       //(which is increamental)
       
       try
       {
           inFile >> pxx0 >> pyy0 >> pxy0;

           insituS.dsxx = pxx0; // it is assumed that the initial stress field is 0
           insituS.dsxy = pxy0; // the increase in stress field is then pxx0, pyy0, pxy0
           insituS.dsyy = pyy0;
           insituS.incres = 1;

           outFile << "In situ stresses:\n"
               << " xx-component of field stress = " << scientific << pxx0 << "\n"
               << " yy-component of field stress = " << scientific << pyy0 << "\n"
               << " xy-component of field stress = " << scientific << pxy0 << std::endl<<endl;
       }
       catch (std::ifstream::failure e)
       {
           std::cerr << "Exception opening/reading/closing file:processInsituStress\n";

       }
       return;
   }




   void processGravity(ifstream& inFile, ostream& outFile)
   {
       float dens_rock, gy, sh_sv_ratio; // y_surf - surface height
       try
       {
           inFile >> dens_rock >> gy >> sh_sv_ratio >> g.y_surf;
           g.sky = dens_rock * gy;
           g.skx = g.sky * sh_sv_ratio;
           outFile << "Gravitational & tectonic stress:\n"
               << " Rock density = " << dens_rock << "\n"
               << " Gravity acceleration = " << gy << "\n"
               << " Horizontal/vertical stress ratio = " << sh_sv_ratio << "\n"
               << " Ground surface height (y) = " << g.y_surf << std::endl;
       }
       catch (std::ifstream::failure e)
       {
           std::cerr << "Exception opening/reading/closing file:processGravity\n";

       }
       return;

   }




   void process_ContactSurface_Properties(string lineData, ofstream& outFile)
   {
       int jmat = 1;
       std::stringstream ss(lineData);

       if (ss)
       {
           ss >> jmat >> s8[jmat].aks0 >> s8[jmat].akn0 >> s8[jmat].phi0 >> s8[jmat].coh0 >>
               s8[jmat].phid0 >> s8[jmat].apert0 >> s8[jmat].apert_r;
       }
        //not sure about li§ne processing here
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
      

           outFile << "   joint material = " << jmat << std::endl
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




   void processRockStrength_properties(ifstream& inFile, ostream& outFile)
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

           outFile << "Intact rock strength: phi = " << rphi_tem << "\n"
               << "                      coh = " << std::scientific << std::fixed <<rcoh_tem << "\n"
               << "                     sigt = " << std::scientific << std::fixed << rst_tem << std::endl;
           rock1[mm].rphi = rock1[mm].rphi / 180.0 * pi;
       }
       catch (std::ifstream::failure e)
       {
           std::cerr << "Exception opening/reading/closing file:"<<
               " in processRockStrength_propert\n";

       }
       return;
   }



   
   void processRectWindow(ifstream& inFile, ofstream& outFile)
   {
       dispwin.ID_win = 1;
       try
       {
           inFile >> dispwin.xll >> dispwin.xur >> dispwin.yll >> dispwin.yur >>
               dispwin.numx >> dispwin.numy;               
          
           outFile << "Plot window (rect):  x = " << dispwin.xll << " " << dispwin.xur << "\n"
               << "                     y = " << dispwin.yll << " " << dispwin.yur << "\n"
               << "  Grid number: x and y = " << dispwin.numx << " " << dispwin.numy << std::endl<<"\n"<< "\n";
       }
       catch (std::ifstream::failure e)
       {
           std::cerr << "Exception opening/reading/closing file: in processRectWindow\n";

       }
       return;
       
   }



   void  processCircWindow(ifstream& inFile, ofstream& outFile)
   {
       dispwin.ID_win = 2;
       try
       {
           inFile >> dispwin.xc0 >> dispwin.yc0 >> dispwin.radium >> dispwin.numr >> dispwin.numa;
           outFile << "Plot window (circle):  centre = " << dispwin.xc0 << " " << dispwin.yc0 << "\n"
               << "                       radium = " << dispwin.radium << "\n"
               << "       rid number: r and seta = " << dispwin.numr << " " << dispwin.numa << std::endl;

           dispwin.xll = dispwin.xc0 - dispwin.radium;
           dispwin.xur = dispwin.xc0 + dispwin.radium;
           dispwin.yll = dispwin.yc0 - dispwin.radium;
           dispwin.yur = dispwin.yc0 + dispwin.radium;
       }
       catch (std::ifstream::failure e)
       {
           std::cerr << "Exception opening/reading/closing file: in processCircWindow\n";

       }
       return;
   }




   void processFracinitWindow(ifstream& inFile, ofstream& outFile)
   {
       int iwin = 1;
       try
       {
           inFile >> s5u.xmin >> s5u.xmax >> s5u.ymin >> s5u.ymax;

           outFile << "Fracture initiation window:  X = " << s5u.xmin << " " << s5u.xmax << "\n"
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
           std::cerr << "Exception opening/reading/closing file: in processCircWindow\n";

       }
       return;
   }




   void processTunnel(ifstream& inFile, ofstream& outFile)
   {
       float xcen, ycen, diam;
       try
       {
           if (!(inFile >> xcen >> ycen >> diam))
           {
               // Handle input error
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



   void processArch(ifstream& inFile, ofstream& outFile)
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

           outFile << "ARC ---- centre position: x,y  = " << xcen << " " << ycen << "\n"
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
                   gradny, num, material, kode);   //arcss(na)=bvs   arcsn(na)=bvn
               break;

           case 2:
           case 12:
               arc_list[na] = Arch(xcen, ycen, diam / 2.0, ang11, ang12, 0, 0, bvs, bvn, gradsy, 
                   gradny, num, material, kode);    //arcdn(na)=bvn   arcds(na)=bvs
               break;

           case 3:
           case 13:
               arc_list[na] = Arch(xcen, ycen, diam / 2.0, ang11, ang12, 0, bvn, bvs, 0, gradsy,
                   gradny, num, material, kode);    //arcsn(na)=bvn    arcds(na)=bvs
               break;

           case 4:
           case 14:
               arc_list[na] = Arch(xcen, ycen, diam / 2.0, ang11, ang12, bvs, 0, 0, bvn, gradsy,
                   gradny, num, material, kode);      //arcdn(na)=bvn arcss(na)=bvs  
             
           }
           na++;
       }
       catch (std::ifstream::failure e)
       {
           std::cerr << "Exception opening/reading/closing file: in processArch\n";

       }
       return;
   }




   void processEllipticalOpens(ifstream& inFile, ofstream& outFile)
   {
       int num = 0;
       int material = 1;
       float xcen, ycen, diam1, diam2, kode, bvs, bvn, gradsy, gradny;
       try
       {
           inFile >> xcen >> ycen >> diam1 >> diam2 >> kode >> bvs >> bvn >> material >> gradsy >> gradny;
               
           nellipse++;
           Elliptical_opening ep(xcen, ycen, diam1, diam2);  //new Ellip opening is defined
           ellip_list[nellipse] = ep;    //add the new one to list of ellipopening

           std::cout << "ARC ---- centre position: x,y  = " << xcen << " " << ycen << "\n"
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




   void input()
   {
       /* for loading a new model without quitting FRACOd2D.EXE */
       std::ifstream inFile;       // file 1 in code for reading user data

       //this is from the winput function
       std::string string1 = "#MOVIE_V2.21#\n";
    
       if (!file7 ||!file2) {
          //// Debug::WriteLine("Path Wrong!!!! in input function");
           exit(EXIT_FAILURE);
       }
       // Write string1 to the file in unformatted mode
       file7.write(string1.c_str(), string1.size());

       //char id[5], tem[21], message[56];      use string instead of array of char
       //Sara! later think about setting nf,no,na and others to 0 for loading new model

       string id;
       string tem;
       string message;
       string lineData;

       //for high toughness differnet file format used
       //string strPathFile = "C:/C++projects/Fracod2/Examples/HighToughness/Example"+
          // to_string(test_id)+ "/Example"+
          // to_string(test_id)+"_TOUK 1E+20.dat";             // /input.dat";
          // Open an existing file   
          //  
       string strPathFile = filepath+ "/Example" +
           to_string(test_id) + ".dat";
       inFile.open(strPathFile.c_str(),std::ios_base::in);
       if (!inFile.is_open())
       {
         ////  Debug::WriteLine("Path Wrong!!!!" );
           exit(EXIT_FAILURE);
       }   
       try
       {
           while (getline(inFile, lineData))   //!inFile.eof()
           {
               if (lineData.empty())			// skip empty lines:
               {
                   continue;
               }
               line++;
               if (lastinput == "endf")
                   return;              
               id = lineData.substr(0,4);
               ToLowerCase(id);

               line++;
               lastinput = id;

               //---------------job title ------------------
               if (id == "titl")
               {
                   getline(inFile,tem);                    
                   file2 << tem << std::endl<<endl;   
               }

               //--------------- symetry --------------------      
               else if (id == "symm")
               {
                   getline(inFile, tem);
                   symmetryProcessing(tem, file2);
               }

                   //--------------- random fracture initiation --------------------
               else if (id == "rand")
               {
                   string linedata;
                   s15.i_rand = 1, s15.l_rand = 1;
                  
                   try
                   {
                       getline(inFile,linedata);
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

                   //--------------- Set fracture initiation cut-off level --------------------
               else if (id == "setf")
                {
                       inFile >> factors.factor_f;
                       file2 << "Fracture initiation cut-off level is set to be " << factors.factor_f
                           << "% of max FoS" << std::endl;
                }

                //--------------- Set the F/Fc check-up value for fracture propagation --------------------
                 //from an elastic fracture 
                //if F / Fc < given value, no checking for elastic fracture growth
               else if (id == "sete")
                {
                       inFile >> factors.factor_e;
                       file2 << "Elastic fracture propagation checking level is set to be " << 
                           factors.factor_e << "Kc" << std::endl;
                }

                   //--------------- Set the tip merging tolerance distance --------------------
               else if (id == "sett")
                   {
                       inFile >> factors.tolerance;
                       file2 << "Fracture tip merging tolerance distance is set to be " << 
                           factors.tolerance << std::endl;
                   }

                   //--------------- boundary fracture initiation --------------------
               else if (id == "boun")
                   {
                       s15.i_bound = 1;
                       file2 << "Fracture initiation at boundaries is allowed" << std::endl;
                   }

                   //--------------- internal fracture initiation --------------------
                   else if (id == "inte") 
                   {
                       s15.i_intern = 1;
                       file2 << "Fracture initiation in intact rock is allowed" << std::endl;
                   }

                   //--------------- fracture initiation element size --------------------
                   else if (id == "isiz")  
                   {
                       inFile >> s15.a_ini;
                       file2 << "Fracture initiation element size is set to be " <<
                           scientific<< std::setprecision(3)<< s15.a_ini << std::endl;                      
                   }

                   //--------------fracture energy-------------------------------------------

                   else if (id == "toug" ) 
                   {
                        getline(inFile, tem);
                        frac_energy_processing(tem, file2);
                       
                   }
                   //!--------------fracture toughness------------------------------------------

                   else if (id == "touk" )
                   {
                        getline(inFile, tem);
                       frac_toughness_process(tem, file2);
                       
                   }
                   // ---------------elastic modulus---------------------------------------------

                   else if (id == "modu" )
                   {
                        getline(inFile, tem);
                        elastic_modul_processing(tem,file2);
                       
                   }
                   //---------------set iteration number ------------------

                   else if (id == "iter" )
                   {
                       inFile >> n_it;
                       file2 << "Iteration cycle number = " << n_it << std::endl;
                   }
                   //---------------Excavation cracks ------------------

                   else if (id == "exca")
                   {
                       processExcavation_Cracks(inFile, file2);
                       
                   }
                   //-------------- - boundary value increment---------------------------------------------------------- 

                   else if (id == "dbou" || id == "darc")
                   {
                       std::ofstream file9("bound.dat");
                       processBoundaries(id, inFile, file2, file9);
                       
                   }

                    //---------------monitoring points ----------------------------
                   else if (id == "moni" )
                   {
                       getline(inFile, tem);
                       processMonitoring_points(tem, file2);
                   }

                    //---------------monitoring lines ----------------------------------

                    else if (id == "monl") 
                    {
                        processMonitoring_lines(inFile,file2);
                          
                    }

                    //-------------- - set numerical stability parameters k_num and d_max--------------
                    else if (id == "nume" )
                    {
                        try
                        {
                            if (!(inFile >> k_num >> d_max)) {
                                // Handle input error
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

                    // !---------------Definite which material region is concrete lining-------------------

                    else if (id == "mlin" ) 
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
                       // ---------------set creep parameters --------------------------------------------------

                       else if (id == "cree" )
                       {
                           processCreepParameters(inFile, file2);
                          
                       }

                       // ---------------insitu stresses---------------------------------------------------------

                       else if (id == "stre" ) 
                       {
                           processInsituStress(inFile, file2);
                          
                       }
                       // ---------------gravity and tectonic stresses-------------------------------------------------

                       else if (id == "grav" )
                       {
                           processGravity(inFile,file2);
                       }

                       //!---------------insitu stress increment ----------------------------------------------
                       else if (id == "dstr")
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
                       //---------------fracture contact properties------------------------------

                       else if (id == "prop")
                       {
                           getline(inFile, tem);
                           process_ContactSurface_Properties(tem, file2);                                      

                       }
                       //---------------Intact rock strength properties----------------------------

                       else if (id == "rock")
                       {
                           processRockStrength_properties(inFile, file2);                          
                       }                      

                       // ---------------plot window (rectangular) for internal stress and displacement-------

                       else if (id == "swin" )
                       {
                           processRectWindow(inFile, file2);                          
                       }
                       // ---------------plot window (circular) for internal stress and displacement---------

                       else if (id == "rwin")
                       {
                           processCircWindow(inFile, file2);

                       }

                       //!---------------fracture initiation window-------------------------------------

                       else if (id == "iwin" )
                       {
                           processFracinitWindow(inFile, file2);
                           
                       }
                       // Define locations, size, orientations and boundary conditions of
                       //boundary elements

                       // !---------------Tunnel for plotting only------------------------------------------

                       else if (id == "tunn")
                       {
                           processTunnel(inFile, file2);
                       }

                       // !---------------arc--------------------------------------------

                       else if (id == "arch")
                       {
                           processArch(inFile, file2);                          

                       }
                       // !---------------elliptical---------------------------------------------
                       else if (id == "elli")
                       {
                           processEllipticalOpens(inFile, file2);
                       }

                       //-------------- - excavation induced random cracks in an elliptical range--------

                           //------------------ - edge, inner or outer straigt bound.line------------
                       else if (id == "edge")
                           processEdge(inFile, file2, 'S', "Straigt boundary");

                       //--------------------- - gost elements for preventing rigid movement------------
                       else if (id == "gost" )
                           processEdge(inFile, file2, 'G', "Gost element");

                       //!-------------------fractures------------
                       else if (id == "frac")
                           processFracture(inFile, file2);

                       //------------------Linear Interface of Multiregions--------------------
                       else if (id == "lint")
                           processLinearInterface(inFile, file2);

                       //-------------- - Arch Interface of Multiregions--------------------------
                       else if (id == "aint")
                           processArcInterface(inFile, file2);


                       //----------------Water pressure------------------
                       else if (id == "wate")
                           processWaterPressure(inFile, file2);


                       //!----------------permeability ------------------
                       else if (id== "perm") 
                           processPermeability(inFile, file2);


                       else if (id== "cycl")
                           start_calculation(lineData,file7, file2);

                        //!------------------save --------------------------------

                      else if (id== "save")
                               saveData(lineData, file2);

                       //!-------------end file -----------------------------

                      else if (id == "endf")    return;

           }               
               endFile(inFile, file2);               
           }
       catch (std::ifstream::failure e) 
        {
            std::cerr << "Exception opening/reading/closing file in input\n";     
           
       }
       return;
   }

