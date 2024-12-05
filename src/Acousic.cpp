#include <common.h>
#include <CommonPara.h>
#include <WinInterface.h>


using namespace WinInterface_h::winvar;
using namespace CommonPara_h::comvar;


void compute_AE_based_on_symm(int mx, int my, int mxk, int myk, int AE_indx, int BE_indx,
    float ra, float r, ofstream & file3, ofstream & file7)
{
    /*
    to make the code short and readable I considered
    Ae[i].x = mx*xsymm - mxk* tempx
    Ae[i].y = my*ysymm - myk* tempy this way all conditions are included
    and just paras in the main code are changed
    AE_indx is AE index numbe+i or so and BE_index is i for boundary element index
    ra,r  are generated random number

    */

    float tempx = elm_list[BE_indx].xm + r * elm_list[BE_indx].a * cosf(ra * 2.0 * pi);
    float tempy = elm_list[BE_indx].ym + r * elm_list[BE_indx].a * sinf(ra * 2.0 * pi);

    AE[AE_indx].x = mx * symm.xsym - mxk * tempx;
    AE[AE_indx].y = my * symm.ysym - myk * tempy;
    AE[AE_indx].m = AE[BE_indx].m;

    file3 << AE[AE_indx].x << AE[AE_indx].y << AE[AE_indx].m << endl;      //check the format for writing Sara!

    file7 << AE[AE_indx].x << " " << AE[AE_indx].y << " " << AE[AE_indx].m << std::endl;

}




void AcousticE() 
{


    int mm = 1;                
    float r, ra;             //random numbers
    string filepath1 = filepath + "/CAE" + to_string(test_id) + ".dat";
    ofstream file3(filepath1);

    // based on symm.ksymm value  write in file Sara!
    if (symm.ksym == 0)
    {
        file3 << numbe << std::endl;
        file7<< numbe << endl;
    }
    else if (symm.ksym == 1 || symm.ksym == 2 || symm.ksym == 3)
    {
        file3 << numbe * 2 << std::endl;
        file7 << numbe * 2 << std::endl;
    }
       
    else if (symm.ksym == 4)
    {
        file3 << numbe * 4 << std::endl;
        file7 << numbe * 4 << std::endl;
    }

    for (int i = 0; i < numbe; ++i)
    {
        mm = elm_list[i].mat_no;
        r = static_cast<float>(rand()) / RAND_MAX; // Random number between 0 and 1  Sara!   generating rand check
        ra = static_cast<float>(rand()) / RAND_MAX; // Random number between 0 and 1

        //The location of AE is at the centre of the fracture
        AE[i].x = elm_list[i].xm + r * elm_list[i].a * cosf(ra * 2.0 * pi);
        AE[i].y = elm_list[i].ym + r * elm_list[i].a * sinf(ra * 2.0 * pi);
                
        float AE_m0 = 0.0;

        //if kod# 5non - fracture element is given with 0 AE_m
        if (elm_list[i].kod == 5)
        {
            //Shear slip along fracture element, elastic is considered too
            if (b_elm[i].jstate == 2 || b_elm[i].jstate == 0)
            {
                //= shear modulus * shear dd * area
                AE_m0 = abs(s4.df[i * 2] - s4.df0[i * 2]) * 
                    rock1[mm].e / (2.0 * (1 + rock1[mm].pr))
                    * 2.0 * elm_list[i].a * 1.0;
                //=lame modulus * tensile dd * area
                AE_m0 += abs(s4.df[i * 2+1] - s4.df0[i * 2+1]) * rock1[mm].e *
                    (rock1[mm].pr / ((1 + rock1[mm].pr) *
                    (1 - 2.0 * rock1[mm].pr))) * 2.0 * elm_list[i].a * 1.0;
            }

            //tensile failure along fracture element
            else if (b_elm[i].jstate == 1)
            {
                //=lame modulus * tensile dd * area
                AE_m0 = abs(s4.df[i * 2+1] - s4.df0[i * 2+1]) * rock1[mm].e * (rock1[mm].pr / ((1 + rock1[mm].pr) *
                    (1 - 2.0 * rock1[mm].pr))) * 2.0 * elm_list[i].a * 1.0;
            }
        }
  //label100
        if (AE_m0 == 0.0)
        {
            AE_m0 = 1.0;
        }

        AE[i].m = 2.0 / 3.0 * log(AE_m0) - 6.0;

        //Sara! write in file
        file3 << AE[i].x << AE[i].y << AE[i].m << endl;

        switch (symm.ksym)
        {
        case(1):            
            compute_AE_based_on_symm(2, 0, 1, -1, numbe + i, i,ra,r, file3, file7);
            break;

        case(2):            
            compute_AE_based_on_symm(0, 2, - 1, 1, numbe + i, i, ra, r, file3, file7);
            break;
        case(3):            
            compute_AE_based_on_symm(2, 2, 1, 1, numbe + i, i, ra, r, file3, file7);
            break;
        case(4):

            compute_AE_based_on_symm(2, 0, 1, -1, numbe + i, i, ra, r, file3, file7);
            compute_AE_based_on_symm(0, 2, -1, 1, 2 * numbe + i, i, ra, r, file3, file7);
            compute_AE_based_on_symm(2, 2, 1, 1, 3*numbe + i, i, ra, r, file3, file7);
            break;

        }

        s4.df0[i * 2 ] = s4.df[i * 2];
        s4.df0[i * 2+1] = s4.df[i * 2+1];
    }
    return;
}
