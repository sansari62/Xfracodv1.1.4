#include<CommonPara.h>
#include<WinInterface.h>
#include<Input.h>
#include<Geoplot.h>
#include <chrono>


using namespace WinInterface_h::winvar;
using namespace CommonPara_h::comvar;



// Declaration of the SendWindowMessage function
//extern "C" {
//    void SendWindowMessage(int message);
//    int ShowStatusMessage(const char* string, const char* message);
//    int OutMessage(const char* string);
//}
//
//extern "C" {
//    void GetWindowText(const char* filename, const char* filename1);
//    int WindowMessageIn();
//    int WindowTextIn(const char* string1, const char* string2);
//    int SetFracodMessage(int message);
//    int WindowTextOut(const char* string);
//}
//
//
//void SendWindowMessage(int message) {
//    // Call the corresponding function from the DLL to send the window message
//    int set = SetFracodMessage(message);
//}
//
//void GetWindowText(const char* string1, const char* string2) {
//    // Call the corresponding function from the DLL to get the window text
//    int result = WindowTextIn(string1, string2);
//}
//
//
//
//void GetWindowMessage(int& message) {
//    // Call the corresponding function from the DLL to get the window message
//    message = WindowMessageIn();
//}
//
//
//
//void SendWindowText(const char* string1) {
//    // Call the corresponding function from the DLL to send the window text
//    int set = WindowTextOut(string1);
//}
//
//
//
//void SendWindowText1(const char* string) {
//    // Call the imported function
//    int set = ShowStatusMessage(string, "Message");
//}
//
//void SendWindowDebug(const char* string) {
//    // Call the imported function
//    int set = OutMessage(string);
//}
//
//void SendWindowRunStatus(int i_am_running) {
//    // You may need to call another function here, since the original DLL function SetFracodStatus is commented out
//    // For now, let's assume you want to send a message with ID 50
//    SendWindowMessage(50);
//}
//
//


//void winput(int id_start) {
//    char filename[256];
//    char filename1[256];
//    char string1[12];
//    char File_ID[17];
//    char tem[17];
//    char Form_ID[3];
//    char tem1[4];
//
//    std::strcpy(File_ID, "FRACOD SAVED FILE");
//
//    //int ncount1, ncount2;
//    int mcyc0;
//    int id;
//    int filterindex;
//    int message;
//    int ID_InitialStart = 0;
//    int ID_NotInitialStart = 1;
//
//    int mcyctem;
//   // std::chrono::system_clock::time_point ncount1 = std::chrono::system_clock::now(); //reading system time
//    auto start = std::chrono::steady_clock::now();   //not sure if it should be real sys clock or just clock
//    if (id_start == ID_NotInitialStart)  SendWindowMessage(ID_Pause);       //let C know I am in pause
//
//
//    // Main loop
//    // -- giving C full control of the plot---------------
//    while (true) {
//        GetWindowMessage(message);
//
//        // -------------Save after 5 minutes-------------
//        if (message == ID_NoMessage) {
//            // Define the duration for waiting (5 minutes in this case)
//            constexpr int minutes = 5;
//            auto duration = std::chrono::minutes(minutes);
//            auto end = start + duration;
//            if (id_start == ID_InitialStart) continue;
//
//            //system_clock(ncount2, n1, n2);
//            //std::chrono::system_clock::time_point ncount2 = std::chrono::system_clock::now();
//            auto current_time = std::chrono::steady_clock::now();
//
//
//            if (current_time >= end) {
//                std::ofstream file10("run.sav", std::ios::binary);
//                save();
//                file10.close();
//                SendWindowText("FRACOD current run is saved to run.sav");
//                //start = (auto)10e10;
//
//                using time_point = std::chrono::steady_clock::time_point;
//
//                // Define the start time explicitly
//                auto start = time_point(std::chrono::nanoseconds(static_cast <long long> (10e10)));
//            }
//        }
//
//        //------------- Quit calculation-------------------
//        if (message == ID_Stop) {
//            StopReturn = true;
//            return;
//        }
//
//        //------------- Load file--------------------------
//        if (id_start == ID_NotInitialStart) {
//            if (message == ID_LoadInputFile || message == ID_LoadSavedFile) {
//                SendWindowMessage(ID_NoLoadFile);
//                SendWindowText("Unable to load file while running");
//                continue;
//            }
//        }
//        else {
//            if (message != ID_LoadInputFile && message != ID_LoadSavedFile) {
//                SendWindowMessage(ID_MustLoadFile);
//                SendWindowText("Must load file to start");
//                continue;
//            }
//        }
//
//        if (message == ID_LoadInputFile)
//        {
//            std::strcpy(filename, "input.dat");
//        label200:  GetWindowText(filename, filename1);
//
//            //inquire(FILE = filename, unformatted = Form_ID)
//            for (int i = 0; i < 80; ++i) {
//                if (filename[i] == '.') {
//                    if ((filename[i + 1] == 's' || filename[i + 1] == 'S') &&
//                        (filename[i + 2] == 'a' || filename[i + 2] == 'A') &&
//                        (filename[i + 3] == 'v' || filename[i + 3] == 'V')) {
//                        std::strcpy(Form_ID, "YES");
//                    }
//                    else {
//                        std::strcpy(Form_ID, "NO");
//                    }
//                    break; // Exit loop after processing '.' character
//                }
//            }
//            if (Form_ID == "NO")
//            {
//                std::ifstream file1(filename);
//                std::ofstream file7(filename1, std::ios::binary);
//                std::strcpy(string1, "#MOVIE_V2.21#");
//                file7 << string1;
//                input();
//                CheckRange();     //specify the range of initiation check points
//                return;
//            }
//            else
//            {
//                std::ifstream file10(filename, std::ios::binary);
//                std::ofstream file7(filename1, std::ios::binary);
//                std::strcpy(string1, "#MOVIE_V2.21#");
//                file7 << string1;
//                file10 >> tem;
//                if (tem != "FRACOD SAVED FILE")
//                {
//                    SendWindowText("File is not a FRACOD saved file");
//                    goto label200;     //Sara! check this
//                }
//                restore();
//                file10.close();
//                numbe--;
//                geoplot(ng, ns);
//                numbe++;
//                return;
//            }
//        }
//
//        //--------------------------Cycle----------------------------
//        if (message == ID_CYCLE)
//        {
//            title = win_exchange.w_title;
//            mcyctem = win_exchange.w_cycle;
//            //if(mcyctem == 0) mcyctem = 0;        Sara! i comment this while it's meaningless        //mcyctem = 1000
//            mcyc0 = mcyc + mcyctem; //(mcyc - 1) + mcyctem; mcyc - 1 - orignal cycle
//            id = 0;
//            return;
//        }
//
//        //!--------------------------save file----------------------------
//        if (message == ID_SaveRun)
//        {
//            filterindex = 0;
//            strcpy(filename, "default.sav");
//            GetWindowText(filename, filename1);
//            std::ifstream file10(filename, std::ios::binary);
//            save();
//            file10.close();
//            continue;
//        }
//
//        //--------------------------add farfield stress------------------
//        if (message == ID_ADDFARS)
//        {
//            insituS.incres = 1;
//            insituS.dsxx = win_exchange.w_dsxx;
//            insituS.dsyy = win_exchange.w_dsyy;
//            insituS.dsxy = win_exchange.w_dsxy;
//        }
//
//        //--------------------------add boundary stress------------------
//        if (message == ID_ADDBOUS)
//        {
//            insituS.incres = 2;
//            insituS.dss = win_exchange.w_dss;
//            insituS.dnn = win_exchange.w_dnn;
//        }
//
//        //--------------------------add stress from file------------------
//        // Handle addfils message
//        if (message == ID_ADDFILS)
//        {
//            insituS.incres = 3;
//            GetWindowText(filename, filename1);
//            std::ifstream file99(filename, std::ios::binary);   //to transfer old data into new data file
//            std::ofstream file9("bound.dat", std::ios::binary);
//            if (!file99 || !file9) {
//                SendWindowText("Incorrect file for stress increment data");
//                SendWindowMessage(500); // Send error code 500
//                return;
//            }
//            double x1, x2, y1, y2;    //not sure about these vars
//            while (!file99)
//            {
//                file99 >> x1 >> x2 >> y1 >> y2 >> insituS.dss >> insituS.dnn;
//
//                file9 << "dbou" << std::endl;
//                file9 << x1 << " " << x2 << " " << y1 << " " << y2 << " " << insituS.dss << " " << insituS.dnn << std::endl;
//            } // End of file reached
//
//            file99.close();
//            file9.close();
//        }
//
//        //------------ Handle setcreep message-----set creep parameters------------
//
//        if (message == ID_SETCREEP)
//        {
//            creep.v1 = win_exchange.w_v1;
//            creep.nn1 = win_exchange.w_nn1;
//            creep.v2 = win_exchange.w_v2;
//            creep.nn2 = win_exchange.w_nn2;
//            creep.totalT = win_exchange.w_totalT;
//            creep.deltaT_min = win_exchange.w_deltaT_min;
//            creep.deltaT_max = win_exchange.w_deltaT_max;
//            //std::cout << creep.totalT << " " << creep.deltaT_min << " " << creep.deltaT_max << std::endl;
//            SendWindowText("Creep parameters have been re-defined, use RUN to restart the calculation");
//            creep.ID_creep = 1;
//        }
//
//        // Handle Resize message----------------Resize plot window ------------------
//        if (message == ID_Resize) {
//            // dispwin.title = win_exchange.w_title;    there is no title defined in dispwin Sara!
//            dispwin.xll = win_exchange.w_xll;
//            dispwin.yll = win_exchange.w_yll;
//            dispwin.xur = win_exchange.w_xur;
//            dispwin.yur = win_exchange.w_yur;
//            dispwin.numx = win_exchange.w_numx;
//            dispwin.numy = win_exchange.w_numy;
//            geoplot(ng, ns);
//            SendWindowMessage(ID_UpdatePlot);
//        }
//    }
//
//    return;
//}