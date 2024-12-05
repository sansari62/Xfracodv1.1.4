// SlimV-Fracod.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include<string>
#include<CommonPara.h>
#include<Source.h>

#include <windows.h>
#include <commdlg.h>



using namespace CommonPara_h::comvar;
using namespace std;





std::wstring openFileDialog() {
    OPENFILENAME ofn;
    wchar_t fileName[MAX_PATH] = L""; // Use wide-character string
    ZeroMemory(&ofn, sizeof(ofn));
    ofn.lStructSize = sizeof(OPENFILENAME);
    ofn.hwndOwner = nullptr;
    ofn.lpstrFilter = L"All Files\0*.*\0Text Files\0*.TXT\0"; // Wide-character filter
    ofn.lpstrFile = fileName;
    ofn.nMaxFile = MAX_PATH;
    ofn.Flags = OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST;

    if (GetOpenFileName(&ofn)) {
        return std::wstring(fileName); // Return as wide string
    }
    else {
        DWORD error = CommDlgExtendedError();
        if (error != 0) {
            std::wcerr << L"Error: Unable to open file dialog. Error Code: " << error << std::endl;
        }
        return L""; // Return empty on failure
    }
}







int main()
{
    std::string dummy;

    /////
   // std::cout << "Please select an input file..." << std::endl;
    /*std::wstring selectedFile = openFileDialog();

    if (!selectedFile.empty()) {
        std::wcout << "You selected: " << selectedFile << std::endl;
    }
    else {
        std::cout << "No file was selected or an error occurred." << std::endl;
    }*/
    ////
    
   
    cout << "The program is running";   

    elm_list.reserve(500);
    b_elm.reserve(500);

    Central_control();
    file2.close();
   

    //int maxThreads = std::thread::hardware_concurrency();  // Get number of hardware threads
//    omp_set_num_threads(maxThreads); // Set OpenMP to use this number of threads
//
//#pragma omp parallel
//    {
//#pragma omp single
//        {
//            Debug::WriteLine(omp_get_num_threads());
//        }
//    }

    
    std::getline(std::cin, dummy);
}
