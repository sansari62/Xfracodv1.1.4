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





#include<Input.h>

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
    input();
    lastinput = "";
    input2();

   // Central_control();
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



//
//
//#include <unordered_map>
//#include <functional>
//
//// Sample functions to call based on keywords
//void processKeywordA(ifstream&f) {
//    std::cout << "Function for Keyword A called!" << std::endl;
//}
//
//void processKeywordB(ifstream&f) {
//    std::cout << "Function for Keyword B called!" << std::endl;
//}
//
//void processKeywordC(ifstream& f) {
//    std::cout << "Function for Keyword C called!" << std::endl;
//}
//
//void unknownKeyword(ifstream& f) {
//    std::cout << "Unknown keyword: " << std::endl;
//}
//void endoffile(ifstream& f) {}
//
//int input2() {
//
//
//    std::ifstream inFile;      
//   
//    string id;
//    string tem;
//    string message;
//    string lineData;
//
//
//    string strPathFile = filepath + "/Example" +
//        to_string(test_id) + ".dat";
//    inFile.open(strPathFile.c_str(), std::ios_base::in);
//    if (!inFile.is_open())
//    {
//        ////  Debug::WriteLine("Path Wrong!!!!" );
//        exit(EXIT_FAILURE);
//    }
//
//    // Map keywords to corresponding functions
//    std::unordered_map<std::string, std::function<void(const std::ifstream&)>> funcMap =
//    {
//        {"titl", processKeywordA},
//        {"symm", processKeywordB)},
//        {"rand", processKeywordC)},
//       {"endf",endoffile(inFile)}
//    };
//
//
//    while (std::getline(inFile, lineData)) {
//
//        if (lineData.empty())			// skip empty lines:
//        {
//            continue;
//        }
//        line++;
//        if (lastinput == "endf")
//            return 1;
//        id = lineData.substr(0, 4);
//       // ToLowerCase(id);
//
//        line++;
//        lastinput = id;
//        // Check if the line contains a keyword
//        auto it = funcMap.find(id);
//        if (it != funcMap.end()) {
//            // Call the corresponding function
//            it->second;
//        }
//        else {
//            // Handle unknown keywords
//           // unknownKeyword(lineData);
//        }
//    }
//
//    inFile.close();
//    return 0;
//}