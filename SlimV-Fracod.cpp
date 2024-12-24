// SlimV-Fracod.cpp : This file contains the 'main' function. Program execution begins and ends there.


//#include <iostream>
//#include<string>
#include<CommonPara.h>
#include<Source.h>

//#include <windows.h>
//#include <commdlg.h>
#include <filesystem> 


//#include<Input.h>
//#include <conio.h>
//#include <process.h>
//#include <stdio.h>





using namespace CommonPara_h::comvar;
using namespace std;





std::wstring openFileDialog() {
    OPENFILENAME ofn;
    wchar_t fileName[MAX_PATH] = L""; 
    ZeroMemory(&ofn, sizeof(ofn));
    ofn.lStructSize = sizeof(OPENFILENAME);
    ofn.hwndOwner = nullptr;
    ofn.lpstrFilter = L"All Files\0*.*\0Text Files\0*.TXT\0"; 
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

   // std::cout << "Please select an input file..." << std::endl;
    std::wstring selectedFile = openFileDialog();

    if (!selectedFile.empty()) {
        std::wcout << "You selected: " << selectedFile << std::endl;
    }
    else {
        std::cout << "No file was selected or an error occurred." << std::endl;
    }
    
    filepath = std::filesystem::path{ selectedFile }.parent_path();
    wstring filename = std::filesystem::path{ selectedFile }.filename();     // "file.txt"


    wcout << L"The program is running\n";   

    elm_list.reserve(500);
    b_elm.reserve(500);
    wstring filename1 = std::filesystem::path{ selectedFile }.stem();     // "file.txt"

    dir = filepath + L"\\" + filename1 + L"_Results";
    
    if (std::filesystem::create_directory(dir) || ERROR_ALREADY_EXISTS == GetLastError()) {
        std::cout << "Directory created successfully.\n";
    }
    else {
        std::cout << "Failed to create directory.\n";
    }
    file2.open(dir + L"/Coutput.dat");

    Central_control(selectedFile);
    file2.close();
          
    std::getline(std::cin, dummy);
}
