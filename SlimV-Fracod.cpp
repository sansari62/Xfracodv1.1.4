// SlimV-Fracod.cpp : This file contains the 'main' function. Program execution begins and ends there.


#include<CommonPara.h>
#include<Source.h>
#include <filesystem> 
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>


using namespace CommonPara_h::comvar;
using namespace std;







void removeCommas(std::string& line) {
    line.erase(std::remove(line.begin(), line.end(), ','), line.end());
}




bool startsWithNumber(const std::string& line) {
    for (char ch : line) {
        if (std::isdigit(ch)) return true;  // Found a digit at the start (ignoring spaces)
        if (!std::isspace(ch)) return false;  // If non-space & not a digit, return false
    }
    return false;
}





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





void  file_preprocesing(const std::wstring& filename)
{

    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    std::ostringstream buffer;
    std::string line;
    while (std::getline(file, line)) {
        if (startsWithNumber(line)) {
            removeCommas(line);  // Remove commas only from numerical lines
        }
        buffer << line << "\n";  // Store processed line in buffer
    }
    file.close();

    // Overwrite the original file with cleaned data
    std::ofstream outFile(filename);
    if (!outFile) {
        std::cerr << "Error writing to file!" << std::endl;
        return;
    }
    outFile << buffer.str();
    outFile.close();
    return;
}





int main()
{
    std::wstring selectedFile = openFileDialog();

    if (!selectedFile.empty()) {
        std::wcout << "You selected: " << selectedFile << std::endl;
    }
    else {
        std::cout << "No file was selected or an error occurred." << std::endl;
    }
    filepath = std::filesystem::path{ selectedFile }.parent_path();
    wstring filename = std::filesystem::path{ selectedFile }.filename();    


    wcout << L"The simulation is running\n";   

    elm_list.reserve(500);
    b_elm.reserve(500);
    wstring filename1 = std::filesystem::path{ selectedFile }.stem();     

    dir = filepath + L"\\" + filename1 + L"_Results";
    
    if (std::filesystem::create_directory(dir) || ERROR_ALREADY_EXISTS == GetLastError()) {
        std::cout << "Result's Directory created successfully.\n";
    }
    else {
        std::cout << "Failed to create directory.\n";
    }
    file2.open(dir + L"/Coutput.dat");

    file_preprocesing(selectedFile);

    inFile.open(selectedFile.c_str(), std::ios_base::in);

    if (!inFile.is_open())
    {
        MessageBox(nullptr, L"File path is wrong!", L"Error", MB_OK);
        exit(EXIT_FAILURE);
    }

    Central_control();
    file2.close();
          
    return 0;
}
