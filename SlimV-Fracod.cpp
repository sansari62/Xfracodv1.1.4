// SlimV-Fracod.cpp : This file contains the 'main' function. Program execution begins and ends there.
#include<stdafx.h>
#include<CommonPara.h>
#include<Source.h>
#include <filesystem> 
#include <algorithm>
#include <cctype>
#include <regex>



#define VERSION "v1.0.2"


using namespace CommonPara_h::comvar;
using namespace std;







void removeCommas(std::string& line) {
    line.erase(std::remove(line.begin(), line.end(), ','), line.end());
}
    



void  fixCommas( std::string& line) {
    std::string result;
    bool lastWasDigit = false;  // Track if last character was a digit

    for (char ch : line) {
        if (ch == ',') {
            if (!result.empty() && lastWasDigit) result += ' ';  // Add space if last was a digit
        }
        else {
            result += ch;
            lastWasDigit = std::isdigit(ch);
        }
    }
    line.swap(result);
    return ;
}




bool startsWithNumber(const std::string& line) {
    std::size_t i = 0;

    // Skip leading spaces
    while (i < line.size() && std::isspace(line[i])) {
        i++;
    }

    // Check if the first non-space character is a minus sign
    if (i < line.size() && line[i] == '-') {
        i++;  // Move to the next character
    }
    return (i < line.size() && std::isdigit(line[i]));
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
        return std::wstring(fileName); 
    }
    else {
        DWORD error = CommDlgExtendedError();
        if (error != 0) {
            std::wcerr << L"Error: Unable to open file dialog. Error Code: " << error << std::endl;
        }
        return L""; 
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
    std::regex tab_only_regex("^\\t+$");

   // std::regex empty_line_regex("^\\s*$");   //("^\\s*$|^\\s*\\*")
    while (std::getline(file, line)) {
        if (startsWithNumber(line)) {
            fixCommas(line);  // Remove commas only from numerical lines
        }
        else if (std::regex_match(line, tab_only_regex))
            continue;
        buffer << line << "\n";  
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
    cout <<
        " XFracod2D v1.0.2 - Fracture Creation and Propagation Code\n" << " Copyright(c) DynaFrax UG LTD.All rights reserved.\n";// <<
    std::wstring selectedFile = openFileDialog();

    if (!selectedFile.empty()) {        

        std::wcout << "\nThe input model: " << selectedFile << std::endl;
    }
    else {
        std::cout << "No file was selected or an error occurred." << std::endl;
    }
    filepath = std::filesystem::path{ selectedFile }.parent_path();
    wstring filename = std::filesystem::path{ selectedFile }.filename();    

    wcout << L"The simulation is running...\n";  
    float dr = 0;
    /*cout << std::setprecision(15)<<"tanf(30.0 / 180.0 * pi) = "<<tanf(30.0 / 180.0 * pi) << ",tanf(0.524)=  " << std::setprecision(15)<<tanf(0.524) << "pi/180.0= " << pi / 180.0 << endl;
    cout << "sinf(85 * pi / 180)= " << std::setprecision(15)<<sinf(85 * pi / 180) << "cosf(5 * pi / 180)= " << std::setprecision(15) << cosf(5 * pi / 180) <<endl;
    cout<<"180. / pi="<< std::setprecision(15) << 180. / pi <<"  cosf(dr)= "<< std::setprecision(15) << cosf(dr)<<"sinf(dr)= "<< sinf(dr)<<endl;
    cout << "180. / pi=" << std::setprecision(15) << 180. / pi << "  cosf(0)= " << std::setprecision(15) << cosf(0) << "sinf(0)= " << std::setprecision(15) << sinf(0) << endl;
    cout << "4.0 * atan(1.0) = "<< std::setprecision(15) << 4.0 * atan(1.0) <<endl;*/

    wstring filename1 = std::filesystem::path{ selectedFile }.stem();     

    dir = filepath + L"\\" + filename1 + L"_Results";
    
    if (std::filesystem::create_directory(dir) || ERROR_ALREADY_EXISTS == GetLastError()) {
        std::cout << "Result's Directory created successfully...\n";
    }
    else {
        std::cout << "Failed to create directory.\n";
    }
    file2.open(dir + L"/Coutput.dat");
    logfile.open(dir + L"/Clog.txt");
    //file57.open(dir + L"/Cpermeability.dat");
    //file9.open(dir+L"/Cbound.dat",  std::ios::in | std::ios::out | std::ios::trunc);   
    if (!logfile.is_open())
    {
        cout << "log file is not opend.\n";       
    }
    logfile<< " The release version: " << VERSION << "\n";

    stress_dir = dir + L"\\" + L"Stress";
    if (std::filesystem::create_directory(stress_dir) || ERROR_ALREADY_EXISTS == GetLastError()) {
    }
    else {
        std::cout << "Failed to create stress directory.\n";
    }
    BE_dir = dir + L"\\" + L"BE";
    if (std::filesystem::create_directory(BE_dir) || ERROR_ALREADY_EXISTS == GetLastError()) {
    }
    else {
        std::cout << "Failed to create BE directory.\n";
    }
    fd_dir = dir + L"\\" + L"Frc_defo";
    if (std::filesystem::create_directory(fd_dir) || ERROR_ALREADY_EXISTS == GetLastError()) {
    }
    else {
        std::cout << "Failed to create frac_defo directory.\n";
    }
    monit_dir  = dir + L"\\" + L"Monitoring";
    if (std::filesystem::create_directory(monit_dir) || ERROR_ALREADY_EXISTS == GetLastError()) {
    }
    else {
        std::cout << "Failed to create monitoring directory.\n";
    }

    file_preprocesing(selectedFile);

    inFile.open(selectedFile.c_str(), std::ios_base::in);

    if (!inFile.is_open())
    {
        MessageBox(nullptr, L"File path is wrong!", L"Error", MB_OK);
        exit(EXIT_FAILURE);
    }

    Central_control();
    file2.close();   
    inFile.close();
    file9.close();
    logfile.close();
          
    return 0;
}
