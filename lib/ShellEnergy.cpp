//
// Created by Zhen Chen on 10/3/18.
//
#include <fstream>
#include <iostream>
#include "ShellEnergy.h"

void ShellEnergy::load_L_list(std::string file_path)
{
    std::ifstream infile(file_path);
    if(!infile)
        return;
    int num;
    infile >> num;
    double d;
    for(int i=0;i<num;i++)
    {
        infile >> d;
        L_list.push_back(d);
    }
}
