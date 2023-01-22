#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <cmath>

const double T = 298.15;
const double enthFormH2O = -285.83;
const double enthFormCO2 = -393.51;

double normalBP(const std::unordered_map<std::string, int> &input, std::unordered_map<std::string, std::unordered_map<std::string, double>> data) {
    double Tb = 0;
    for (const auto &x : input) {
        Tb += x.second*data[x.first]["Tb"];
    }
    return Tb+198;
}

double enthForm(const std::unordered_map<std::string, int> &input, std::unordered_map<std::string, std::unordered_map<std::string, double>> data) {
    double Hform = 0;
    for (const auto &x : input) {
        Hform += x.second*data[x.first]["Hform"];
    }
    return Hform+68.29;
}

double critTemp(const std::unordered_map<std::string, int> &input, std::unordered_map<std::string, std::unordered_map<std::string, double>> data) {
    double Tc = 0;
    for (const auto &x : input) {
        Tc += x.second*data[x.first]["Tc"];
    }
    return normalBP(input, data)/(0.584+0.965*Tc-pow(Tc, 2));
}

double critPres(const std::unordered_map<std::string, int> &input, std::unordered_map<std::string, std::unordered_map<std::string, double>> data) {
    double Pc = 0, N = 0;
    for (const auto &x : input) {
        Pc += x.second*data[x.first]["Pc"];
        N += x.second*data[x.first]["N"];
    }
    return pow(0.113+0.0032*N-Pc, -2);
}

double enthVapTnbp(const std::unordered_map<std::string, int> &input, std::unordered_map<std::string, std::unordered_map<std::string, double>> data) {
    return 1.092*8.3144598*normalBP(input, data)*(log(critPres(input, data))-1.013)/(0.93-normalBP(input, data)/critTemp(input, data))/1000;
}

double enthVap(const std::unordered_map<std::string, int> &input, std::unordered_map<std::string, std::unordered_map<std::string, double>> data) {
    return enthVapTnbp(input, data)*(pow((critTemp(input, data)-T)/(critTemp(input, data)-normalBP(input, data)), 0.38));
}

double enthCombAlkanol(const std::unordered_map<std::string, int> &input, std::unordered_map<std::string, std::unordered_map<std::string, double>> data) {
    double Nc = 0;
    for (const auto &x : input) {
        Nc += x.second*data[x.first]["Nc"];
    }
    return -enthForm(input, data)+(Nc*enthFormCO2)+((Nc+1)*enthFormH2O)+enthVap(input, data);
}

int main() {
    std::string line, val, group;
    std::vector<std::string> types;
    std::unordered_map<std::string, int> input;
    std::unordered_map<std::string, std::unordered_map<std::string, double>> data;

    std::fstream finInput("../input.csv", std::ios::in);
    if (finInput.is_open()) {
        while (getline(finInput, line)) {
            std::stringstream str(line);
            int i = 0;
            while(getline(str, val, ',')) {
                if (i == 0) {
                    group = val;
                } else {
                    input[group] = std::stoi(val);
                }
                i++;
            }
        }
    }

    std::fstream finData("../data.csv", std::ios::in);
    if (finData.is_open()) {
        int i = 0;
        while (getline(finData, line)) {
            std::stringstream str(line);
            int j = 0;
            while (getline(str, val, ',')) {
                if (i == 0) {
                    types.push_back(val);
                } else {
                    if (j == 0) {
                        group = val;
                    } else {
                        data[group][types[j]] = std::stod(val);
                    }
                }
                j++;
            }
            i++;
        }
    } else {
        std::cout << "Cannot open data.csv" << "\n";
    }

    std::cout << enthCombAlkanol(input, data) << " kJ mol^-1\n";
    return 0;
}
