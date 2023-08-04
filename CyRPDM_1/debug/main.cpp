#include <iostream>

#include "film.hpp"


int main()
{   
    int dig;
    string res_path;

    Film F;

    int Nx_0;
    double h;
    double L_b, L_e, dL;

    int check = 0;
    string par_file = "./params.txt";
    std::ifstream par_input;
    par_input.open(par_file, std::ios::in);
    if (!par_input) {
        std::cerr << "Unable to open file " << par_file << '\n';
        exit(1);
    }
    for (string line; getline(par_input, line); ) {
        std::stringstream input_line(line);
        string word1, word2;
        input_line >> word1 >> word2;
        if (word1 == "Nx") { F.set_Nx(stoi(word2)); Nx_0 = stoi(word2); ++check; }
        else if (word1 == "n_save") { F.set_n_save(stoi(word2)); ++check; }
        else if (word1 == "res_path") { res_path = word2; ++check; }
        else if (word1 == "phi_ext") { F.set_phi_ext(stod(word2)); ++check; }
        else if (word1 == "L") { F.set_L(stod(word2)); ++check; }
        else if (word1 == "T") { F.set_T(stod(word2)); ++check; }
        else if (word1 == "D_MV") { F.set_D_MV(stod(word2)); ++check; }
        else if (word1 == "D_OV") { F.set_D_OV(stod(word2)); ++check; }
        else if (word1 == "k1_0") { F.set_k1_0(stod(word2)); ++check; }
        else if (word1 == "k2_0") { F.set_k2_0(stod(word2)); ++check; }
        else if (word1 == "k3_0") { F.set_k3_0(stod(word2)); ++check; }
        else if (word1 == "k4_0") { F.set_k4_0(stod(word2)); ++check; }
        else if (word1 == "k5_0") { F.set_k5_0(stod(word2)); ++check; }
        else if (word1 == "an") { F.set_an(stod(word2)); ++check; }
        else if (word1 == "e_f") { F.set_e_f(stod(word2)); ++check; }
        else if (word1 == "e_dl") { F.set_e_dl(stod(word2)); ++check; }
        else if (word1 == "e_cdl") { F.set_e_cdl(stod(word2)); ++check; }
        else if (word1 == "d_dl") { F.set_d_dl(stod(word2)); ++check; }
        else if (word1 == "d_cdl") { F.set_d_cdl(stod(word2)); ++check; }
        else if (word1 == "Temp") { F.set_Temp(stod(word2)); ++check; }
        else if (word1 == "Ff") { F.set_Ff(stod(word2)); ++check; }
        else if (word1 == "e_0") { F.set_e_0(stod(word2)); ++check; }
        else if (word1 == "R") { F.set_R(stod(word2)); ++check; }
        else if (word1 == "A_k") { F.set_A_k(stod(word2)); ++check; }
        else if (word1 == "A_D") { F.set_A_D(stod(word2)); ++check; }
        else if (word1 == "A_L") { F.set_A_L(stod(word2)); ++check; }
        else if (word1 == "dig") { dig = stoi(word2); ++check; }
        else {
            std::cerr << "Error: Argument " << word1 <<  " not recognized.\n";
            exit(-1);
        }       
    }
    par_input.close();

    if (check != 27) {
        std::cerr << "Error: Not all parameters were set.\n";
        exit(check);
    }

    int b = F.init();
    std::cout << b << '\n';

    F.solve();

    for (int i = 1; i < 20; i += 5) {
        F.set_Nx(Nx_0 + i*5);
        F.init();
        F.solve();
    }


    return 0;
}