#include "scene_file.h"
#include <iostream>
#include <fstream>
#include <sstream>

namespace LWS {

    std::string getDirectoryFromPath(std::string str) {
        using namespace std;
        vector<string> parts;
        splitString(str, parts, '/');

        int nParts = parts.size();
        if (nParts == 1) return "./";
        
        string path = "";

        for (int i = 0; i < nParts - 1; i++) {
            path = path + parts[i] + "/";
        }

        return path;
    }

    void processLine(SceneData &data, std::string dir_root, std::vector<std::string> &parts) {
        using namespace std;
        string key = parts[0];

        if (key == "curve") {
            if (parts.size() != 2) {
                std::cerr << "Incorrect arguments to curve" << std::endl;
                exit(1);
            }
            data.curve_filename = dir_root + parts[1];
        }

        else if (key == "repel_curve") {
            if (parts.size() <= 1) {
                data.tpe_alpha = 3;
                data.tpe_beta = 6;
            }
            else if (parts.size() == 3) {
                data.tpe_alpha = stod(parts[1]);
                data.tpe_beta = stod(parts[2]);
            }
            else if (parts.size() == 4) {
                data.tpe_alpha = stod(parts[1]);
                data.tpe_beta = stod(parts[2]);
                data.tpe_weight = stod(parts[3]);
            }
            else {
                std::cerr << "Incorrect arguments to repel_curve" << std::endl;
                exit(1);
            }
        }
    }

    SceneData ParseSceneFile(std::string filename) {
        using namespace std;
        SceneData sceneData;
        string directory = getDirectoryFromPath(filename);
        std::cout << "Directory: " << directory << std::endl;

        ifstream inFile;
        inFile.open(filename);

        if (!inFile) {
            cerr << "Could not open file " << filename << endl;
            exit(1);
        }
    
        string line;
        std::vector<std::string> parts;
        while (inFile >> line) {
            std::cout << line << endl;
            splitString(line, parts, ' ');
            processLine(sceneData, directory, parts);
        }

        inFile.close();
        return sceneData;
    }

}