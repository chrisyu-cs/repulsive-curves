#include "lws_options.h"

namespace LWS {
    int LWSOptions::iterationNum;
    int LWSOptions::frameNum = 0;
    bool LWSOptions::outputFrames = false;
    
    bool LWSOptions::showWindow;
    bool LWSOptions::runLWSFlow;
    bool LWSOptions::runBoundaryFlow;
    bool LWSOptions::runTPE;
    bool LWSOptions::useSobolev;
    bool LWSOptions::useMultigrid = true;
}
