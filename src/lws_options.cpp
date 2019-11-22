#include "lws_options.h"

namespace LWS {
    int LWSOptions::iterationNum;
    int LWSOptions::frameNum = 0;
    bool LWSOptions::outputFrames = false;
    
    bool LWSOptions::showWindow;
    bool LWSOptions::runLWSFlow;
    bool LWSOptions::runBoundaryFlow;
    bool LWSOptions::runTPE;
    bool LWSOptions::useSobolev = true;
    bool LWSOptions::useMultigrid = true;
    bool LWSOptions::useBarnesHut = true;
    bool LWSOptions::normalizeView = false;
    
    double LWSOptions::tpeAlpha = 3;
    double LWSOptions::tpeBeta = 6;
}
