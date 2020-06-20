#include "lws_options.h"

namespace LWS {
    bool LWSOptions::outputFrames = false;
    
    bool LWSOptions::showWindow;
    bool LWSOptions::runLWSFlow;
    bool LWSOptions::runBoundaryFlow;
    bool LWSOptions::runTPE;
    bool LWSOptions::useSobolev = true;
    bool LWSOptions::useMultigrid = false;
    bool LWSOptions::useBarnesHut = true;
    bool LWSOptions::normalizeView = false;
    
    double LWSOptions::tpeAlpha = 3;
    double LWSOptions::tpeBeta = 6;
}
