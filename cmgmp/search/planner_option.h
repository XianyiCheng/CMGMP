#ifndef UTILS_H
#define UTILS_H
    #include "../utilities/utilities.h"
#endif

#define METHOD_QPCC 0
#define METHOD_LCP 1
#define METHOD_ALL 2
#define METHOD_NONE 3

class RRTPlannerOptions{

public:
    int method = METHOD_QPCC;
    int sampleSO3 = true;
    double goal_biased_prob = 0.8;
    int max_samples = 5;
    Vector3d rotation_sample_axis;

    RRTPlannerOptions(){};
    RRTPlannerOptions(const RRTPlannerOptions& opts){
        
        this->method = opts.method;
        this->sampleSO3 = opts.sampleSO3;
        this->goal_biased_prob = opts.goal_biased_prob;
        this->max_samples = opts.max_samples;
        this->rotation_sample_axis = opts.rotation_sample_axis;
    }

};