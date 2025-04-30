#ifndef Work_h
#define Work_h



void safety_check();
void third_correction_run(int &it);
void run_check(int mode);
void calc_bound_stress(int it);
void work0(int mode);
void work1(int mode);
void increment();
void write_monitor_data(int file_type, int i, int mcyc, int it,const float xmon,
    const float ymon, float sigxx, float sigyy, float sigxy, float uxneg, float uyneg);
//void water();

void monitoring_point(float xp, float yp, float& sigxx, float& sigyy, float& sigxy,
    float& ux, float& uy);




#endif 

