#ifndef Work_h
#define Work_h

void safety_check();
void third_correction_run();
void run_check();
void calc_bound_stress(int it);
void work0(int mode);
void work1(int mode);
void increment();
void water();

void monitoring_point(float xp, float yp, float& sigxx, float& sigyy, float& sigxy,
    float& ux, float& uy);

void point(float xp, float yp, float& sig1, float& sig2, float& bet, float& sig12,
    float& set, float& disp, float& zet);


#endif // !Work_h

