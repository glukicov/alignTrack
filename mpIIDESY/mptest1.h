#ifndef MPTEST1
#define MPTEST1

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <random>

#include <Mille.h>

int track_count; // Number of tracks

int plane_count; // Number of planes

float plane_x_begin; // x-value of first plane
float plane_x_sep; // distance between planes
float plane_thickness; // thickness of planes
float plane_height; // height of planes
float plane_eff; // efficiency of planes
float meas_sigma; // standard deviation of measurements

float displ_sigma;
float drift_sigma;

float[plane_count] plane_pos_devs; // array of plane position deviations (alignment parameter)
float[plane_count] drift_vel_dev; // array of drift velocity deviations (calibration parameter)

float[plane_count] true_plane_effs; // array of efficiencies for each plane
float[plane_count] true_meas_sigmas; // array of measurement sigmas for each plane
float[plane_count] hit_sigmas;


int[plane_count] x_hits;
int[plane_count] i_hits;
int[plane_count] y_hits;
int[plane_count] y_drifts;

void mptest1();

std::vector<float> genlin();


#endif
