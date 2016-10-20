/* Copyright (C) 2014 - 2016  James Balamuta, Stephane Guerrier, Roberto Molinari
 *
 * This file is part of GMWM R Methods Package
 *
 * The `gmwm` R package is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The `gmwm` R package is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *  
 */

#include <RcppArmadillo.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdexcept>
#include <string.h>
#include <iostream>
#include <cmath>

#include "rtoarmadillo.h"

// Format for sensor information
struct imu_info{
  std::string name;
  int time_type;
  int data_type;
  int header_size;
  double scale_gyro;
  double scale_acc;
};

// Helper function to create the data structure needed
imu_info get_imu_info(std::string imu_type){
  
  // Transform imu_type to capitals
  transform(imu_type.begin(), imu_type.end(), imu_type.begin(), ::toupper);
  
  // Create class definition
  imu_info imu;
  
  // opted for if/else vs. map to save in constructing all imu types. so, O(n) instead of O(log(n)) =(
  if(imu_type == "IMAR"){
    imu.name          = "IMAR";
    imu.time_type     = 8; // double
    imu.data_type     = 4; // int
    imu.header_size   = 0;
    imu.scale_gyro    = 0.10000000*M_PI/180.0/3600.0; // Scale gyro to rad
    imu.scale_acc     = 0.00152588/1000.0;            // scale accel to m/s
  }else if(imu_type == "LN200"){
    imu.name          = "LN200";
    imu.time_type     = 8; // double
    imu.data_type     = 4; // int
    imu.header_size   = 0;
    imu.scale_gyro    = 1.0/2097152.0; // Scale gyro to rad
    imu.scale_acc     = 1.0/16384.0;	 // scale accel to m/s
  }else if(imu_type == "LN200IG"){
    imu.name          = "LN200IG";
    imu.time_type     = 8; // double
    imu.data_type     = 4; // int
    imu.header_size   = 0;
    imu.scale_gyro    = 1.0/524288.0;    // Scale gyro to rad
    imu.scale_acc     = 1.0/16384.0;		 // scale accel to m/s
  }else if(imu_type == "IXSEA"){
    imu.name          = "IXSEA";
    imu.time_type     = 8; // double
    imu.data_type     = 8; // double
    imu.header_size   = 0;
    imu.scale_gyro    = M_PI/180.0/3600.0; // Scale gyro to rad
    imu.scale_acc     = 0.001;		      	 // scale accel to m/s
  }else if(imu_type == "NAVCHIP_FLT"){
    imu.name          = "NAVCHIP_FLT";
    imu.time_type     = 8; // double
    imu.data_type     = 8; // double
    imu.header_size   = 0;
    imu.scale_gyro    = ((1.0/3600.0)/360.0)*2.0*M_PI; // Scale gyro to rad
    imu.scale_acc     = 0.001;			                   // scale accel to m/s
  }else if(imu_type == "NAVCHIP_INT"){
    imu.name          = "NAVCHIP_INT";
    imu.time_type     = 8; // double
    imu.data_type     = 4; // int
    imu.header_size   = 0;
    imu.scale_gyro    = 0.00000625;       // Scale gyro to rad
    imu.scale_acc     = 0.0000390625;			// scale accel to m/s
  }else{
    throw std::runtime_error("The IMU type "+ imu_type + " is not supported");
  }
  
  return imu;
}

//' @title Read an IMU Binary File into R
//' 
//' @description
//' The function will take a file location in addition to the type of sensor it
//' came from and read the data into R.
//' 
//' @param file_path A \code{string} that contains the full file path.
//' @param imu_type A \code{string} that contains a supported IMU type given below.
//' @details
//' Currently supports the following IMUs:
//' \itemize{
//' \item IMAR
//' \item LN200
//' \item LN200IG
//' \item IXSEA
//' \item NAVCHIP_INT
//' \item NAVCHIP_FLT
//' }
//' 
//' We hope to soon be able to support delimited files.
//' @return A matrix with dimensions N x 7, where the columns represent:
//' \describe{
//' \item{Col 0}{Time}
//' \item{Col 1}{Gyro 1}
//' \item{Col 2}{Gyro 2}
//' \item{Col 3}{Gyro 3}
//' \item{Col 4}{Accel 1}
//' \item{Col 5}{Accel 2}
//' \item{Col 6}{Accel 3}
//' }
//' @references
//' Thanks goes to Philipp Clausen of Labo TOPO, EPFL, Switzerland, topo.epfl.ch, Tel:+41(0)21 693 27 55
//' for providing a matlab function that reads in IMUs.
//' The function below is a heavily modified port of MATLAB code into Armadillo/C++. 
//' 
//' @examples
//' \dontrun{
//' read_imu(file_path = "F:/Desktop/short_test_data.imu", imu_type = "IXSEA")
//' }
//' @keywords internal
// [[Rcpp::export]]
arma::field<arma::mat> read_imu(std::string file_path, std::string imu_type) {
  
  // -- File Operations
  
  // Split the file name from the path.
  size_t found = file_path.find_last_of("/\\");
  
  std::string file_loc = file_path.substr(0,found);
  std::string file_name = file_path.substr(found+1);
  
  // Open data file - need "r" to read, and "b" to indicate binary (o.w. fails on Windows)
  FILE *fid = fopen((char*)file_path.c_str(),"rb");
  
  if (fid == NULL){
    throw std::runtime_error("Cannot open the " + file_name + " at " + file_loc);
  }
  
  // -- Check IMU Type requested
  
  // Get data structure
  imu_info imu = get_imu_info(imu_type);
  
  // BitsPerEpoch
  unsigned int BitsPerEpoch = imu.time_type + 6*imu.data_type;
  
  // Set cursor at end of file
  fseek(fid, 0, SEEK_END);
  
  // ftell returns the file position (end location)
  double lsize = ftell(fid); // probably could get away with a float here.
  
  // Count epochs and control it
  double nEpochs = (lsize-imu.header_size)/BitsPerEpoch;
  
  // Is nEpochs an integer? 
  if(trunc(nEpochs) != nEpochs){
    throw std::runtime_error("The file does not have the expected file size. Check the type of IMU or if the file is corrupted.");
  } 
  
  // display info to command window
  Rcpp::Rcout << file_name <<  " contains " << (int)nEpochs << " epochs " << std::endl << "Reading ..." << std::endl;
    
   
  // -- Read time

  // Initialize data matrix
  arma::mat data(nEpochs,7);
  
  // Set cursor at begining of data (e.g. skip header if it exists)
  fseek(fid, imu.header_size, SEEK_SET);
    
  
  // Fill the data matrix
  if(imu.data_type == 8){ // Data is only doubles
    
    double data_buffer[7];
  
    for(unsigned int i = 0; i < nEpochs; i++){
      if(!fread(&data_buffer, 8, 7, fid)) // double
        break; 
      
      for(int j = 0; j < 7; j++){
        data(i,j) = data_buffer[j];
      }
    }
    
  }else{ // Data is a mix of double then 6 ints
    
    double time_buffer[1];
    int32_t data_buffer[6];
    
    for(unsigned int i = 0; i < nEpochs; i++){
      if(!fread(&time_buffer, 8, 1, fid)) // double
        break;
      
      if(!fread(&data_buffer, 4, 6, fid)) // int
        break; 
    
      data(i,0) = time_buffer[0];
      
      for(int j = 1; j <= 6; j++){
        data(i,j) = data_buffer[j-1];
      }
    }
    
  } // end if
  
  // Close the file connection
  fclose(fid);

  // Data Rate
  double fIMU = 1.0/mean_diff(data.col(0));
    
  fIMU = round(fIMU); 
  
  //printf("(data @ %.2f Hz, sGr %f sAc %f) ...", fIMU, imu.scale_gyro, imu.scale_acc); 
    
  // Scale data
  data.cols(1,3) *= fIMU * imu.scale_gyro;  
  data.cols(4,6) *= fIMU * imu.scale_acc; 
  
  arma::vec stats(3);
  stats(0) = fIMU;
  stats(1) = imu.scale_gyro;
  stats(2) = imu.scale_acc;
  
  arma::field<arma::mat> out(2);
  out(0) = data;
  out(1) = stats;
  return out;
}

