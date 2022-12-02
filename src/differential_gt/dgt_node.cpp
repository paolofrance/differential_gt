#include "ros/ros.h"
#include <differential_gt/cgt.h>

std::vector<double> range(double min, double max, double dt) {
    std::vector<double> range;
    for(int i=0; i<max/dt; i++) {
        range.push_back(min + i*dt);
    }
    return range;
}


int main(int argc, char **argv)
{
  ros::init(argc, argv, "d_mpc");
  ros::NodeHandle n;
  
  int n_dofs=1;
  double dt = 0.01;
  
  std::vector<double> time = range(0.0,10.0,dt);
  
  
  Eigen::VectorXd ref_h; ref_h.resize(time.size());
  Eigen::VectorXd ref_r; ref_r.resize(time.size());
  
  for (int i = 0;i<time.size();i++)
  {
    ref_h(i) = std::sin(time[i]);
    ref_r(i) = 0.5*std::sin(time[i]);
  }
  
  
  Eigen::MatrixXd Ac;
  Eigen::MatrixXd Bc;
  Eigen::MatrixXd Cc;
  Ac.resize(2*n_dofs,2*n_dofs);
  Bc.resize(2*n_dofs,n_dofs);
  Cc.resize(n_dofs,2*n_dofs);
  
  double m,c,k;
  
  m=10;
  k=0;
  c=25;

  Ac << 0, 1,
        -k/m, -c/m;
        
  Bc << 0,
        1/m;
       
  Cc << 1, 0;
  
  CoopGT cgt(n_dofs,dt);
  cgt.setSysParams(Ac,Bc);
  
  Eigen::VectorXd X=Eigen::VectorXd::Zero(2*n_dofs);
  Eigen::VectorXd dX=Eigen::VectorXd::Zero(2*n_dofs);
  cgt.setCurrentState(X);
  
  Eigen::MatrixXd Qhh; Qhh.resize(2,2); 
  Eigen::MatrixXd Qhr; Qhr.resize(2,2);
  Eigen::MatrixXd Qrr; Qrr.resize(2,2);
  Eigen::MatrixXd Qrh; Qrh.resize(2,2); 
  
  Qhh <<1,0,
       0,0.0001;
  Qhr <<0,0,
       0,0.0001;
       
  Qrr <<1,0,
       0,0.0001;
  Qrh <<0.0001,0,
       0,0;
  Eigen::MatrixXd Rh; Rh.resize(1,1); Rh<< .0001;
  Eigen::MatrixXd Rr; Rr.resize(1,1); Rr<< .0001;
  double alpha = 0.8;
  
  cgt.setAlpha(alpha);
  
  cgt.setCostsParams(Qhh,Qhr,Qrh,Qrr,Rh,Rr);
  
  cgt.computeCooperativeGains();
  Eigen::VectorXd rh = ref_h.segment(0,1);
  Eigen::VectorXd rr = ref_r.segment(0,1);
  cgt.setPosReference(rh,rr);
  Eigen::VectorXd ref = cgt.getReference();
  Eigen::MatrixXd Kgt = cgt.getCooperativeGains();
  ROS_INFO_STREAM("Kgt: "<<Kgt);
  ROS_INFO_STREAM("ref: "<<ref.transpose());
  
  Eigen::MatrixXd control = cgt.computeControlInputs();
  
  ROS_INFO_STREAM("control: " << control);
  

  for (int i = 0;i<10;i++)  
  {
    auto mid = std::chrono::steady_clock::now();
    
    rh = ref_h.segment(i,1);
    rr = ref_r.segment(i,1);
    
//     Eigen::VectorXd state = cgt.getCurrentState();
//     cgt.step(state,rh,rr);
//     ROS_INFO_STREAM("state: "<<state.transpose());

    cgt.setPosReference(rh,rr);
    cgt.setCurrentState(X);
    cgt.computeCooperativeGains();
    Kgt = cgt.getCooperativeGains();
    ROS_INFO_STREAM("Kgt: \n"<<Kgt);
    ref = cgt.getReference();
    ROS_INFO_STREAM("ref: "<<ref.transpose());
  
    Eigen::VectorXd control = cgt.computeControlInputs();
    dX = Ac*X + Bc*control.segment(0,n_dofs) + Bc*control.segment(n_dofs,n_dofs);
    X = X + dX*dt;
    ROS_INFO_STREAM("control: " << control.transpose());
    ROS_INFO_STREAM("state: " << X.transpose());

    auto end = std::chrono::steady_clock::now();
    ROS_INFO_STREAM("time to compute: "<<std::chrono::duration_cast<std::chrono::microseconds>(end- mid).count());
  }
  
  
  
  
  
  
  return 0;
}





