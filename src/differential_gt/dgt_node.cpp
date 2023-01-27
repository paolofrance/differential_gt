#include "ros/ros.h"
#include <differential_gt/cgt.h>
#include <differential_gt/ncgt.h>

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
    ref_h(i) = 1;//std::sin(time[i]);
    ref_r(i) = 0.5;//*std::sin(time[i]);
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
  
  NonCoopGT ncgt(n_dofs,dt);
  ncgt.setSysParams(Ac,Bc);
  
  CoopGT cgt(n_dofs,dt);
  cgt.setSysParams(Ac,Bc);
  
  Eigen::VectorXd X=Eigen::VectorXd::Zero(2*n_dofs);
  Eigen::VectorXd dX=Eigen::VectorXd::Zero(2*n_dofs);

  cgt.setCurrentState(X);
  ncgt.setCurrentState(X);
  
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
  Eigen::MatrixXd Rh; Rh.resize(1,1); Rh<< .0005;
  Eigen::MatrixXd Rr; Rr.resize(1,1); Rr<< .0001;
  double alpha = 0.8;
  
  cgt.setAlpha(alpha);
  
  cgt.setCostsParams(Qhh,Qhr,Qrh,Qrr,Rh,Rr);
  
  Eigen::MatrixXd Qh; Qh.resize(2,2);
  Eigen::MatrixXd Qr; Qr.resize(2,2);
  cgt.getCostMatrices(Qh,Qr,Rh,Rr);
  
  ncgt.setCostsParams(Qh,Qr,Rh,Rr);
  
  ROS_INFO_STREAM("Qh: \n"<<Qh);
  ROS_INFO_STREAM("Qr: \n"<<Qr);
  ROS_INFO_STREAM("Rh: \n"<<Rh);
  ROS_INFO_STREAM("Rr: \n"<<Rr);
  
  
  cgt.computeCooperativeGains();
  
  ROS_INFO_STREAM("cgt: \n"<<cgt.getCooperativeGains());
  
  ncgt.computeNonCooperativeGains();
  
  Eigen::MatrixXd Kh,Kr;
  ncgt.getNonCooperativeGains(Kh,Kr);
  ROS_INFO_STREAM("Kh: \n"<<Kh);
  ROS_INFO_STREAM("Kr: \n"<<Kr);
  
  
  std::cin.get();
  
  
  Eigen::VectorXd rh = ref_h.segment(0,1);
  Eigen::VectorXd rr = ref_r.segment(0,1);
  cgt.setPosReference(rh,rr);
  ncgt.setPosReference(rh,rr);  

  for (int i = 0;i<10;i++)  
  {
    
    rh = ref_h.segment(i,1);
    rr = ref_r.segment(i,1);
    
    Eigen::VectorXd cgt_state = cgt.getCurrentState();
    cgt.step(cgt_state ,rh,rr);
    ROS_INFO_STREAM("cgt_state : "<<cgt_state .transpose());
    
    Eigen::VectorXd ncgt_state = ncgt.getCurrentState();
    ncgt.step(ncgt_state ,rh,rr);
    ROS_INFO_STREAM("ncgt_state : "<<ncgt_state .transpose());

//     cgt.setPosReference(rh,rr);
//     cgt.setCurrentState(X);
//     cgt.computeCooperativeGains();
//     Eigen::MatrixXd Kgt = cgt.getCooperativeGains();
//     Eigen::VectorXd control = cgt.computeControlInputs();
//     dX = Ac*X + Bc*control.segment(0,n_dofs) + Bc*control.segment(n_dofs,n_dofs);
//     X = X + dX*dt;
//     ROS_INFO_STREAM("control: " << control.transpose());
//     ROS_INFO_STREAM("state: " << X.transpose());

  }
  
  
  
  
  
  
  return 0;
}





