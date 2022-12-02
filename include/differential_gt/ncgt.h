#pragma once

#include <ros/ros.h>
#include <eigen3/Eigen/Dense>

class NonCoopGT
{
public:

  NonCoopGT(const int& n_dofs, const double& dt);
  
  void setSysParams(const Eigen::MatrixXd& A,
                    const Eigen::MatrixXd& B,
                    const Eigen::MatrixXd& C);
  void setSysParams(const Eigen::MatrixXd& A,
                    const Eigen::MatrixXd& B);
  
  bool getSysParams(Eigen::MatrixXd& A,
                    Eigen::MatrixXd& B,
                    Eigen::MatrixXd& C);

  void setCostsParams(const Eigen::MatrixXd& Q11,
                      const Eigen::MatrixXd& Q12,
                      const Eigen::MatrixXd& Q21,
                      const Eigen::MatrixXd& Q22,
                      const Eigen::MatrixXd& R1,
                      const Eigen::MatrixXd& R2);
  
  void setCostsParams(const Eigen::MatrixXd& Q1,
                      const Eigen::MatrixXd& Q2,
                      const Eigen::MatrixXd& R1,
                      const Eigen::MatrixXd& R2);
  
  bool getCostMatrices(Eigen::MatrixXd& Q1,
                       Eigen::MatrixXd& Q2,
                       Eigen::MatrixXd& R1,
                       Eigen::MatrixXd& R2);
  
  bool setCurrentState(const Eigen::VectorXd& x);
  Eigen::VectorXd getCurrentState();
  
  bool setPosReference(const Eigen::VectorXd& ref_1, const Eigen::VectorXd& ref_2);
  bool setReference(const Eigen::VectorXd& ref_1, const Eigen::VectorXd& ref_2);
  void getReference(Eigen::VectorXd& ref_1, Eigen::VectorXd& ref_2);

  bool computeControlInputs(Eigen::VectorXd& u1, Eigen::VectorXd& u2);
  Eigen::VectorXd step(const Eigen::VectorXd& x, const Eigen::VectorXd& ref_1, const Eigen::VectorXd& ref_2);
  Eigen::VectorXd step(const Eigen::VectorXd& ref_1, const Eigen::VectorXd& ref_2);
  
  void computeNonCooperativeGains();
  void getNonCooperativeGains(Eigen::MatrixXd& K1, Eigen::MatrixXd& K2);

protected:

  
Eigen::MatrixXd solveRiccati(const Eigen::MatrixXd &A,
                             const Eigen::MatrixXd &B,
                             const Eigen::MatrixXd &Q,
                             const Eigen::MatrixXd &R, Eigen::MatrixXd &P);

void solveNashEquilibrium( const Eigen::MatrixXd &A,
                           const Eigen::MatrixXd &B1,
                           const Eigen::MatrixXd &B2,
                           const Eigen::MatrixXd &Q1,
                           const Eigen::MatrixXd &Q2,
                           const Eigen::MatrixXd &R1,
                           const Eigen::MatrixXd &R2, 
                           const Eigen::MatrixXd &R12,
                           const Eigen::MatrixXd &R21, 
                           Eigen::MatrixXd &P1,Eigen::MatrixXd &P2);
  
  Eigen::MatrixXd A_;
  Eigen::MatrixXd B_;
  Eigen::MatrixXd C_;
  Eigen::VectorXd X_;
  Eigen::VectorXd dX_;
    
  Eigen::MatrixXd Q1_;
  Eigen::MatrixXd Q2_;
  Eigen::MatrixXd R1_;
  Eigen::MatrixXd R2_;
  Eigen::MatrixXd R12_;
  Eigen::MatrixXd R21_;
  
  Eigen::MatrixXd K_1_;
  Eigen::MatrixXd K_2_;
  
  Eigen::VectorXd ref_1_;
  Eigen::VectorXd ref_2_;

  bool state_ok_;
  bool reference_ok_;
  bool sys_params_set_;
  bool cost_params_set_;
  bool gains_set_;
  
  int n_dofs_;
  double dt_;
};

