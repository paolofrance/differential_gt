#pragma once

#include <ros/ros.h>
#include <eigen3/Eigen/Dense>

class CoopGT
{
public:

  CoopGT(const int& n_dofs, const double& dt);
  
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
  
  bool getCostMatrices(Eigen::MatrixXd& Q1,
                       Eigen::MatrixXd& Q2,
                       Eigen::MatrixXd& R1,
                       Eigen::MatrixXd& R2);
  
  bool updateGTMatrices();  
  void updateGTMatrices(const double& alpha );
  
  bool setAlpha(const double& alpha);
  bool setCurrentState(const Eigen::VectorXd& x);
  
  bool setPosReference(const Eigen::VectorXd& ref_1, const Eigen::VectorXd& ref_2);
  bool setReference(const Eigen::VectorXd& ref_1, const Eigen::VectorXd& ref_2);
  Eigen::VectorXd getReference();
  
  Eigen::VectorXd computeControlInputs();
  Eigen::VectorXd step(const Eigen::VectorXd& x, const Eigen::VectorXd& ref_1, const Eigen::VectorXd& ref_2);
  Eigen::VectorXd step(const Eigen::VectorXd& ref_1, const Eigen::VectorXd& ref_2);
  Eigen::VectorXd getCurrentState();
  
  void computeCooperativeGains(const double& alpha);
  void computeCooperativeGains();
  void computeCooperativeGains(const Eigen::MatrixXd& Q, const Eigen::MatrixXd& R);
  void computeNonCooperativeGains(const Eigen::MatrixXd& Q, const Eigen::MatrixXd& R, Eigen::MatrixXd& K1_ngt, Eigen::MatrixXd& K2_ngt);
  Eigen::MatrixXd getCooperativeGains();

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
  
  Eigen::MatrixXd theta_;
  Eigen::MatrixXd psi_; 
  
  Eigen::MatrixXd Q1_;
  Eigen::MatrixXd R1_;
  Eigen::MatrixXd Q2_;
  Eigen::MatrixXd R2_;
  
  Eigen::MatrixXd Q11_;
  Eigen::MatrixXd Q12_;
  Eigen::MatrixXd Q21_;
  Eigen::MatrixXd Q22_;
  
  Eigen::MatrixXd Q_gt_;
  Eigen::MatrixXd R_gt_;
  
  Eigen::MatrixXd K_cgt_;
  
  Eigen::VectorXd reference_;

  bool state_ok_;
  bool reference_ok_;
  bool sys_params_set_;
  bool cost_params_set_;
  bool gains_set_;
  bool alpha_set_;
  
  int n_dofs_;
  double dt_;
  double alpha_;
};

