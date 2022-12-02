#include <differential_gt/cgt.h>


CoopGT::CoopGT(const int& n_dofs, const double& dt): n_dofs_(n_dofs), dt_(dt)
{
  A_.resize(2*n_dofs_,2*n_dofs_);
  B_.resize(2*n_dofs_,n_dofs_);
  C_.resize(n_dofs_,2*n_dofs_);
  X_.resize(2*n_dofs_);
  dX_.resize(2*n_dofs_);
  
  Q11_.resize(2*n_dofs_,2*n_dofs_); Q11_.setZero();
  Q12_.resize(2*n_dofs_,2*n_dofs_); Q12_.setZero();
  Q21_.resize(2*n_dofs_,2*n_dofs_); Q21_.setZero();
  Q22_.resize(2*n_dofs_,2*n_dofs_); Q22_.setZero();
  
  Q1_  .resize(2*n_dofs_,2*n_dofs_);Q1_  .setZero(); 
  Q2_  .resize(2*n_dofs_,2*n_dofs_);Q2_  .setZero(); 
  Q_gt_.resize(2*n_dofs_,2*n_dofs_);Q_gt_.setZero(); 
  
  R1_.resize(n_dofs_,n_dofs_); R1_.setZero();
  R2_.resize(n_dofs_,n_dofs_); R2_.setZero();
                                          
  R_gt_.resize(2*n_dofs_,2*n_dofs_); R_gt_.setZero(); 
  
  reference_.resize(n_dofs_); reference_.setZero();
  
  K_cgt_.resize(2*n_dofs_,2*n_dofs_); K_cgt_.setZero();
  
  state_ok_       = false;
  reference_ok_   = false;
  sys_params_set_ = false;
  cost_params_set_= false;
  gains_set_      = false;
  alpha_set_      = false;
}
void CoopGT::setSysParams( const Eigen::MatrixXd& A,
                        const Eigen::MatrixXd& B)
{
  Eigen::MatrixXd C; C.resize(n_dofs_,2*n_dofs_); C.setZero();
  C.block(0,0,n_dofs_,n_dofs_) = Eigen::MatrixXd::Identity(n_dofs_,n_dofs_);
  setSysParams(A,B,C);
}
void CoopGT::setSysParams( const Eigen::MatrixXd& A,
                        const Eigen::MatrixXd& B,
                        const Eigen::MatrixXd& C)
{
  A_ = A;
  B_ = B;
  C_ = C;
  
  
  ROS_DEBUG_STREAM("A:\n"<<A_);
  ROS_DEBUG_STREAM("B:\n"<<B_);
  ROS_DEBUG_STREAM("C:\n"<<C_);
    
  sys_params_set_ = true;
}

bool CoopGT::getSysParams(Eigen::MatrixXd& A,
                       Eigen::MatrixXd& B,
                       Eigen::MatrixXd& C)
{
  if(!sys_params_set_)
  {
    ROS_ERROR("system params not yet set. return");
    return false;
  }
  
  A = A_;
  B = B_;
  C = C_;
  
  return true;
}

void CoopGT::setCostsParams(const Eigen::MatrixXd& Q11,
                            const Eigen::MatrixXd& Q12,
                            const Eigen::MatrixXd& Q21,
                            const Eigen::MatrixXd& Q22,
                            const Eigen::MatrixXd& R1,
                            const Eigen::MatrixXd& R2)
{
  Q11_ = Q11;
  Q12_ = Q12;
  Q21_ = Q21;
  Q22_ = Q22;
  R1_  = R1;
  R2_  = R2;

  updateGTMatrices();
  
  cost_params_set_ = true;
}
bool CoopGT::getCostMatrices(Eigen::MatrixXd& Q1,
                             Eigen::MatrixXd& Q2,
                             Eigen::MatrixXd& R1,
                             Eigen::MatrixXd& R2)
{
  if (!cost_params_set_)
  {
    ROS_ERROR("Cost params not yet set");
    return false;
  }
  
  Q1 = Q1_;
  Q2 = Q2_;
  R1 = R1_;
  R2 = R2_;
  
  return true;
}


bool CoopGT::setAlpha(const double& alpha)
{
  if(alpha>1 || alpha <0)
  {
    ROS_ERROR_STREAM("weight alpha must be 0 < alpha < 1 . Current value of alpha: "<<alpha);
  }
  alpha_ = alpha;
  alpha_set_ = true;
  return true;
}

bool CoopGT::setCurrentState(const Eigen::VectorXd& x)
{  
  if (x.size() != 2*n_dofs_)
  {
    ROS_ERROR_STREAM("State size is not correct. got: "<< x.size()<<", required: "<< 2*n_dofs_);
    return false;
  }

  X_ = x;
  state_ok_=true;
  return true;
}

bool CoopGT::updateGTMatrices()
{
  if(!alpha_set_)
  {
    ROS_ERROR("parameter alpha not yet set! ");
    return false;
  }
    
  updateGTMatrices(alpha_);
  return true;
}

void CoopGT::updateGTMatrices(const double& alpha )
{
  Q1_ = alpha*Q11_ + (1-alpha)*Q21_;
  Q2_ = alpha*Q12_ + (1-alpha)*Q22_;
  
  Q_gt_ = Q1_ + Q2_;

  R_gt_.topLeftCorner    (R1_.rows(),R1_.cols()) = alpha*R1_;
  R_gt_.bottomRightCorner(R1_.rows(),R1_.cols()) = (1-alpha)*R2_;
}

Eigen::VectorXd CoopGT::getCurrentState(){return X_;};

void CoopGT::computeCooperativeGains(const double& alpha)
{
  setAlpha(alpha);
  if(!updateGTMatrices())
    ROS_ERROR("comething wrong in updating matrices");
  computeCooperativeGains(Q_gt_, R_gt_);
}

void CoopGT::computeCooperativeGains()
{
  computeCooperativeGains(Q_gt_, R_gt_);
}


void CoopGT::computeCooperativeGains(const Eigen::MatrixXd& Q, const Eigen::MatrixXd& R)
{
  Eigen::MatrixXd B_gt; B_gt.resize(B_.rows(),2*B_.cols());
  B_gt << B_,B_;
  Eigen::MatrixXd P_cgt;
  K_cgt_ = solveRiccati(A_, B_gt, Q, R, P_cgt);
  
  ROS_DEBUG_STREAM("K_cgt\n: "<<K_cgt_);

  gains_set_ = true;
}


Eigen::MatrixXd CoopGT::getCooperativeGains()
{
  if(!gains_set_)
    ROS_WARN_STREAM("gains have not yet been computed ! ");
  return K_cgt_;
}

bool CoopGT::setPosReference(const Eigen::VectorXd& ref_1, const Eigen::VectorXd& ref_2)
{
  if(ref_1.size()<n_dofs_ || ref_2.size()<n_dofs_)
  {
    ROS_ERROR_STREAM("reference vectors have wrong length. Expected: "<<n_dofs_<<", got ref_1: "<<ref_1.size()<<" and ref_2: "<<ref_2.size() );
    return false;
  }
  
  Eigen::VectorXd r1; r1.resize(2*n_dofs_); r1.setZero();
  Eigen::VectorXd r2; r2.resize(2*n_dofs_); r2.setZero();
  
  r1.segment(0,n_dofs_) = ref_1;
  r2.segment(0,n_dofs_) = ref_2;
  if(!setReference(r1,r2))
    return false;
    
  return true;
}

bool CoopGT::setReference(const Eigen::VectorXd& ref_1, const Eigen::VectorXd& ref_2)
{
  if(ref_1.size()<2*n_dofs_ || ref_2.size()<2*n_dofs_)
  {
    ROS_ERROR_STREAM("reference vectors have wrong length. Expected: "<<2*n_dofs_<<", got ref_1: "<<ref_1.size()<<" and ref_2: "<<ref_2.size() );
    return false;
  }
    
  reference_ = Q_gt_.inverse() * (Q1_*ref_1 + Q2_*ref_2);
  
  ROS_INFO_STREAM("reference_: "<<reference_.transpose());
  ROS_INFO_STREAM("ref_1: "<<ref_1.transpose());
  ROS_INFO_STREAM("ref_2: "<<ref_2.transpose());
  
  reference_ok_ = true;
  return true;
}

Eigen::VectorXd CoopGT::getReference()
{
  return reference_;
}




Eigen::MatrixXd CoopGT::solveRiccati(const Eigen::MatrixXd &A,
                                  const Eigen::MatrixXd &B,
                                  const Eigen::MatrixXd &Q,
                                  const Eigen::MatrixXd &R, Eigen::MatrixXd &P)
{
  const uint dim_x = A.rows();
  const uint dim_u = B.cols();

  Eigen::MatrixXd Ham = Eigen::MatrixXd::Zero(2 * dim_x, 2 * dim_x);
  Ham << A, -B * R.inverse() * B.transpose(), -Q, -A.transpose();

  Eigen::EigenSolver<Eigen::MatrixXd> Eigs(Ham);

  Eigen::MatrixXcd eigvec = Eigen::MatrixXcd::Zero(2 * dim_x, dim_x);
  int j = 0;
  for (int i = 0; i < 2 * dim_x; ++i) {
    if (Eigs.eigenvalues()[i].real() < 0.) {
      eigvec.col(j) = Eigs.eigenvectors().block(0, i, 2 * dim_x, 1);
      ++j;
    }
  }

  Eigen::MatrixXcd Vs_1, Vs_2;
  Vs_1 = eigvec.block(0, 0, dim_x, dim_x);
  Vs_2 = eigvec.block(dim_x, 0, dim_x, dim_x);
  P = (Vs_2 * Vs_1.inverse()).real();
  
  return R.inverse()*B.transpose()*P;
}

void CoopGT::solveNashEquilibrium( const Eigen::MatrixXd &A,
                                const Eigen::MatrixXd &B1,
                                const Eigen::MatrixXd &B2,
                                const Eigen::MatrixXd &Q1,
                                const Eigen::MatrixXd &Q2,
                                const Eigen::MatrixXd &R1,
                                const Eigen::MatrixXd &R2, 
                                const Eigen::MatrixXd &R12,
                                const Eigen::MatrixXd &R21, 
                                Eigen::MatrixXd &P1,Eigen::MatrixXd &P2)
{
  Eigen::MatrixXd S1  = B1 * R1.inverse() * B1.transpose();
  Eigen::MatrixXd S2  = B2 * R2.inverse() * B2.transpose();
  Eigen::MatrixXd S12 = B1 * R1.inverse() * R21 * R1.inverse() * B1.transpose();
  Eigen::MatrixXd S21 = B2 * R2.inverse() * R12 * R2.inverse()* B2.transpose();

  solveRiccati(A,B1,Q1,R1,P1);
  solveRiccati(A,B2,Q2,R2,P2);
  
  Eigen::MatrixXd P1_prev = P1;
  Eigen::MatrixXd P2_prev = P2;
  double err_1 = 1;
  double err_2 = 1;
  double toll = 0.00001;

  while (err_1>toll && err_2>toll)
  {    
    Eigen::MatrixXd A1 = A - S2*P2;
    Eigen::MatrixXd A2 = A - S1*P1;
    
    Eigen::MatrixXd Q_1 = Q1 + P1*S21*P1;
    solveRiccati(A1,B1,Q_1,R1,P1);
    Eigen::MatrixXd Q_2 = Q2 + P2*S12*P2;
    solveRiccati(A2,B2,Q_2,R2,P2);
  
    err_1 = (P1-P1_prev).norm();
    err_2 = (P2-P2_prev).norm();
    
    P1_prev = P1;
    P2_prev = P2;
  }
  
  return;
}

Eigen::VectorXd CoopGT::computeControlInputs()
{
  if (!state_ok_)
  {
    ROS_WARN_STREAM("State is not updated. computing gains on the last state received: " << X_.transpose());
  }
  if (!reference_ok_)
  {
    ROS_WARN_STREAM("Reference is not updated. computing gains on the last reference received: " << reference_.transpose());
  }
  
  state_ok_     = false;
  reference_ok_ = false;
  
  Eigen::VectorXd control = -K_cgt_ * X_ + K_cgt_ * reference_;
  
  return control;
}

Eigen::VectorXd CoopGT::step(const Eigen::VectorXd& x, const Eigen::VectorXd& ref_1, const Eigen::VectorXd& ref_2)
{
  if (x.size() != 2*n_dofs_)
  {
    ROS_ERROR_STREAM("State size is not correct. got: "<< x.size()<<", required: "<< 2*n_dofs_);
  }
  if (!sys_params_set_)
  {
    ROS_ERROR_STREAM("System parameters not set. use setC2DSysParams or setSysParams to set the parameters before !");
  }
  
  if(ref_1.size()!=ref_2.size())
    ROS_ERROR_STREAM("references are not the same size ! ");
  if(ref_1.size()==n_dofs_ && ref_2.size()==n_dofs_)
    setPosReference(ref_1,ref_2);
  else if(ref_1.size()==2*n_dofs_ && ref_2.size()==2*n_dofs_)
    setReference(ref_1,ref_2);
  else
    ROS_ERROR("references have an incorrect length .");
  
  Eigen::VectorXd u = computeControlInputs();
  
  
  setCurrentState(x);
  dX_ = A_*X_ + B_*u.segment(0,n_dofs_) + B_*u.segment(n_dofs_,n_dofs_);
  X_ = X_ + dX_*dt_;
  
  return X_;
}

Eigen::VectorXd CoopGT::step(const Eigen::VectorXd& ref_1, const Eigen::VectorXd& ref_2)
{
  return step(X_,ref_1,ref_2);
}


