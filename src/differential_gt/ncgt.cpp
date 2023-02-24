#include <differential_gt/ncgt.h>


NonCoopGT::NonCoopGT(const int& n_dofs, const double& dt): n_dofs_(n_dofs), dt_(dt)
{
  A_.resize(2*n_dofs_,2*n_dofs_);
  B_.resize(2*n_dofs_,n_dofs_);
  C_.resize(n_dofs_,2*n_dofs_);
  X_.resize(2*n_dofs_);
  dX_.resize(2*n_dofs_);
    
  Q1_  .resize(2*n_dofs_,2*n_dofs_);Q1_  .setZero(); 
  Q2_  .resize(2*n_dofs_,2*n_dofs_);Q2_  .setZero(); 
  
  R1_.resize(n_dofs_,n_dofs_); R1_.setZero();
  R2_.resize(n_dofs_,n_dofs_); R2_.setZero();
  
  R12_.resize(n_dofs_,n_dofs_); R12_.setZero();
  R21_.resize(n_dofs_,n_dofs_); R21_.setZero();
  
  ref_1_.resize(2*n_dofs_); ref_1_.setZero();
  ref_2_.resize(2*n_dofs_); ref_2_.setZero();
  
  K_1_.resize(n_dofs_,2*n_dofs_); K_1_.setZero();
  K_2_.resize(n_dofs_,2*n_dofs_); K_2_.setZero();
  
  state_ok_       = false;
  reference_ok_   = false;
  sys_params_set_ = false;
  cost_params_set_= false;
  gains_set_      = false;
  init_P12_       = true;
}
void NonCoopGT::setSysParams( const Eigen::MatrixXd& A,
                              const Eigen::MatrixXd& B)
{
  Eigen::MatrixXd C; C.resize(n_dofs_,2*n_dofs_); C.setZero();
  C.block(0,0,n_dofs_,n_dofs_) = Eigen::MatrixXd::Identity(n_dofs_,n_dofs_);
  setSysParams(A,B,C);
}
void NonCoopGT::setSysParams( const Eigen::MatrixXd& A,
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

bool NonCoopGT::getSysParams(Eigen::MatrixXd& A,
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

void NonCoopGT::setCostsParams(const Eigen::MatrixXd& Q1,
                               const Eigen::MatrixXd& Q2,
                               const Eigen::MatrixXd& R1,
                               const Eigen::MatrixXd& R2,
                               const Eigen::MatrixXd& R12,
                               const Eigen::MatrixXd& R21)
{
  Q1_ = Q1;
  Q2_ = Q2;
  R1_  = R1;
  R2_  = R2;
  R21_ = R21;
  R12_ = R12;
  
  cost_params_set_ = true;
}
void NonCoopGT::setCostsParams(const Eigen::MatrixXd& Q1,
                               const Eigen::MatrixXd& Q2,
                               const Eigen::MatrixXd& R1,
                               const Eigen::MatrixXd& R2)
{
  Eigen::MatrixXd R12 = Eigen::MatrixXd::Zero(n_dofs_, n_dofs_);
  Eigen::MatrixXd R21 = Eigen::MatrixXd::Zero(n_dofs_, n_dofs_);
  setCostsParams(Q1,Q2,R1,R2,R12,R21);
}
bool NonCoopGT::getCostMatrices(Eigen::MatrixXd& Q1,
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

bool NonCoopGT::setCurrentState(const Eigen::VectorXd& x)
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

Eigen::VectorXd NonCoopGT::getCurrentState(){return X_;};

void NonCoopGT::computeNonCooperativeGains()
{
  Eigen::MatrixXd P1,P2;
  solveNashEquilibrium(A_,B_,B_,Q1_,Q2_,R1_,R2_,R12_,R21_,P1,P2);
  K_1_ = R1_.inverse()*B_.transpose()*P1;
  K_2_ = R2_.inverse()*B_.transpose()*P2;  
  gains_set_ = true;
}


void NonCoopGT::getNonCooperativeGains(Eigen::MatrixXd& K1, Eigen::MatrixXd& K2)
{
  if(!gains_set_)
    ROS_WARN_STREAM("gains have not yet been computed ! ");
  K1 = K_1_;
  K2 = K_2_;
}

bool NonCoopGT::setPosReference(const Eigen::VectorXd& ref_1, const Eigen::VectorXd& ref_2)
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

bool NonCoopGT::setReference(const Eigen::VectorXd& ref_1, const Eigen::VectorXd& ref_2)
{
  if(ref_1.size()<2*n_dofs_ || ref_2.size()<2*n_dofs_)
  {
    ROS_ERROR_STREAM("reference vectors have wrong length. Expected: "<<2*n_dofs_<<", got ref_1: "<<ref_1.size()<<" and ref_2: "<<ref_2.size() );
    return false;
  }
  
  ref_1_ = ref_1;
  ref_2_ = ref_2;
  
  reference_ok_ = true;
  return true;
}

void NonCoopGT::getReference(Eigen::VectorXd& ref_1, Eigen::VectorXd& ref_2)
{
  ref_1 = ref_1_;
  ref_2 = ref_2_;
}




Eigen::MatrixXd NonCoopGT::solveRiccati(const Eigen::MatrixXd &A,
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

void NonCoopGT::solveNashEquilibrium( const Eigen::MatrixXd &A,
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
  
  if (init_P12_)
  {
    solveRiccati(A,B1,Q1,R1,P1);
    solveRiccati(A,B2,Q2,R2,P2);
    P1_prev_ = P1;
    P2_prev_ = P2;
    init_P12_ = false;
  }
  else
  {
    P1=P1_prev_;
    P2=P2_prev_;

  }

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
  
    err_1 = (P1-P1_prev_).norm();
    err_2 = (P2-P2_prev_).norm();
    
    P1_prev_ = P1;
    P2_prev_ = P2;
  }
  
  return;
}

Eigen::VectorXd  NonCoopGT::computeControlInputs()
{
  if (!state_ok_)
    ROS_WARN_STREAM("State is not updated. computing gains on the last state received: " << X_.transpose());
  if (!reference_ok_)
    ROS_WARN_STREAM("Reference is not updated. computing gains on the last reference received. ");
  
  state_ok_     = false;
  reference_ok_ = false;
  
  Eigen::VectorXd u1 = -K_1_ * (X_ - ref_1_);
  Eigen::VectorXd u2 = -K_2_ * (X_ - ref_2_);
  
  Eigen::VectorXd control; control.resize(2*n_dofs_);
  control << u1,
             u2;
             
  if(n_dofs_>3)
  {
    control(9) = 0;
    control(10)= 0;
  }

             
  return control;
}

Eigen::VectorXd NonCoopGT::step(const Eigen::VectorXd& x, const Eigen::VectorXd& ref_1, const Eigen::VectorXd& ref_2)
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
  
  Eigen::VectorXd u1,u2,u; u1.resize(n_dofs_);u2.resize(n_dofs_);u.resize(2*n_dofs_);
  u = computeControlInputs();
  u1=u.segment(0,n_dofs_);
  u2=u.segment(n_dofs_,n_dofs_);
  
  setCurrentState(x);
  dX_ = A_*X_ + B_*u1 + B_*u2;
  X_ = X_ + dX_*dt_;
  
  return X_;
}

Eigen::VectorXd NonCoopGT::step(const Eigen::VectorXd& ref_1, const Eigen::VectorXd& ref_2)
{
  return step(X_,ref_1,ref_2);
}


